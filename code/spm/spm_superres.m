function nii_out = spm_superres(img,modality,channel,opt)

% Sanity-check for super-resolution
%--------------------------------------------------------------------------
if ~ismember(modality,{'MRI'})
    error('~ismember(modality, {''MRI''})')
end

% Parameters
%--------------------------------------------------------------------------
C         = numel(img);
vx1       = opt.vx;
proj_mat  = opt.proj_mat;
no_lambda = opt.no_lambda;

% Estimate model parameters (sd,lambda)
%--------------------------------------------------------------------------
sd    = cell(1,C);
lambda = zeros(1,C);
for c=1:C
    N       = numel(img{c});
    sd{c}   = zeros(1,N);
    lambda1 = zeros(1,N);
    for n=1:N
        vx0 = spm_misc('vxsize',img{c}(n).mat);
        if 0%any(round(vx0,3)>1)
            scl = log(max(vx0)); % Ad-hoc fudge-factor
        else
            scl = 1;
        end
        
        sd{c}(n) = estimate_sd(img{c}(n).dat.fname,modality);
        sd{c}(n) = scl*sd{c}(n);
        
        lambda1(n) = estimate_lambda(img{c}(n).dat.fname,modality);
    end
    lambda(c) = mean(lambda1);
end

% For some MR contrasts it is difficult to estimate lambda (e.g. IR), here a simple
% fix is to replace those estimates with the average of the other estimates
nlambda = zeros(size(lambda));
for c=1:C
    if any(strcmp(no_lambda,channel{c}))
        ix         = ~(1:C==c);
        nlambda(c) = mean(lambda(ix));
    end
end

for c=1:C
    if any(strcmp(no_lambda,channel{c}))        
        lambda(c) = nlambda(c);
    end
end

% Get LR images
%--------------------------------------------------------------------------
X   = cell(1,C);
msk = cell(1,C);
for c=1:C % Loop over channels
    N = numel(img{c}); % Loop over within-channel images
    for n=1:N   
        X{c}{n}             = single(img{c}(n).dat(:,:,:)); 
        msk{c}{n}           = spm_misc('msk_modality',X{c}{n},modality); 
        X{c}{n}(~msk{c}{n}) = NaN;        
    end
end

% Compute super-resolved images' orientation matrix and dimensions
%--------------------------------------------------------------------------
is2d = size(X{1}{1},3)==1;

if ~is2d
    [mat1,dm1] = max_bb_orient(img,vx1);
else
    error('is2d');
%     M0  = V{1}{1}.mat;
%     vx0 = spm_misc('vxsize',M0);    
%     dm0 = V{1}{1}.dim;
%     d   = vx0./vx1;
%     D   = diag([d 1]);
%     
%     M1     = M0/D;
%     dm1    = floor(D(1:3,1:3)*dm0')';
%     dm1(3) = 1;
end

% Initialise projection matrices (A)
%--------------------------------------------------------------------------
dat = cell(1,C);
for c=1:C
    dat{c}.mat      = mat1;
    dat{c}.dm       = dm1;
    dat{c}.proj_mat = proj_mat;    
    dat{c}          = init_proj_mat(img{c},dat{c});         
end

% Get a starting estimate using linear interpolation (Y0)
%--------------------------------------------------------------------------
Y0 = get_starting_estimate(X,dat,vx1,dm1,mat1);
  
% Run algorithm
%--------------------------------------------------------------------------
Y = admm_superres(X,Y0,dat,sd,lambda,dm1,vx1,opt.admm);

% Clean FOV
%--------------------------------------------------------------------------
for c=1:C    
    Y{c} = clean_fov(Y{c},dat{c});
end

% Make integer
%--------------------------------------------------------------------------
for c=1:C
    Y{c} = round(Y{c});
    if strcmpi(modality,'MRI')
        Y{c}(Y{c}<0) = 0;
    end
end

% Write results
%--------------------------------------------------------------------------
fnames  = {};
nii_out = nifti;
dtype   = img{1}(1).dat.dtype;
for c=1:C
    N = numel(img{c});
    for n=1:N
        fname           = img{c}(n).dat.fname;
        fnames{end + 1} = fname;    
    end
    
    [pth,nam,ext]   = fileparts(fname);
    nfname          = fullfile(pth,['sr_' nam ext]);
    spm_misc('create_nii',nfname,Y{c},mat1,dtype,'super-resolved');            
    nii_out(c)      = nifti(nfname);        
    fnames{end + 1} = nfname;
end

if opt.verbose
    % Visualise result
    spm_check_registration(char(fnames))
end

% Delete input data
for c=1:C
    N = numel(img{c});
    for n=1:N
        fname = img{c}(n).dat.fname;
        delete(fname);
    end
end
%==========================================================================

%==========================================================================
function Y = admm_superres(X,Y,dat,sd,lambda,dm,vx,pars)

% Set parameters
% -------------------------------------------------------------------------
niter     = pars.niter;
tol       = pars.tol;
rho       = single(pars.rho);
verbose   = pars.verbose;
mu        = pars.mu;
alpha     = pars.alpha;
cgs_tol   = pars.cgs_tol;
cgs_niter = pars.cgs_niter;
est_rho   = pars.est_rho;

C   = numel(lambda); % Number of channels
ndm = 3; % Number of dimensions
if dm(3)==1, 
    ndm = 2; 
end 

nonneg = true(1,C);
for c=1:C
    N = numel(X{c});
    for n=1:N
        X{c}{n}(~isfinite(X{c}{n})) = 0;    
        if sum(sum(sum(X{c}{n}<0)))
           nonneg(c) = false;
        end
    end
end

% Get total number of observed voxels
nm = 0;
for c=1:C
    N = numel(X{c});
    for n=1:N
        nm = nm + numel(X{c}{n});
    end
end

tau = cell(1,C);
for c=1:C
    tau{c} = 1./(sd{c}.^2);
end

% Define Laplace prior (in Fourier space)
% -------------------------------------------------------------------------
prior = laplace_prior(dm,vx);

% Initialize algorithm variables
% -------------------------------------------------------------------------
W = cell(C,ndm);
U = cell(C,ndm);
for c=1:C
    for d=1:ndm
        W{c,d} = zeros(dm,'single');
        U{c,d} = zeros(dm,'single');
    end
end
  
% Calculate gradients (gx,gy,gz)
%-------------------------------------------------------------------------
D = cell(C,ndm);
for c=1:C
    [D{c,1:ndm}] = spm_imbasics('grad',Y{c},vx);
end

% Main loop
% -------------------------------------------------------------------------
ll = -Inf;
for iter=1:niter

    % ---------------------------------------------------------------------
    % Sub-problem U
    if est_rho
        oU = U;
    end
    
    Unorm = single(eps);
    for c=1:C
        for d=1:ndm
            U{c,d} = lambda(c)*D{c,d}+(1/rho)*W{c,d};
            Unorm  = Unorm + U{c,d}.^2;
        end
    end
    Unorm = sqrt(Unorm);
    scale = max(Unorm - 1/rho,0)./Unorm;
    clear Unorm
    
    scale(~isfinite(scale)) = 0;            
    for c=1:C
        for d=1:ndm
            U{c,d} = U{c,d}.*scale;
        end
    end
    clear scale
    
    % Sub-problem Y
    % ---------------------------------------------------------------------
    for c=1:C
        % Compute divergence
        if ndm==3
            DtU = lambda(c)*spm_imbasics('dive',U{c,:},vx);
            DtW = lambda(c)*spm_imbasics('dive',W{c,:},vx);
        else
            DtU = lambda(c)*spm_imbasics('dive',U{c,:},[],vx);
            DtW = lambda(c)*spm_imbasics('dive',W{c,:},[],vx);
        end

        LHS = @(Y) AtA(Y,@(Y) lambda(c)^2*prior(Y),tau{c}/rho,1,dat{c});
        RHS = DtU - (1/rho)*DtW + (1/rho)*At(X{c},tau{c},dat{c});
        clear DtU DtW
        
        % Solve least-squares problem
        [Y{c},cg_j,cg_d,cg_tol] = cg_image_solver(LHS,RHS,Y{c},cgs_niter,cgs_tol);
        clear RHS LHS
        
        if verbose
            fprintf('%2d | %8.7f %8.7f %8.7f\n',iter,cg_j,cg_d,cg_tol);        
        end
            
        if nonneg(c)
            % Non-negativity 'constraint'
            Y{c}(Y{c}<0) = 0;
        end
    end

    % Sub-problem W
    % ---------------------------------------------------------------------
    for c=1:C
        [D{c,1:ndm}] = spm_imbasics('grad',Y{c},vx);
        for d=1:ndm
            W{c,d} = W{c,d} + rho*(lambda(c)*D{c,d} - U{c,d});
        end
    end

    % Check convergence
    % ---------------------------------------------------------------------
    ll1 = log_likelihood(X,Y,D,dat,tau,lambda);    
    ll  = [ll,ll1];    
    
    d1 = abs((ll(end - 1)*(1 + 10*eps) - ll(end))/ll(end));
    d2 = ll(end) - ll(end - 1);
    
    if verbose
        fprintf('%2d | %8.7f %8.7f %8.7f %8.7f %8.7f\n',iter,ll1,d1,d2,1e-4*nm,rho);
        
        show_img(X,Y,lambda,sd,'MRI','(SPM) Super-resolution');
    end    
    
%     if iter>=20 && d1<tol
%         % Finished
%         break; 
%     end        
   
    if est_rho
        % Compute primal and dual residuals
        %------------------------------------------------------------------
        primal_res = 0;
        dual_res   = 0;
        for c=1:C
            tmp = cell(1,ndm);
            for d=1:ndm
                primal_res = primal_res + sum(sum(sum((D{c,d} + U{c,d}).^2)));
                tmp{d}     = U{c,d} - oU{c,d};
            end
            if ndm==3
                tmp = rho*spm_imbasics('dive',tmp{:},vx);
            else
                tmp = rho*spm_imbasics('dive',tmp{:},[],vx);
            end
            dual_res = dual_res + sum(sum(sum(tmp.^2)));
        end
        clear tmp

        primal_res = sqrt(primal_res);
        dual_res   = sqrt(dual_res);

        % Compute varying penalty parameter
        %------------------------------------------------------------------
        scale = 1;
        if primal_res>mu*dual_res
            rho   = alpha*rho;
            scale = 1/alpha;
        elseif dual_res>mu*primal_res
            rho   = rho/alpha;
            scale = alpha;
        end
        if scale~=1
            for c=1:C
                for d=1:ndm
                    W{c,d} = W{c,d}*scale;
                end
            end
        end  
    end
end
%==========================================================================

%==========================================================================
function ll = log_likelihood(X,Y,D,dat,tau,lambda)
C   = numel(X);
ndm = numel(size(X{1}));
ll1 = 0;
llr = 0;
for c=1:C
    N  = numel(X{c});
    AY = A(Y{c},dat{c});
    for n=1:N
        I     = numel(X{c}{n});
        sd    = sqrt(1/tau{c}(n));
        diff1 = (AY{n} - X{c}{n});
        diff1 = diff1(:)'*diff1(:);
        ll1   = ll1 + 0.5*(2*I*log(sd) + I*log(2*pi) + 1/sd^2*diff1);
    end
    
    for d=1:ndm
        llr = llr + (lambda(c)^2)*(D{c,d}.^2);
    end
end
ll = -ll1 - sum(sum(sum(sqrt(llr))));
%==========================================================================

%==========================================================================
function X = A(Y,dat)
if strcmp(dat.proj_mat,'sinc') 
    Y(~isfinite(Y)) = 0;
    Y               = fftn(Y);    
end

X = cell(1,dat.N);
for n=1:dat.N
    T = dat.mat\dat.A(n).mat;
    y = make_deformation(T,dat.A(n).dm);
    
    if strcmp(dat.proj_mat,'sinc')       
        tmp  = real(ifftn(Y.*dat.A(n).S));                      
    elseif strcmp(dat.proj_mat,'smo')
        tmp = zeros(size(Y),'single');
        spm_smooth(Y,tmp,dat.A(n).S);
    end
    
    X{n} = samp0(tmp,y); 
    
    vx1  = spm_misc('vxsize',dat.mat);
    vx0  = spm_misc('vxsize',dat.A(n).mat);
    scl  = prod(vx1./vx0);
    X{n} = scl*X{n};
end
%==========================================================================  

%==========================================================================
function Y = At(X,tau,dat)  
Y = single(0);   
for n=1:dat.N      
    T = dat.mat\dat.A(n).mat;
    y = make_deformation(T,dat.A(n).dm);
    
    tmp = spm_diffeo('push',X{n},y,dat.dm);     
    tmp(~isfinite(tmp)) = 0;
    
    if strcmp(dat.proj_mat,'sinc')               
        tmp = real(ifftn(fftn(tmp).*dat.A(n).S));                    
    elseif strcmp(dat.proj_mat,'smo')
        spm_smooth(tmp,tmp,dat.A(n).S); 
    end 
     
    Y = Y + tau(n).*tmp;           
end
%==========================================================================

%==========================================================================
function Y1 = AtA(Y,prior,tau,lambda,dat)  
Y1 = single(0);
if strcmp(dat.proj_mat,'sinc') 
    F               = Y;        
    F(~isfinite(F)) = 0;
    F               = fftn(F);    
end
for n=1:dat.N
    T = dat.mat\dat.A(n).mat;
    y = make_deformation(T,dat.A(n).dm);
        
    if strcmp(dat.proj_mat,'sinc')      
        tmp = real(ifftn(F.*dat.A(n).S));    
    elseif strcmp(dat.proj_mat,'smo') 
        tmp = zeros(size(Y),'single');
        spm_smooth(Y,tmp,dat.A(n).S); 
    end
    
    tmp = samp0(tmp,y); 
    vx1 = spm_misc('vxsize',dat.mat);
    vx0 = spm_misc('vxsize',dat.A(n).mat);
    scl = prod(vx1./vx0);
    tmp = scl*tmp;
    
    tmp                 = spm_diffeo('push',tmp,y,dat.dm);          
    tmp(~isfinite(tmp)) = 0;
        
    if strcmp(dat.proj_mat,'sinc')                                          
        tmp = real(ifftn(fftn(tmp).*dat.A(n).S));    
    elseif strcmp(dat.proj_mat,'smo')
        spm_smooth(tmp,tmp,dat.A(n).S);    
    end    
    
    Y1 = Y1 + tau(n).*tmp;
end

Y1 = Y1 + lambda*prior(Y);
%==========================================================================

%==========================================================================
function [Y,j,d,tol] = cg_image_solver(LHS,RHS,Y,niter,tol,verbose)
if nargin<3, Y       = zeros(size(RHS),'single'); end
if nargin<4, niter   = 32; end
if nargin<5, tol     = 1e-3; end
if nargin<6, verbose = false; end

% Initilisation  
%--------------------------------------------------------------------------
normRHS = sqrt(sum(RHS(:).*RHS(:))); % Norm of RHS
R       = RHS - LHS(Y);              % Residual RHS - LHS(x)
clear RHS

normR = sum(R(:).*R(:));           % R'R
P     = R;                         % Initial conjugate directions P
beta  = 0;                         % Initial search direction for new P

if verbose
    fprintf('%g %g\n', normRHS, sqrt(normR));
end

% Run algorithm
%--------------------------------------------------------------------------
j = 1;
while sqrt(normR) > tol*normRHS,
    % Calculate conjugate directions P which defines the direction of descent
    %----------------------------------------------------------------------
    P = R + beta*P;

    % Finds the step size of the conj. gradient descent
    %----------------------------------------------------------------------
    AtAP  = LHS(P);
    alpha = normR / sum(P(:).*AtAP(:));
    
    % Perform conj. gradient descent, obtaining updated X and R, using the calculated
    % P and alpha
    %----------------------------------------------------------------------
    Y = Y + alpha *P; 
    R = R - alpha*AtAP;
    clear alpha AtAP
    
    % Finds the step size for updating P
    %----------------------------------------------------------------------    
    normR = sum(R(:).*R(:));       
    if verbose
        fprintf('%g\n', sqrt(normR));
    end
    
    RtRp = normR;    
    beta = normR / RtRp;            
    clear RtRp
    
    % Check if converged
    %----------------------------------------------------------------------
    if j>=niter, 
        % Finished
        break; 
    end;                   

    j = j + 1; 
end

d   = sqrt(normR);
tol = tol*normRHS;
%==========================================================================

%==========================================================================
function [mat,dim] = max_bb_orient(img,vx)
C = numel(img);

% Concatenate all NIfTIs
Nii = nifti;
cnt = 1;
for c=1:C
    N = numel(img{c});
    for n=1:N
       Nii(cnt) = img{c}(n);
       cnt      = cnt + 1;
    end
end

% Calculate bounding-box orientation and dimensions
mn = [ Inf  Inf  Inf]';
mx = [-Inf -Inf -Inf]';
for i=1:numel(Nii),
    d = size(Nii(i).dat);
    if numel(d)==2, d(3) = 1; end
        
    t = uint8(0:7);
    c = diag(d+1)*double([bitshift(bitand(t,bitshift(uint8(1),1-1)),1-1)
                          bitshift(bitand(t,bitshift(uint8(1),2-1)),1-2)
                          bitshift(bitand(t,bitshift(uint8(1),3-1)),1-3)]);
    c = bsxfun(@plus,Nii(i).mat(1:3,1:3)*c,Nii(i).mat(1:3,4));
    mx = max(mx,max(c,[],2));
    mn = min(mn,min(c,[],2));
end
mat = spm_matrix(mn-1)*diag([vx 1])*spm_matrix(-[1 1 1]);
dim = ceil((mat\[mx'+1 1]')');
dim = dim(1:3);
%==========================================================================

%==========================================================================
function dat = init_proj_mat(img,dat)
N  = numel(img);
Mj = dat.mat;
dj = dat.dm;
    
dat.N = N;
for n=1:N
    dat.A(n).mat = img(n).mat;    
    dat.A(n).dm  = size(img(n).dat(:,:,:));
    
    if strcmp(dat.proj_mat,'sinc')
        Mi         = dat.A(n).mat;
        M          = Mj\Mi;
        R          = (M(1:3,1:3)/diag(sqrt(sum(M(1:3,1:3).^2))))';
        dat.A(n).S = blur_fun(dj,R,[sqrt(sum(M(1:3,1:3).^2))]);
    elseif strcmp(dat.proj_mat,'smo')
        vxj        = abs(det(dat.mat(1:3,1:3)))^(1/3);
        vxi        = sqrt(sum(dat.A(n).mat(1:3,1:3).^2));        
        dat.A(n).S = sqrt(max(vxi.^2 - vxj.^2,0));
    end
end
%==========================================================================

%==========================================================================
function f = blur_fun(d,M,n)
if nargin<1, d = [64 64]; end
if nargin<2, M = eye(numel(d)); end
if nargin<3, n = ones(1,numel(d)); end
if any(size(M)~=numel(d)) || numel(n)~=numel(d), error('Incompatible dimensions.'); end

r    = cell(1,numel(d));
X    = cell(1,numel(d));
for i=1:numel(d), r{i} = single([0:ceil(d(i)/2-1) -floor(d(i)/2):-1]'*pi/d(i)); end
[X{:}] = ndgrid(r{:});

Y  = cell(size(X));
for i=1:numel(d),
    Y{i} = single(0);
    for j=1:numel(d), Y{i} = Y{i} + M(i,j)*X{j}; end
end
clear X

f  = single(0);
for i=1:numel(d), f = f + Y{i}.^2; end;
f  = ((cos(min(f,pi^2/4)*4/pi)+1)/2);  % Some sort of window function

for i=1:numel(d),
    tmp = sin((n(i))*Y{i})./(Y{i}.*cos(Y{i}/pi^(1/2)));
    tmp(~isfinite(tmp)) = n(i);
    f = f.*tmp;
end
%==========================================================================

%==========================================================================
function y = make_deformation(M,dm)
[x0,y0,z0] = ndgrid(single(1:dm(1)),...
                    single(1:dm(2)),...
                    single(1:dm(3)));
y          = cat(4,M(1,1)*x0 + M(1,2)*y0 + M(1,3)*z0 + M(1,4), ...
                   M(2,1)*x0 + M(2,2)*y0 + M(2,3)*z0 + M(2,4), ...
                   M(3,1)*x0 + M(3,2)*y0 + M(3,3)*z0 + M(3,4));
%==========================================================================   

%==========================================================================
function prior = laplace_prior(dm,vx)
L = zeros(dm);
if dm(1)>=2
    tmp        = 1/(vx(1)^2);
    L(  1,1,1) = L(  1,1,1) + tmp*2;
    L(  2,1,1) = L(  2,1,1) - tmp;
    L(end,1,1) = L(end,1,1) - tmp;
end
if dm(2)>=2
    tmp        = 1/(vx(2)^2);
    L(1,  1,1) = L(1,  1,1) + tmp*2;
    L(1,  2,1) = L(1,  2,1) - tmp;
    L(1,end,1) = L(1,end,1) - tmp;
end
if dm(3)>=2
    tmp        = 1/(vx(3)^2);
    L(1,1,  1) = L(1,1,  1) + tmp*2;
    L(1,1,  2) = L(1,1,  2) - tmp;
    L(1,1,end) = L(1,1,end) - tmp;
end
L     = single(fftn(L));
prior = @(Y) real(ifftn(fftn(Y).*L));
%==========================================================================

%==========================================================================
function Y0 = get_starting_estimate(X,dat,vx1,dm1,mat1)
C  = numel(X);
Y0 = cell(1,C);
for c=1:C
    Y0{c} = zeros(dm1,'single');
    for n=1:dat{c}.N
        vx0 = spm_misc('vxsize',dat{c}.A(n).mat);
        
        T = dat{c}.A(n).mat\mat1;
        
        
        [x0,y0,z0] = ndgrid(1:dm1(1),...
                            1:dm1(2),...
                            1:dm1(3));

        x1 = T(1,1)*x0 + T(1,2)*y0 + T(1,3)*z0 + T(1,4);
        y1 = T(2,1)*x0 + T(2,2)*y0 + T(2,3)*z0 + T(2,4);
        z1 = T(3,1)*x0 + T(3,2)*y0 + T(3,3)*z0 + T(3,4);
        
        C   = spm_bsplinc(X{c}{n},[4 4 4 0 0 0]);
        img = spm_bsplins(C,x1,y1,z1,[4 4 4 0 0 0]);
        img(~isfinite(img) | img<0) = 0;
                
%         scl = prod(vx1./vx0);
%         img = scl*img;
    
        Y0{c} = Y0{c} + single(img);  
    end
    Y0{c} =  Y0{c}/dat{c}.N;
end
%==========================================================================

%==========================================================================
function Y = clean_fov(Y,dat)
mat1 = dat.mat;      
dm1  = dat.dm; 
msk  = cell(dat.N,3);
for n=1:dat.N
    mat0 = dat.A(n).mat;  
    dm0  = dat.A(n).dm;

    % Get the mapping from Mref to Mmod
    T = mat0\mat1;

    % Use ndgrid to give an array of voxel indices
    [x0,y0,z0] = ndgrid(single(1:dm1(1)),...
                        single(1:dm1(2)),...
                        single(1:dm1(3)));

    % Transform these indices to the indices that they point to in the reference image
    D = cat(4,T(1,1)*x0 + T(1,2)*y0 + T(1,3)*z0 + T(1,4), ...
              T(2,1)*x0 + T(2,2)*y0 + T(2,3)*z0 + T(2,4), ...
              T(3,1)*x0 + T(3,2)*y0 + T(3,3)*z0 + T(3,4));

    % Mask according to whether these are < 1 or > than the dimensions of the reference image.        
    for i=1:3
        msk{n,i} = D(:,:,:,i) >= 1 & D(:,:,:,i) <= dm0(i);
    end
end  

% Generate cleaned up image
for n=1:dat.N
    for i=1:3
        Y = msk{n,i}.*Y;
    end
end
%==========================================================================