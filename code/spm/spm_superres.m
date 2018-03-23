function obj = spm_superres(obj)

% Parameters
%--------------------------------------------------------------------------
vx1       = obj.preproc.superres.vx;
V         = obj.scans;
N         = numel(V);
pars_sr   = obj.preproc.superres;
pars_admm = pars_sr.admm;
proj_mat  = pars_sr.proj_mat;
trunc_ct  = pars_sr.trunc_ct;

% Initialise
%--------------------------------------------------------------------------

% Get LR images
X   = {};
msk = {};
for n=1:N
    I = numel(V{n});
    for i=1:I
        Nii                   = nifti(V{n}{i}.fname);        
        X{n}{i}               = single(Nii.dat(:,:,:)); % Observed data
        msk{n}{i}             = msk_modality(X{n}{i},obj.modality,trunc_ct); % Mask
        X{n}{i} (~msk{n}{i} ) = NaN;        
    end
end

% Compute super-resolved images' orientation matrix and dimensions
if V{1}{1}.dim(3)>1
    pth = {};
    cnt = 1;
    for n=1:N
        I = numel(V{n});
        for i=1:I
            pth{cnt} = V{n}{i}.fname;
            cnt      = cnt + 1;
        end
    end
    [M1,dm1] = max_bb_orient(nifti(char(pth)),vx1);
else
    % Image is 2D (for testing)
    M0  = V{1}{1}.mat;
    vx0 = vxsize(M0);
    dm0 = V{1}{1}.dim;
    d   = vx0./vx1;
    D   = diag([d 1]);
    
    M1     = M0/D;
    dm1    = floor(D(1:3,1:3)*dm0')';
    dm1(3) = 1;
end

% Initialise projection matrices (A)
dat = {};
for n=1:N
    dat{n}.mat      = M1;
    dat{n}.dm       = dm1;
    dat{n}.proj_mat = proj_mat;
    
    dat{n} = init_A(V{n},dat{n});         
end

% Estimate model parameters (tau,lambda)
[tau,lambda] = estimate_model_parameters(V,obj.modality);

% Run algorithm
%--------------------------------------------------------------------------
Y = admm_superres(X,dat,tau,lambda,dm1,vx1,pars_admm);

% Make integer
for n=1:N
    Y{n}             = round(Y{n});
    if strcmp(obj.modality,'MRI')
        Y{n}(Y{n}<0) = 0;
    end
end

% Write results
%--------------------------------------------------------------------------
fnames = {};
cnt    = 1;
for n=1:N
    [pth,nam,ext] = fileparts(V{n}{1}.fname);
    nfname        = fullfile(pth,['sr_' nam ext]);
    create_nii(nfname,Y{n},M1,V{n}{1}.dt,'super-resolved');
    
    V1(n) = spm_vol(nfname);
    
    for i=1:numel(V{n})
        fnames{cnt} = V{n}{i}.fname;
        cnt         = cnt + 1;
    end
    
    fnames{cnt} = nfname;
    cnt         = cnt + 1;
end

if pars_sr.verbose
    % Visualise result
    if V{1}{1}.dim(3)>1
        spm_check_registration(char(fnames))
    else               
        figure(666);
        Nii = nifti(fnames{1});
        img = Nii.dat(:,:,:);
        subplot(121);
        if strcmp(obj.modality,'CT')
            imagesc(img',[0 100]); colormap(gray); axis off xy
        else
            imagesc(img'); colormap(gray); axis off xy
        end         
        
        Nii = nifti(fnames{2});
        img = Nii.dat(:,:,:);
        subplot(122); 
        if strcmp(obj.modality,'CT')
            imagesc(img',[0 100]); colormap(gray); axis off xy
        else
            imagesc(img'); colormap(gray); axis off xy
        end   
    end
end

% Delete input data
for n=1:N
    obj.scans{n}    = {};
    obj.scans{n}{1} = V1(n);
    
    I = numel(V{n});
    for i=1:I 
        delete(V{n}{i}.fname);
    end
end
%==========================================================================

%==========================================================================
function Y = admm_superres(X,dat,tau,lambda,dm,vx,pars)

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

N   = numel(lambda); % Number of modalities
ndm = 3; if dm(3)==1, ndm = 2; end % Number of dimensions

nonneg = true(1,N);
for n=1:N
    I = numel(X{n});
    for i=1:I
        X{n}{i}(~isfinite(X{n}{i})) = 0;    
        if sum(sum(sum(X{n}{i}<0)))
           nonneg(n) = false;
        end
    end
end

% Define Laplace prior (in Fourier space)
% -------------------------------------------------------------------------
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

% Initialize algorithm variables
% -------------------------------------------------------------------------
Y = cell(N,1);
W = cell(N,ndm);
U = cell(N,ndm);
for n=1:N
    Y{n} = zeros(dm,'single');
    for d=1:ndm
        W{n,d} = zeros(dm,'single');
        U{n,d} = zeros(dm,'single');
    end
end
  
% Calculate gradients (gx,gy,gz)
%-------------------------------------------------------------------------
D = cell(N,ndm);
for n=1:N
    [D{n,1:ndm}] = spm_imbasics('grad',Y{n},vx);
end

% Main loop
% -------------------------------------------------------------------------
ll = -Inf;
for iter=1:niter

    % ---------------------------------------------------------------------
    % Sub-problem U
    oU    = U;
    Unorm = single(eps);
    for n=1:N
        for d=1:ndm
            U{n,d} = lambda(n)*D{n,d}+(1/rho)*W{n,d};
            Unorm  = Unorm + U{n,d}.^2;
        end
    end
    Unorm = sqrt(Unorm);
    scale = max(Unorm - 1/rho,0)./Unorm;
    clear Unorm
    
    scale(~isfinite(scale)) = 0;            
    for n=1:N
        for d=1:ndm
            U{n,d} = U{n,d}.*scale;
        end
    end
    clear scale
    
    % Sub-problem Y
    % ---------------------------------------------------------------------
    for n=1:N
        % Compute divergence
        if ndm==3
            DtU = lambda(n)*spm_imbasics('dive',U{n,:},vx);
            DtW = lambda(n)*spm_imbasics('dive',W{n,:},vx);
        else
            DtU = lambda(n)*spm_imbasics('dive',U{n,:},[],vx);
            DtW = lambda(n)*spm_imbasics('dive',W{n,:},[],vx);
        end

        LHS = @(Y) AtA(Y,@(Y) lambda(n)^2*prior(Y),tau{n}/rho,1,dat{n});
        RHS = DtU - (1/rho)*DtW + (1/rho)*At(X{n},tau{n},dat{n});
        clear DtU DtW
        
        % Solve least-squares problem
        Y{n} = cg_image_solver(LHS,RHS,Y{n},cgs_niter,cgs_tol);
        clear RHS LHS
        
        if nonneg(n)
            % Non-negativity 'constraint'
            Y{n}(Y{n}<0) = 0;
        end
    end

    % Sub-problem W
    % ---------------------------------------------------------------------
    for n=1:N
        [D{n,1:ndm}] = spm_imbasics('grad',Y{n},vx);
        for d=1:ndm
            W{n,d} = W{n,d} + rho*(lambda(n)*D{n,d} - U{n,d});
        end
    end

    % Check convergence
    % ---------------------------------------------------------------------
    ll1 = log_likelihood(X,Y,D,dat,tau,lambda);    
    ll  = [ll,ll1];    
    
    d1 = abs((ll(end - 1)*(1 + 10*eps) - ll(end))/ll(end));
    
    if verbose
        fprintf('%2d | %8.7f %8.7f %8.7f\n',iter,ll1,d1,rho);
    end    
    
    if iter>=20 && d1<tol
        % Finished
        break; 
    end        
   
    if est_rho
        % Compute primal and dual residuals
        %------------------------------------------------------------------
        primal_res = 0;
        dual_res   = 0;
        for n=1:N
            tmp = cell(1,ndm);
            for d=1:ndm
                primal_res = primal_res + sum(sum(sum((D{n,d} + U{n,d}).^2)));
                tmp{d}     = U{n,d} - oU{n,d};
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
            for n=1:N
                for d=1:ndm
                    W{n,d} = W{n,d}*scale;
                end
            end
        end  
    end
end
%==========================================================================

%==========================================================================
function ll = log_likelihood(X,Y,D,dat,tau,lambda)
N     = numel(X);
ndm   = numel(size(X{1}));
likel = 0;
prior = 0;
for n=1:N
    AY = A(Y{n},X{n},dat{n});
    for i=1:numel(X{n})
        likel = likel + 0.5*tau{n}(i)*sum(sum(sum((AY{i} - X{n}{i}).^2)));     
    end
    clear AY
    
    for d=1:ndm
        prior = prior + (lambda(n)^2)*(D{n,d}.^2);
    end
end
ll = -likel - sum(sum(sum(sqrt(prior))));
%==========================================================================

%==========================================================================
function X = A(Y,X,dat)
if strcmp(dat.proj_mat,'sinc') 
    Y(~isfinite(Y)) = 0;
    Y               = fftn(Y);    
end
for i=1:dat.I   
    T = dat.mat\dat.A(i).mat;
    y = make_deformation(T,dat.A(i).dm);
    
    if strcmp(dat.proj_mat,'sinc')       
        tmp  = real(ifftn(Y.*dat.A(i).S));                      
    elseif strcmp(dat.proj_mat,'smo')
        tmp = zeros(size(Y),'single');
        spm_smooth(Y,tmp,dat.A(i).S);
    end
    
    X{i} = samp0(tmp,y); 
end
%==========================================================================  

%==========================================================================
function Y = At(X,tau,dat)  
Y = single(0);   
for i=1:dat.I        
    T = dat.mat\dat.A(i).mat;
    y = make_deformation(T,dat.A(i).dm);
    
    tmp = spm_diffeo('push',X{i},y,dat.dm);     
    tmp(~isfinite(tmp)) = 0;
    
    if strcmp(dat.proj_mat,'sinc')               
        tmp                 = real(ifftn(fftn(tmp).*dat.A(i).S));                    
    elseif strcmp(dat.proj_mat,'smo')
        spm_smooth(tmp,tmp,dat.A(i).S); 
    end 
    
    Y = Y + tau(i).*tmp;    
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
for i=1:dat.I        
    T = dat.mat\dat.A(i).mat;
    y = make_deformation(T,dat.A(i).dm);
        
    if strcmp(dat.proj_mat,'sinc')      
        tmp = real(ifftn(F.*dat.A(i).S));    
    elseif strcmp(dat.proj_mat,'smo') 
        tmp = zeros(size(Y),'single');
        spm_smooth(Y,tmp,dat.A(i).S); 
    end
    
    tmp                 = samp0(tmp,y); 
    tmp                 = spm_diffeo('push',tmp,y,dat.dm);          
    tmp(~isfinite(tmp)) = 0;
        
    if strcmp(dat.proj_mat,'sinc')                                          
        tmp = real(ifftn(fftn(tmp).*dat.A(i).S));    
    elseif strcmp(dat.proj_mat,'smo')
        spm_smooth(tmp,tmp,dat.A(i).S);    
    end    
    
    Y1 = Y1 + tau(i).*tmp;
end

Y1 = Y1 + lambda*prior(Y);
%==========================================================================

%==========================================================================
function Y = cg_image_solver(LHS,RHS,Y,niter,tol,verbose)
if nargin<3, Y       = zeros(size(RHS),'single'); end
if nargin<4, niter   = 32; end
if nargin<5, tol     = 1e-3; end
if nargin<6, verbose = false; end

% Initilisation  
%--------------------------------------------------------------------------
normRHS = sqrt(sum(RHS(:).*RHS(:))); % Norm of RHS
R       = RHS - LHS(Y);              % Residual RHS - LHS(x)
normR   = sum(R(:).*R(:));           % R'R
P       = R;                         % Initial conjugate directions P
beta    = 0;                         % Initial search direction for new P

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

    % Finds the step size for updating P
    %----------------------------------------------------------------------
    RtRp  = normR;
    normR = sum(R(:).*R(:));   
    
    if verbose
        fprintf('%g\n', sqrt(normR));
    end
    
    beta = normR / RtRp;            

    % Check if converged
    %----------------------------------------------------------------------
    if j>=niter, 
        % Finished
        break; 
    end;                   

    j = j + 1; 
end
%==========================================================================

%==========================================================================
function [mat,dim] = max_bb_orient(Nii,vx)
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
function dat = init_A(V,dat)
I  = numel(V);
Mj = dat.mat;
dj = dat.dm;
    
dat.I = I;
for i=1:I
    dat.A(i).mat = V{i}.mat;
    dat.A(i).dm  = V{i}.dim;
    
    if strcmp(dat.proj_mat,'sinc')
        Mi         = dat.A(i).mat;
        M          = Mj\Mi;
        R          = (M(1:3,1:3)/diag(sqrt(sum(M(1:3,1:3).^2))))';
        dat.A(i).S = blur_fun(dj,R,[sqrt(sum(M(1:3,1:3).^2))]);
    elseif strcmp(dat.proj_mat,'smo')
        vxj        = abs(det(dat.mat(1:3,1:3)))^(1/3);
        vxi        = sqrt(sum(dat.A(i).mat(1:3,1:3).^2));        
        dat.A(i).S = sqrt(max(vxi.^2 - vxj.^2,0));
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
function y = make_deformation(T,dm)
[x0,y0,z0] = ndgrid(single(1:dm(1)),...
                    single(1:dm(2)),...
                    single(1:dm(3)));
y          = cat(4,T(1,1)*x0 + T(1,2)*y0 + T(1,3)*z0 + T(1,4), ...
                   T(2,1)*x0 + T(2,2)*y0 + T(2,3)*z0 + T(2,4), ...
                   T(3,1)*x0 + T(3,2)*y0 + T(3,3)*z0 + T(3,4));
%==========================================================================                  