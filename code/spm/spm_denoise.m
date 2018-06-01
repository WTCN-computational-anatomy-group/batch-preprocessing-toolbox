function nii_out = spm_denoise(nii,modality,opt)

% Sanity-check for denoising
%--------------------------------------------------------------------------
vx = spm_misc('vxsize',nii(1).mat);
dm = size(nii(1).dat(:,:,:));
dm = [dm 1];
C  = numel(nii);
for c=2:C
    vx1 = spm_misc('vxsize',nii(c).mat);   
    if ~isequal(vx,vx1)
        error('isequal(vx,vx1)')
    end

    dm1 = size(nii(c).dat(:,:,:));
    if ~isequal(dm,dm1)
        error('isequal(dm,dm1)')
    end
end

if ~ismember(modality,{'CT','MRI'})
    error('~ismember(modality, {''CT'',''MRI''})')
end
    
% Estimate model parameters (tau,lambda)
%--------------------------------------------------------------------------
sd     = zeros(1,C);
lambda = zeros(1,C);
for c=1:C    
    sd(c)     = estimate_sd(nii(c).dat.fname,modality);       
    lambda(c) = estimate_lambda(nii(c).dat.fname,modality,opt.lambda_ct);
end

% Get image data
%--------------------------------------------------------------------------
X   = cell(1,C);
msk = cell(1,C);
for c=1:C   
   X{c}          = single(nii(c).dat(:,:,:)); % Observed data
   msk{c}        = spm_misc('msk_modality',X{c},modality); % Mask
   X{c}(~msk{c}) = NaN;
end

% Run algorithm
%--------------------------------------------------------------------------
Y = admm_denoising(X,sd,lambda,vx,modality,opt.admm);

% Re-apply mask
%--------------------------------------------------------------------------
for c=1:C
   Y{c}(~msk{c}) = NaN;
end

% Make integer again
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
for c=1:C
    fname           = nii(c).dat.fname;
    fnames{end + 1} = fname;    

    [pth,nam,ext]   = fileparts(fname);
    nfname          = fullfile(pth,['den_' nam ext]);
    spm_misc('create_nii',nfname,Y{c},nii(c).mat,nii(c).dat.dtype,'denoised');            
    nii_out(c)      = nifti(nfname);        
    fnames{end + 1} = nfname;
end

% Some verbose
%--------------------------------------------------------------------------
if opt.verbose  
    spm_check_registration(char(fnames))
end

% Delete input data
%--------------------------------------------------------------------------
for c=1:C
    fname = nii(c).dat.fname;
    delete(fname);
end
%==========================================================================

%==========================================================================
function Y = admm_denoising(X,sd,lambda,vx,modality,pars)

% Set parameters
%--------------------------------------------------------------------------
niter   = pars.niter;
tol     = pars.tol;
rho     = single(pars.rho);
verbose = pars.verbose;
mu      = pars.mu;
alpha   = pars.alpha;
est_rho = pars.est_rho;
tau     = 1./(sd.^2);

dm = size(X{1}); % Dimensions
if numel(dm)==2, 
    dm(3) = 1; 
end
ndm = 3; % Number of dimensions
if dm(3)==1, 
    ndm = 2; 
end 
C = numel(lambda); % Number of channels

nonneg = true; % Image contains negative values? If not, impose non-negativity constraint
for c=1:C
    X{c}(~isfinite(X{c})) = 0;    
    if sum(sum(sum(X{c}<0)))
       nonneg = false;
    end
end

% Get total number of observed voxels
nm = 0;
for c=1:C
    nm = nm + numel(X{c});
end

% Define Laplacian filter (in Fourier space)
%--------------------------------------------------------------------------
L = laplace_kernel(dm,vx);

% Initialize algorithm variables
%--------------------------------------------------------------------------
Y = cell(C,1);
W = cell(C,ndm);
U = cell(C,ndm);
for c=1:C
    Y{c} = zeros(dm,'single');
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
%--------------------------------------------------------------------------
ll = -Inf;
for iter=1:niter

    %----------------------------------------------------------------------
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
    %----------------------------------------------------------------------
    for c=1:C
        % Compute divergence
        if ndm==3
            DtU = lambda(c)*spm_imbasics('dive',U{c,:},vx);
            DtW = lambda(c)*spm_imbasics('dive',W{c,:},vx);
        else
            DtU = lambda(c)*spm_imbasics('dive',U{c,:},[],vx);
            DtW = lambda(c)*spm_imbasics('dive',W{c,:},[],vx);
        end

        RHS = DtU - (1/rho)*DtW + (tau(c)/rho)*X{c};    
        clear DtU DtW
        
        % Solve least-squares problem
        Y{c} = real(ifftn(fftn(RHS)./(L*(lambda(c)^2) + (tau(c)/rho))));
        clear RHS
        
        if nonneg
            % Non-negativity 'constraint'
            Y{c}(Y{c}<0) = 0;
        end
    end

    % Sub-problem W
    %----------------------------------------------------------------------
    for c=1:C
        [D{c,1:ndm}] = spm_imbasics('grad',Y{c},vx);
        for d=1:ndm
            W{c,d} = W{c,d} + rho*(lambda(c)*D{c,d} - U{c,d});
        end
    end

    % Check convergence
    %----------------------------------------------------------------------
    ll1 = log_likelihood(X,Y,D,tau,lambda);    
    ll  = [ll,ll1];    
    
    d1 = abs((ll(end - 1)*(1 + 10*eps) - ll(end))/ll(end));
    d2 = ll(end) - ll(end - 1);
        
    if verbose
        fprintf('%2d | %8.7f %8.7f %8.7f %8.7f %8.7f\n',iter,ll1,d1,d2,1e-4*nm,rho);
                                
        show_img(X,Y,lambda,sd,modality,'Denoising');
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
function ll = log_likelihood(X,Y,D,tau,lambda)
C   = numel(X);
ndm = numel(size(X{1}));
ll1 = 0;
llr = 0;
for c=1:C
    I     = numel(X{c});
    sd    = sqrt(1/tau(c));
    diff1 = (Y{c} - X{c});
    diff1 = diff1(:)'*diff1(:);
    ll1   = ll1 + 0.5*(2*I*log(sd) + I*log(2*pi) + 1/sd^2*diff1);
    
    for d=1:ndm
        llr = llr + (lambda(c)^2)*(D{c,d}.^2);
    end
end
ll = -ll1 - sum(sum(sum(sqrt(llr))));
%==========================================================================

%==========================================================================
function L = laplace_kernel(dm,vx)
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
L = single(fftn(L));
%==========================================================================     