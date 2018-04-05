function obj = spm_denoise(obj)

% Parameters
%--------------------------------------------------------------------------
V         = obj.scans;
N         = numel(V);
vx        = vxsize(V{1}{1}.mat);
pars_den  = obj.preproc.denoise;
pars_admm = pars_den.admm;

% Estimate model parameters (tau,lambda)
%--------------------------------------------------------------------------
tau    = zeros(1,N);
lambda = zeros(1,N);
for n=1:N
    tau(n)    = estimate_tau(V{n}{1}.fname,obj.modality);
    lambda(n) = estimate_lambda(V{n}{1}.fname,obj.modality);
end

% Get image data
%--------------------------------------------------------------------------
X   = cell(1,N);
msk = cell(1,N);
for n=1:N
   X{n}          = single(V{n}{1}.private.dat(:,:,:)); % Observed data
   msk{n}        = msk_modality(X{n},obj.modality); % Mask
   X{n}(~msk{n}) = NaN;
end

% Run algorithm
%--------------------------------------------------------------------------
Yhat = admm_denoising(X,tau,lambda,vx,pars_admm);

% Re-apply mask
for n=1:N
   Yhat{n}(~msk{n}) = NaN;
end

% Make integer again
for n=1:N
    Yhat{n} = round(Yhat{n});
    if strcmp(obj.modality,'MRI')
        Yhat{n}(Yhat{n}<0) = 0;
    end
end

% Write results
%--------------------------------------------------------------------------
fnames = {};
cnt    = 1;
for n=1:N
    [pth,nam,ext] = fileparts(V{n}{1}.fname);
    nfname        = fullfile(pth,['den_' nam ext]);
    create_nii(nfname,Yhat{n},V{n}{1}.mat,V{n}{1}.dt,'denoised');
    
    V1(n) = spm_vol(nfname);
    
    for i=1:numel(V{n})
        fnames{cnt} = V{n}{i}.fname;
        cnt         = cnt + 1;
    end
    
    fnames{cnt} = nfname;
    cnt         = cnt + 1;
end

if pars_den.verbose
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
function Y = admm_denoising(X,tau,lambda,vx,pars)

% Set parametersWrite results
%--------------------------------------------------------------------------
niter   = pars.niter;
tol     = pars.tol;
rho     = single(pars.rho);
verbose = pars.verbose;
mu      = pars.mu;
alpha   = pars.alpha;
est_rho = pars.est_rho;

dm = size(X{1});
if numel(dm)==2, dm(3) = 1; end

N   = numel(lambda); % Number of modalities
ndm = 3; if dm(3)==1, ndm = 2; end % Number of dimensions

nonneg = true;
for n=1:N
    X{n}(~isfinite(X{n})) = 0;    
    if sum(sum(sum(X{n}<0)))
       nonneg = false;
    end
end

% Define Laplacian filter (in Fourier space)
%--------------------------------------------------------------------------
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

% Initialize algorithm variables
%--------------------------------------------------------------------------
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
%--------------------------------------------------------------------------
ll = -Inf;
for iter=1:niter

    %----------------------------------------------------------------------
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
    %----------------------------------------------------------------------
    for n=1:N
        % Compute divergence
        if ndm==3
            DtU = lambda(n)*spm_imbasics('dive',U{n,:},vx);
            DtW = lambda(n)*spm_imbasics('dive',W{n,:},vx);
        else
            DtU = lambda(n)*spm_imbasics('dive',U{n,:},[],vx);
            DtW = lambda(n)*spm_imbasics('dive',W{n,:},[],vx);
        end

        RHS = DtU - (1/rho)*DtW + (tau(n)/rho)*X{n};    
        clear DtU DtW
        
        % Solve least-squares problem
        Y{n} = real(ifftn(fftn(RHS)./(L*(lambda(n)^2) + (tau(n)/rho))));
        clear RHS
        
        if nonneg
            % Non-negativity 'constraint'
            Y{n}(Y{n}<0) = 0;
        end
    end

    % Sub-problem W
    %----------------------------------------------------------------------
    for n=1:N
        [D{n,1:ndm}] = spm_imbasics('grad',Y{n},vx);
        for d=1:ndm
            W{n,d} = W{n,d} + rho*(lambda(n)*D{n,d} - U{n,d});
        end
    end

    % Check convergence
    %----------------------------------------------------------------------
    ll1 = log_likelihood(X,Y,D,tau,lambda);    
    ll  = [ll,ll1];    
    
    d1 = abs((ll(end - 1)*(1 + 10*eps) - ll(end))/ll(end));
    
    if verbose
        fprintf('%2d | %8.7f %8.7f %8.7f %8.7f\n',iter,ll1,d1,rho,tau);
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
function ll = log_likelihood(X,Y,D,tau,lambda)
N     = numel(X);
ndm   = numel(size(X{1}));
likel = 0;
prior = 0;
for n=1:N
    likel = likel + 0.5*tau(n)*sum(sum(sum((Y{n} - X{n}).^2)));     
    for d=1:ndm
        prior = prior + (lambda(n)^2)*(D{n,d}.^2);
    end
end
ll = -likel - sum(sum(sum(sqrt(prior))));
%==========================================================================