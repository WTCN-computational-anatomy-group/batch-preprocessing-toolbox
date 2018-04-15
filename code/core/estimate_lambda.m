function lambda = estimate_lambda(fname,modality,verbose)
if nargin<3, verbose = false; end

if strcmp(modality,'CT')
    lambda = 1.6*1e-3; % Seems empirically like a good value
elseif strcmp(modality,'MRI')
    n = nifti(fname);
    X = single(n.dat(:,:,:));

    k   = 0.4878; % Approximated from IXI high-res images
    msk = isfinite(X) | X ~=0;
    mu  = mean(X(msk));
    
    sd     = k*mu;
    b      = sqrt(sd.^2./2);
    lambda = 1./b;
end

if verbose
    fprintf('lambda = %4.4f\n',lambda);
end
%==========================================================================