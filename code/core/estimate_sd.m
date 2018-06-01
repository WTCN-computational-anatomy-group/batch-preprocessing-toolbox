function [sd,X] = estimate_sd(fname,modality)
Nii = nifti(fname);
X   = Nii.dat(:,:,:);

if strcmpi(modality,'CT')    
    msk     = X>=0 & X<=50;
    X(~msk) = 0;    
    noise   = evar(X);
    sd      = sqrt(noise);
elseif strcmpi(modality,'MRI')
    sd  = my_spm_noise_estimate(fname);     
end
%==========================================================================

%==========================================================================
function sd = calc_sd(x)
x  = x(:);
N  = numel(x);
m  = 1/N*sum(x);
sd = sqrt(1/(N - 1)*sum(abs(x - m).^2));
%==========================================================================