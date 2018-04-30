function tau = estimate_tau(fname,modality)
Nii = nifti(fname);
X   =  single(Nii.dat(:,:,:));
vx  = spm_misc('vxsize',Nii.mat);

if strcmp(modality,'CT')
    msk = X>-1010 & X<-990;
    sd  = std(X(msk));
elseif strcmp(modality,'MRI')
    sd  = my_spm_noise_estimate(fname);     
end

tau = 1/(2*sd^2)*1/prod(vx);
%==========================================================================

%==========================================================================
function SmoSuf = ComputeSmoSuf(img,vx)
SmoSuf    = zeros(1,8);
tmp       = img;
msk       = isfinite(tmp(:));
SmoSuf(1) = sum(msk(:));
SmoSuf(2) = sum(tmp(msk(:)).^2);

d         = size(img);
[id{1:3}] = ndgrid(single(1:d(1)),single(1:d(2)),single(1:d(3)));
id        = cat(4,id{:});

% fc           = spm_diffeo('bsplinc',img,[1 1 1 1 1 1]);
% [~,Dx,Dy,Dz] = spm_diffeo('bsplins',fc,id,[1 1 1 1 1 1]);

[Dx,Dy,Dz] = spm_imbasics('grad',tmp,vx);

tmp       = Dx;
msk       = isfinite(tmp(:));
SmoSuf(3) = sum(msk(:));
SmoSuf(4) = sum(tmp(msk(:)).^2);

tmp       = Dy;
msk       = isfinite(tmp(:));
SmoSuf(5) = sum(msk(:));
SmoSuf(6) = sum(tmp(msk(:)).^2);

tmp       = Dz;
msk       = isfinite(tmp(:));
SmoSuf(7) = sum(msk(:));
SmoSuf(8) = sum(tmp(msk(:)).^2);
%==========================================================================

%==========================================================================
function SmoSuf = ComputeSmoSuf_slice(img,vx)
d         = size(img);
d         = [d 1];
[id{1:3}] = ndgrid(single(1:d(1)),single(1:d(2)),single(1));
id        = cat(4,id{:});
SmoSuf    = zeros(1,8);
for z=1:d(3)
    tmp = img(:,:,z);
    
    msk       = isfinite(tmp(:));
    SmoSuf(1) = SmoSuf(1) + sum(msk(:));
    SmoSuf(2) = SmoSuf(2) + sum(tmp(msk(:)).^2);
    
%     fc        = spm_diffeo('bsplinc',tmp,[1 1 0 1 1 1]);
%     [~,Dx,Dy] = spm_diffeo('bsplins',fc,id,[1 1 0 1 1 1]);

    [Dx,Dy] = spm_imbasics('grad',tmp,vx);
    
    tmp       = Dx;
    msk       = isfinite(tmp(:));
    SmoSuf(3) = SmoSuf(3) + sum(msk(:));
    SmoSuf(4) = SmoSuf(4) + sum(tmp(msk(:)).^2);

    tmp       = Dy;
    msk       = isfinite(tmp(:));
    SmoSuf(5) = SmoSuf(5) + sum(msk(:));
    SmoSuf(6) = SmoSuf(6) + sum(tmp(msk(:)).^2);
end
%==========================================================================

%==========================================================================
function fwhm = compute_fwhm(SmoSuf)
fwhm = sqrt(4*log(2)*(SmoSuf(2)/SmoSuf(1))/(sum(SmoSuf([4 6 8]))/sum(SmoSuf([3 5 7]))));
fwhm = sqrt(max(fwhm^2 - 0.5, 4*log(2)/pi)); 
%==========================================================================

%==========================================================================
function [tau,sd] = compute_tau(fwhm,vx)
if nargin<2, vx = [1 1 1]; end

fwhm = fwhm + mean(vx); 
sd   = fwhm/sqrt(8*log(2));  % Standard deviation
tau  = prod(4*pi*(sd./vx).^2 + 1)^(-1/2); 
%==========================================================================