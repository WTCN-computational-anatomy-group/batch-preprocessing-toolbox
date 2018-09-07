function nii_res = reslice_labels(V_reference,nii_labels)

V    = spm_vol;
V(1) = V_reference;
dm0  = V_reference.dim;

labels = nii_labels.dat(:,:,:);
lkp    = unique(labels)';
Kl     = numel(lkp);

fname         = nii_labels.dat.fname;
[pth,nam,ext] = fileparts(fname);

nlabels = zeros([dm0 Kl]);
cnt     = 1;
for k=lkp
    bin_labels = single(labels==k);
    
    nfname = fullfile(pth,['tmp-' nam ext]);    
    spm_misc('create_nii',nfname,bin_labels,nii_labels.mat,nii_labels.dat.dtype,nii_labels.descrip,nii_labels.dat.offset,nii_labels.dat.scl_slope,nii_labels.dat.scl_inter);        
    
    % Reslice
    V(2)  = spm_vol(nfname);
    rV    = spm_impreproc('reslice',V,1,1);
    rlabels = rV(2).private.dat(:,:,:);
    nlabels(:,:,:,cnt) = cat(4,rlabels);
    
    cnt = cnt + 1;
end
delete(fname);
% figure; imshow3D(squeeze(nlabels(:,:,floor(dm0(3)/2) + 1,:)))

[~,nlabels] = max(nlabels,[],4);
nlabels     = nlabels - 1;

nii_res            = nifti(rV(2).fname);
nii_res.dat(:,:,:) = nlabels;
%==========================================================================