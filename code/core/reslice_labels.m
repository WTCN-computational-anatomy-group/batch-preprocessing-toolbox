function nii_res = reslice_labels(V_reference,nii_labels)

V    = spm_vol;
V(1) = V_reference;
dm0  = V_reference.dim;

labels = nii_labels.dat(:,:,:);

% figure; imshow3D(labels)

lkp           = unique(labels)';
lkp(lkp == 0) = [];
Kl            = numel(lkp);

fname         = nii_labels.dat.fname;
[pth,nam,ext] = fileparts(fname);

nlabels = zeros([dm0 Kl],'single');
cnt     = 1;
for k=lkp
    bin_labels = single(labels==k);
%     bin_labels = bin_labels + 1e-2*rand(size(bin_labels));
%     spm_imbasics('smooth_img_in_mem',bin_labels,1);
    
    nfname = fullfile(pth,['tmp-' nam ext]);    
    spm_misc('create_nii',nfname,bin_labels,nii_labels.mat,nii_labels.dat.dtype,nii_labels.descrip,nii_labels.dat.offset,nii_labels.dat.scl_slope,nii_labels.dat.scl_inter);        
    clear bin_labels
    
    % Reslice
    V(2)  = spm_vol(nfname);
    rV    = spm_impreproc('reslice',V,0,1);
    rlabels = rV(2).private.dat(:,:,:);
    
    nlabels(:,:,:,cnt) = cat(4,rlabels);
    
    cnt = cnt + 1;
end
delete(fname);
% figure; imshow3D(squeeze(nlabels(:,:,floor(dm0(3)/2) + 1,:)))

dm      = size(nlabels);
nlabels = reshape(nlabels,[prod(dm(1:3)) dm(4)]);
msk     = sum(nlabels,2) > 0;

% figure(3); imshow3D(reshape(msk,dm(1:3)))

nlabels1 = zeros([nnz(msk) dm(4)],'single');
for k=1:dm(4)
    nlabels1(:,k) = nlabels(msk,k);
end

[~,mllabels] = max(nlabels1,[],2);
% mllabels     = mllabels - 1;

mllabels1      = zeros([prod(dm(1:3)) 1],'single');
mllabels1(msk) = mllabels;

mllabels1 = reshape(mllabels1,dm(1:3));

% figure; imshow3D(mllabels1)

nii_res            = nifti(rV(2).fname);
nii_res.dat(:,:,:) = mllabels1;
%==========================================================================