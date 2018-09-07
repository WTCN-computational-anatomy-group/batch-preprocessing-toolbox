function collapse_labels(nii,part)
K       = numel(part);
labels  = nii.dat(:,:,:);
nlabels = zeros(size(labels));
for k=1:K
    p            = part{k}';
    msk          = ismember(labels,p);
    nlabels(msk) = k - 1;
end
nii.dat(:,:,:) = nlabels;
%==========================================================================