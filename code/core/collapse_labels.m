function collapse_labels(nii,part)
K       = numel(part);
labels  = round(nii.dat(:,:,:));
if ischar(part)
    p   = unique(labels);
    nlabels = ismember(labels,p(2:end));
else
    nlabels = zeros(size(labels));
    for k=1:K
        if iscell(part)
            p = part{k}';
        else
            p = part(k)';
        end
        msk          = ismember(labels,p);
        nlabels(msk) = k - 1;
    end
end
nii.dat(:,:,:) = nlabels;
%==========================================================================