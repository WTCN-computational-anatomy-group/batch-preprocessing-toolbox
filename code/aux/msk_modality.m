function msk = msk_modality(f,modality,trunc_ct)
if nargin<3, trunc_ct = []; end

if strcmp(modality,'MRI'),    
    msk = isfinite(f) & (f~=0);
elseif strcmp(modality,'CT'), 
    if isempty(trunc_ct)
        msk = isfinite(f) & (f~=min(f(:))) & (f~=0);         
    else
        msk = isfinite(f) & (f>trunc_ct(1)) & (f<trunc_ct(2));          
    end    
    
    msk = imfill(msk,'holes'); % Because there might be 0 values voxels within the brain that gets masked out above
end
%==========================================================================

