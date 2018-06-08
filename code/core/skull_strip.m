function dat = skull_strip(dat)

% For segmented data we assume... 
m = 1; % ...one modality...
n = 1; % ...and one image per channel   

tissue = {'GM','WM','CSF'};
K      = numel(tissue);
type   = 'c';
V0     = spm_vol;
t      = dat.segmentation_map(type);
for i=1:K
    k     = dat.segmentation{t}.class_map(tissue{i});
    V0(k) = spm_vol(dat.segmentation{t}.class{k}.nii.dat.fname);
end

msk   = zeros(V0(1).dim,'single');
for k=1:K
    Nii  = nifti(V0(k).fname);
    resp = single(Nii.dat(:,:,:)); 
    msk  = msk + resp;
end
clear resp

for z=1:V0(1).dim(3) % loop over slices
    msk(:,:,z) = imgaussfilt(msk(:,:,z),1);    % Smooth
    msk(:,:,z) = msk(:,:,z)>0.5;               % Threshold
    msk(:,:,z) = imfill(msk(:,:,z),4,'holes'); % Fill holes
end

% Mask out voxels based on SPM TPM size
pth_tpm = fullfile(spm('dir'),'tpm','TPM.nii,');
V1      = spm_vol(pth_tpm);

M0  = V0(1).mat;      
dm0 = V0(1).dim; 
M1  = V1(1).mat;  
dm1 = V1(1).dim;

T = M1\M0; % Get the mapping from M0 to M1

% Use ndgrid to give an array of voxel indices
[x0,y0,z0] = ndgrid(single(1:dm0(1)),...
                    single(1:dm0(2)),...
                    single(1:dm0(3)));

% Transform these indices to the indices that they point to in the reference image
D = cat(4,T(1,1)*x0 + T(1,2)*y0 + T(1,3)*z0 + T(1,4), ...
          T(2,1)*x0 + T(2,2)*y0 + T(2,3)*z0 + T(2,4), ...
          T(3,1)*x0 + T(3,2)*y0 + T(3,3)*z0 + T(3,4));
clear x0 y0 z0

% Mask according to whether these are < 1 or > than the dimensions of the reference image.        
msk1 = cell(1,3);
ix   = [1 1 20];
for i=1:3
    msk1{i} = D(:,:,:,i) >= ix(i) & D(:,:,:,i) <= dm1(i);
end
clear D

% Generate masked image
for i=1:3
    msk = msk1{i}.*msk;
end

% Overwrite image data with skull-stripped version
if isfield(dat.modality{m},'channel')            
    C = numel(dat.modality{m}.channel);
    for c=1:C            
        Nii            = dat.modality{m}.channel{c}.nii(n);
        img            = single(Nii.dat(:,:,:));
        img(~msk)      = NaN;  
        Nii.dat(:,:,:) = img;                    
    end
else        
    Nii            = dat.modality{m}.nii(n);
    img            = single(Nii.dat(:,:,:));
    img(~msk)      = NaN;  
    Nii.dat(:,:,:) = img;  
end 
%==========================================================================        