function Niio = reslice_labels(V_ref,Nii)

% Parameters
deg     = 1;
ref     = 1;
dt      = [spm_type('float32') spm_platform('bigend')];

% Create spm_vol object and put reference at index 1
V    = spm_vol;
V(1) = V_ref;
dm   = V_ref.dim;

% Get labels, etc
labels        = Nii.dat(:,:,:);
lkp           = unique(labels)';
K             = numel(lkp); % Number of labels
fname         = Nii.dat.fname;
[pth,nam,ext] = fileparts(fname);

% Iterate over each label and create a resliced label image (nlabels)
labelso = zeros([dm K],'single');
cnt     = 1;
for k=lkp
    labels_k = single(labels == k);

    % Smooth
    labels_k = convn(labels_k,reshape([0.25 0.5 0.25],[3,1,1]),'same');
    labels_k = convn(labels_k,reshape([0.25 0.5 0.25],[1,3,1]),'same');
    labels_k = convn(labels_k,reshape([0.25 0.5 0.25],[1,1,3]),'same');
                
    fname_k = fullfile(pth,['n' nam ext]);    
    spm_misc('create_nii',fname_k,labels_k,Nii.mat,dt,'Resliced labels');        
        
    % Reslice
    V(2)     = spm_vol(fname_k);
    Vo       = spm_impreproc('reslice',V,deg,ref);
    labels_k = single(Vo(2).private.dat(:,:,:));

    labelso(:,:,:,cnt) = cat(4,labels_k);
    
    cnt = cnt + 1;
end
clear labels labels_k
delete(fname);

% Get MLs of resliced labels
[~,ml] = max(labelso,[],4);
ml     = ml - 1;

% Write output
Niio            = nifti(Vo(2).fname);
Niio.dat(:,:,:) = ml;
%==========================================================================