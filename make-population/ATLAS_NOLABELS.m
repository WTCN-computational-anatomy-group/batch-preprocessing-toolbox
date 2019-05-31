clear; clc;

addpath('/data/mbrud/dev/auxiliary-functions/');

nf = '/home/mbrud/Data/images/ATLAS/ATLAS-NOLABELS/';
f  = '/home/mbrud/Data/images/ATLAS/native_1/';

files = spm_select('FPListRec',f,'^.*\.nii.gz$');
S     = size(files,1);
s0    = 0;

%% Unzip
for s=1:S
    fname = strtrim(files(s,:));        
    gunzip(fname,nf);
end

%% Get names of all subjects
files = spm_select('FPList',nf,'^.*\.nii$');
S     = size(files,1);
names = {};
for s=1:S
    fname   = strtrim(files(s,:));
    [~,nam] = fileparts(fname);
    name    = nam(1:13);
    
    if ~any(strcmp(names,name))        
        names{end + 1} = name;
    end
end

% %% Sum up labels
% S = numel(names);
% pths = cell(S,2);
% for s=1:S
%     name = names{s};
%     
%     pth_img = fullfile(nf,[name '.nii']);
%     nii_img = nifti(pth_img);
%     dm_img  = size(nii_img.dat(:,:,:));
%     labels  = zeros(dm_img);
%     
%     pths{s,1} = pth_img;
%     
%     % Create label image
%     pth_lab = fullfile(nf,[name '_labels' '.nii']);  
%     nii_lab = nifti(fullfile(nf,[name '_LesionSmooth' '.nii']));
%     spm_misc('create_nii',pth_lab,labels,nii_lab.mat,nii_lab.dat.dtype,nii_lab.descrip,nii_lab.dat.offset,nii_lab.dat.scl_slope,nii_lab.dat.scl_inter);
%     
%     pths{s,2} = pth_lab;
%     
%     % Get names{s} labels
%     pth_labels = spm_select('FPList',nf,['^' name '_LesionSmooth.*']);
%     L          = size(pth_labels,1);
%     for l=1:L
%         fname   = strtrim(pth_labels(l,:));
%         nii_lab = nifti(fname);        
%         dm_lab  = size(nii_lab.dat(:,:,:));
%         
%         if isequal(dm_img,dm_lab)
%             labels = labels + nii_lab.dat(:,:,:);
%         end
%         
%         delete(fname);
%     end
%     
%     labels = labels>0;
%     
%     if sum(labels(:))==0
%         s
%         warning('sum(labels(:))==0')
%     end
%     
%     nii_lab            = nifti(pth_lab);   
%     nii_lab.dat(:,:,:) = labels;
% end

%%
for s=1:S
    % Image
    fname       = strtrim(files(s,:));
    [~,nam,ext] = fileparts(fname);
    
    a            = struct;
    a.name       = names{s};
    a.population = 'ATLAS';
    a.modality   = 'MRI';
    a.channel    = 'T1';
    a.lesion     = true;    
    a.pth        = [nam ext];

    a = orderfields(a);

    pth_json = fullfile(nf,[nam '.json']);
    spm_jsonwrite(pth_json,a);
end

%%
dat = spm_json_manager('init_dat',nf);

s  = 100;
f1 = dat{s}.modality{1}.channel{1}.nii.dat.fname;

spm_check_registration(char({f1}))