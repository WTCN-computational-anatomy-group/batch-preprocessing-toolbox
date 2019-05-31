clear; clc;

f  = '/home/mbrud/Desktop/CROMIS-CAL-LABELS';
nf = '/data/mbrud/populations/original/CROMIS-LABELS-CAL';

if exist(nf,'dir') == 7, rmdir(nf,'s'); end; mkdir(nf); 

files = spm_select('FPList',f,'^.*\.nii$');
S     = size(files,1);

%%
for s=1:3:S
    fprintf('.')
    
    f_cal = strtrim(files(s,:));
    f_img = strtrim(files(s + 1,:));
    f_les = strtrim(files(s + 2,:));
    
    % Set affine of f_cal to f_img
    Nii = nifti(f_img);
    mat = Nii.mat;    
    spm_get_space(f_cal,mat);   
    
    % Sum up mask
    Nii_cal = nifti(f_cal);
    Nii_les = nifti(f_les);
    
    msk      = Nii_cal.dat(:,:,:) > 0;
    img      = Nii_les.dat(:,:,:);
    img(msk) = 2;
    
    Nii_les.dat(:,:,:) = img;
    
    delete(f_cal);
end
fprintf('\n')

%%
files = spm_select('FPList',f,'^.*\.nii$');
S     = size(files,1);

for s=1:2:S
    fprintf('.')
    
    f_img = strtrim(files(s,:));
    f_les = strtrim(files(s + 1,:));

    % Image
    fname        = f_img;
    [~,nam0,ext] = fileparts(fname);
    nfname       = fullfile(nf,[nam0 ext]);
    copyfile(fname,nfname);
    
    a            = struct;
    a.name       = nam0;
    a.population = 'CROMIS-LABELS';
    a.modality   = 'CT';
    a.healthy    = false;    
    a.pth        = [nam0 ext];

    a = orderfields(a);

    pth_json = fullfile(nf,[nam0 '.json']);
    spm_jsonwrite(pth_json,a);

    % Labels
    fname        = f_les;
    [~,nam1,ext] = fileparts(fname);
    nfname       = fullfile(nf,[nam1 ext]);
    copyfile(fname,nfname);
    
    a              = struct;
    a.name         = nam0;     
    a.rater        = 'unknown';
    a.population   = 'CROMIS-LABELS';
    a.nam_modality = 'CT';
    a.ix_img       = 1;
    a.pth          = [nam1 ext];
    
    a = orderfields(a);
                
    pth_json = fullfile(nf,[nam1 '.json']);
    spm_jsonwrite(pth_json,a);
end
fprintf('\n')

%%
addpath('/home/mbrud/dev/ToolBoxes/auxiliary-functions');

dat = spm_json_manager('init_dat',nf);

s  = 30;
f1 = dat{s}.modality{1}.nii.dat.fname;
f2 = dat{s}.label{1}.nii.dat.fname;

spm_check_registration(char({f1,f2}))