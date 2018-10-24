clear; clc;

addpath('/data/mbrud/dev/auxiliary-functions/');

f  = '/home/mbrud/Data/labels/CROMIS-Parashkev-manual/converted_all';

files = spm_select('FPList',f,'^.*\.nii$');
S     = size(files,1);

%%
for s=1:2:S
    % Image
    fname       = strtrim(files(s,:));
    [~,nam0,ext] = fileparts(fname);
    
    a            = struct;
    a.name       = nam0;
    a.population = 'CROMIS-LABELS';
    a.modality   = 'CT';
    a.healthy    = false;    
    a.pth        = [nam0 ext];

    a = orderfields(a);

    pth_json = fullfile(f,[nam0 '.json']);
    spm_jsonwrite(pth_json,a);
    
    % Labels
    fname       = strtrim(files(s + 1,:));
    [~,nam,ext] = fileparts(fname);
    
    a              = struct;
    a.name         = nam0;     
    a.rater        = 'unknown';
    a.population   = 'CROMIS-LABELS';
    a.nam_modality = 'CT';
    a.ix_img       = 1;
    a.pth          = [nam ext];
    
    a = orderfields(a);
                
    pth_json = fullfile(f,[nam '.json']);
    spm_jsonwrite(pth_json,a);
end

%%
dat = spm_json_manager('init_dat',f);

s  = 30;
f1 = dat{s}.modality{1}.nii.dat.fname;
f2 = dat{s}.label{1}.nii.dat.fname;

spm_check_registration(char({f1,f2}))