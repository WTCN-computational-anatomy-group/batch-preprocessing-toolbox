clear; clc;

addpath('/data/mbrud/dev/auxiliary-functions/');

f0 = '/home/mbrud/Data/labels/CROMIS-Parashkev-manual/labels';
d  = dir(f0);
d  = d(3:end);
S  = numel(d);

fn = '/data/mbrud/populations/original/CROMIS-LABELS/';

%%
for s=1:S
    
    ds    = fullfile(f0,d(s).name);    
    files = spm_select('FPList',ds,'^.*\.nii$');
    
    % Image
    fname        = strtrim(files(1,:));
    [~,nam0,ext] = fileparts(fname);
    nfname       = fullfile(fn,[nam0 ext]);
    copyfile(fname,nfname);
    
    a            = struct;
    a.name       = nam0;
    a.population = 'CROMIS-LABELS';
    a.modality   = 'CT';
    a.healthy    = false;    
    a.pth        = [nam0 ext];

    a = orderfields(a);

    pth_json = fullfile(fn,[nam0 '.json']);
    spm_jsonwrite(pth_json,a);
    
    % Labels
    fname        = strtrim(files(2,:));
    [~,nam1,ext] = fileparts(fname);
    nfname       = fullfile(fn,[nam1 ext]);
    copyfile(fname,nfname);
    
    a              = struct;
    a.name         = nam0;     
    a.rater        = 'unknown';
    a.population   = 'CROMIS-LABELS';
    a.nam_modality = 'CT';
    a.ix_img       = 1;
    a.pth          = [nam1 ext];
    
    a = orderfields(a);
                
    pth_json = fullfile(fn,[nam1 '.json']);
    spm_jsonwrite(pth_json,a);
end

%%
dat = spm_json_manager('init_dat',fn);

%%
s  = 30;
f1 = dat{s}.modality{1}.nii.dat.fname;
f2 = dat{s}.label{1}.nii.dat.fname;

spm_check_registration(char({f1,f2}))