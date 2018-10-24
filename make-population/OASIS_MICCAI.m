clear; clc;

f  = '/home/mbrud/Data/labels/OASIS-MICCAI';
nf = '/data/mbrud/populations/original/OASIS-MICCAI';

files = spm_select('FPList',f,'^.*\.nii$');
S     = size(files,1);

%%
for s=1:2:S
    fname       = strtrim(files(s,:));        
    [~,nam,ext] = fileparts(fname);
    
    name = nam(1:4);
    
    % Image
    nfname = fullfile(nf,[nam ext]);
    copyfile(fname,nfname);
    
    a            = struct;
    a.name       = name;
    a.population = 'OASIS-MICCAI';
    a.modality   = 'MRI';
    a.channel    = 'T1';
    a.healthy    = false;
    a.pth        = [nam ext];

    a = orderfields(a);

    pth_json = fullfile(nf,[name '.json']);
    spm_jsonwrite(pth_json,a);
    
    % Labels
    fname_lab   = strtrim(files(s + 1,:));   
    [~,nam,pth] = fileparts(fname_lab);
    
    nfname = fullfile(nf,[nam ext]);
    copyfile(fname_lab,nfname);
    
    a              = struct;
    a.name         = name;
    a.rater        = 'unknown';
    a.population   = 'OASIS-MICCAI';
    a.nam_modality = 'MRI';
    a.nam_channel  = 'T1';
    a.ix_img       = 1;
    a.pth          = [nam ext];

    a = orderfields(a);

    pth_json = fullfile(nf,[name '-seg' '.json']);
    spm_jsonwrite(pth_json,a);
end