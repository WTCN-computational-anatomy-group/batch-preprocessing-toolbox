clear; clc;

[files,dirs] = spm_select('FPListRec','/data/mbrud/populations/original/ROB','^.*\.nii$');
ndir         = '/data/mbrud/populations/original/ROB-new';

%%
S = size(files,1);
for s=1:S 
    fprintf('%i ',s)
    
    fname0        = strtrim(files(s,:));
    [fol,nam,ext] = fileparts(fname0);   
    fname1        = fullfile(ndir,[nam ext]);
    
    copyfile(fname0,fname1);
end
fprintf('\n')

%%
files = spm_select('FPList',ndir,'^.*\.nii$');
S     = size(files,1);

for s=1:S
    fprintf('.')
    
    pth           = strtrim(files(s,:));    
    [fol,nam,ext] = fileparts(pth);        

    a.pth        = pth;
    a.name       = nam;
    a.modality   = 'CT';
    a.population = 'ROB';
    a.healthy    = true;
    
    a = orderfields(a);
    
    pth_json = fullfile(ndir,[nam '.json']);
    spm_jsonwrite(pth_json,a);
end
fprintf('\n')