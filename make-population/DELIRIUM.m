clear; clc;

[files,dirs] = spm_select('FPListRec','/data/mbrud/populations/original/DELIRIUM/imaging-data','^.*\.nii$');
S            = size(files,1);

%%
for s=1:S    
    fname0 = strtrim(files(s,:));
    
    for s1=1:S
       if s==s1, continue; end
       
       fname1 = strtrim(files(s1,:));
       
       if strcmp(fname0,fname1)
           warning('strcmp(fname0,fname1)')
       end
    end
end

%%
ndir = '/data/mbrud/populations/original/temp';

for s=1:S 
    fprintf('%i ',s)
    
    fname0        = strtrim(files(s,:));
    [fol,nam,ext] = fileparts(fname0);   
    fname1        = fullfile(ndir,[nam ext]);
    
    copyfile(fname0,fname1);
end
fprintf('\n')

%%
pth_meta_data = '/data/mbrud/populations/original/DELIRIUM-old/meta-data/deliriumtargets.mat';
meta_data     = load(pth_meta_data);

files = spm_select('FPListRec','/data/mbrud/populations/original/DELIRIUM/','^.*\.nii$');
S     = size(files,1);

ix_all = [];
for s=1:S
    fprintf('.')
    
    pth           = strtrim(files(s,:));    
    [fol,nam,ext] = fileparts(pth);        
   
    ix   = strfind(nam,'_');
    name = nam(1:ix - 1);
    ix   = str2double(nam(1:ix - 1));
    
    ix_all = [ix_all ix];
    
    age = meta_data.imageage{ix};
    ix  = strfind(age,'Y');
    age = str2double(age(1:ix - 1));
    sex = meta_data.imagesex{ix};
            
    a.pth        = pth;
    a.name       = name;
    a.modality   = 'CT';
    a.age        = age;
    a.sex        = sex;
    a.population = 'DELIRIUM';
    a.healthy    = false;
    
    a = orderfields(a);
    
    pth_json = fullfile(fol,[nam '.json']);
    spm_jsonwrite(pth_json,a);
end

if ~(length(ix_all) == length(unique(ix_all)))
    warning('There are replicates...')
end

fprintf('\n')