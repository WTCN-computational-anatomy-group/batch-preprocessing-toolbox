clear; clc;

addpath('/data/mbrud/dev/auxiliary-functions/');

f  = '/data/mbrud/images-to-be-organised/OASIS/Longitudinal-MRI';
nf = '/data/mbrud/populations/original/OASIS-LONG';


[~,dirs] = spm_select('FPList',f);
D        = size(dirs,1);

read = cell(1,D);
for i1=1:D    
    d1 = strtrim(dirs(i1,:));
    
    name = strsplit(d1,filesep);
    name = name{end}(1:9);    
    
    if any(strcmp(read,name))
        read{i1} = name;
        continue;
    else
        read{i1} = name;
    end
     
    d2    = fullfile(d1,'RAW');
    files = spm_select('FPList',d2,'^mpr-1.*\.img');   
    fname = deblank(files(1,:));
    
    [~,nam] = fileparts(fname);
    nfname  = fullfile(nf,[name '-' nam '.nii']);
    
    V       = spm_vol(fname);
    ima     = spm_read_vols(V);
    V.fname = nfname;
    spm_write_vol(V,ima);
        
    a            = struct;
    a.pth        = [name '-' nam '.nii'];
    a.name       = name;
    a.modality   = 'MRI';
    a.channel    = 'T1';
    a.population = 'OASIS-LONG';
    a.healthy    = false;

    a = orderfields(a);

    pth_json = fullfile(nf,[name '-' nam '.json']);
    spm_jsonwrite(pth_json,a);
end