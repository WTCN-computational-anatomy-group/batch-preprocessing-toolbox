clear; clc;

[~,dirs] = spm_select('FPListRec','/home/mbrud/Data/challenges/MRBrainS18/organised','^.*\.nii.gz$');
ndir     = '/data/mbrud/populations/original/MSBrainS18';

%%
D = size(dirs,1);
for i=1:D         
    dir_name = strtrim(dirs(i,:));    
    files    = spm_select('FPList',dir_name,'^.*\.nii.gz$');
    
    S = size(files,1);
    for s=1:S
        pth           = strtrim(files(s,:));
        [fol,nam,ext] = fileparts(pth);   
        fname         = fullfile(ndir,[num2str(i) '-' nam ext]);            
        
        nam  = strsplit(nam,'.');
        nam0 = nam{1};
        
        copyfile(pth,fname);
        
        fname1        = gunzip(fname);
        [fol,nam,ext] = fileparts(fname1{1});   
        delete(fname);
        
        a            = struct;
        a.pth        = fname1{1};       
        a.name       = num2str(i);
        a.population = 'MSBrainS18';        
        if strcmp(nam0,'segm')                                    
            a.rater        = 'unknown';              
            a.nam_modality = 'MRI';
            a.nam_channel  = 'FLAIR';
            a.ix_img       = 1;
        else            
            a.modality   = 'MRI';
            a.channel    = nam0;
            a.healthy    = false;
        end

        a = orderfields(a);
    
        pth_json = fullfile(ndir,[num2str(i) '-' nam0 '.json']);
        spm_jsonwrite(pth_json,a);                    
    end
end
fprintf('\n')