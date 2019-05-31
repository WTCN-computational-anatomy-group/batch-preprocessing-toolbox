clear; clc;

addpath('/home/mbrud/dev/mbrud/code/matlab/auxiliary-functions/');

DirData   = '/scratch/Projects/Mo/OASIS/';
DirOut    = DirData;

FileNames = spm_select('FPList',DirData,'^.*\.nii$');
NumSubj   = size(FileNames,1);

%%
for s=1:NumSubj
    fprintf('.');
    
    FileName = deblank(FileNames(s,:));
    
    [~,FileName] = fileparts(FileName);
    FileParts    = strsplit(FileName,'_');
    Name         = FileParts{1};
    Channel      = FileParts{2};
            
    j            = struct;
    j.name       = Name;
    j.population = 'MO_OASIS';
    j.pth        = [FileName '.nii'];                                     
    j.modality   = 'MRI';
    j.channel    = Channel;
    j.healthy    = true;         

    j = orderfields(j);

    pth_json = fullfile(DirOut,[FileName '.json']);
    spm_jsonwrite(pth_json,j);    
end
fprintf('\nDone!\n');

%%

spm_json_manager('init_dat',DirOut,fullfile(DirOut,'dat.mat'));