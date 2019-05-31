clear; clc;

DirData   = '/home/mbrud/Data/Mo/best_session';
DirOut    = '/data/mbrud/populations/original/MO';
FileNames = '/home/mbrud/Data/Mo/only_heads_T1_T2_FLAIR.mat';
var       = load(FileNames);
heads     = var.heads;

sort_heads = cellstr(heads);
sort_heads = sort(sort_heads);
S          = numel(sort_heads);

for s=1:S
    fprintf('.');
    
    FileName     = sort_heads{s};
    SpltFileName = strsplit(FileName,'\');
    
    cFileName = fullfile(DirData,SpltFileName{1},SpltFileName{2},SpltFileName{3});
    nFileName = fullfile(DirOut,[SpltFileName{1} '_' SpltFileName{2} '_' SpltFileName{3}]);
    
    Name     = SpltFileName{1};
    Sequence = SpltFileName{2};
    copyfile(cFileName,nFileName);
    
    [~,nam] = fileparts(nFileName);
    
    j            = struct;
    j.name       = Name;
    j.population = 'MO';
    j.pth        = [nam '.nii'];                                     
    j.modality   = 'MRI';
    j.channel    = Sequence;
    j.healthy    = true;         

    j = orderfields(j);

    pth_json = fullfile(DirOut,[nam '.json']);
    spm_jsonwrite(pth_json,j);
end
fprintf('\nDone!\n');

spm_json_manager('init_dat',DirOut,fullfile(DirOut,'dat.mat'));