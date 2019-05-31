clear; clc; close all;

DirData = '/home/mbrud/Data/challenges/ABCD-default/fmriresults01/image03/testing2'; % Path to ziped image files
DirOut  = '/home/mbrud/Data/challenges/ABCD-testing2';

if exist(DirOut,'dir') == 7, rmdir(DirOut,'s'); end; mkdir(DirOut); 

% Read .tgz files
FileNames = spm_select('FPList',DirData,'^.*\.tgz$');
NumSubj   = size(FileNames,1);

for s=1:NumSubj
    fprintf('%i ',s)
    
    dn = strtrim(FileNames(s,:));
    dn = untar(dn,DirOut);
    dn = dn{1};
    
    name = strsplit(dn,filesep);
    name = name{end - 1};
        
    FileName = spm_select('FPListRec',dn,'^.*\gz$');
    
    for i=1:size(FileName,1)
        fn   = strtrim(FileName(i,:));
        
        if contains(fn,'brain')
            fnuz = gunzip(fn,DirOut);

            [pth,nam,ext] = fileparts(fnuz{1});
            nfn           = [name '_' nam];
            fn            = fullfile(pth,[nfn ext]);
            movefile(fnuz{1},fn);

            j            = struct;
            j.name       = name;
            j.population = 'ABCD';
            j.pth        = [nfn '.nii'];                                     
            j.modality   = 'MRI';
            j.channel    = 'T1';
            j.healthy    = true;         

            j = orderfields(j);

            pth_json = fullfile(DirOut,[nfn '.json']);
            spm_jsonwrite(pth_json,j);
        end
    end
    
    if exist(dn,'dir') == 7, rmdir(dn,'s'); end
end
fprintf('\n')

spm_json_manager('init_dat',DirOut,fullfile(DirOut,'dat.mat'));