clear; close all; clc;

addpath(genpath('../code'))

modality = 'CT';

dir_dicom = '/home/mbrud/Data/DCM/CROMIS-full/';
dir_nii   = '/home/mbrud/Data/DCM/CROMIS-nii/';
dir_final = '/data/mbrud/images/CT/CROMIS/';

%% Convert DCMs
%--------------------------------------------------------------------------

search_and_convert_dcm(dir_dicom,dir_nii);

%% Copy NIfTIs so to have correct directory structure
%--------------------------------------------------------------------------

if exist(dir_final,'dir')
    rmdir(dir_final,'s');
end
mkdir(dir_final);

files = dir(fullfile(dir_nii,'*.nii'));
S     = numel(files);

for s=1:S
    fprintf('.')
    
    fname = fullfile(dir_nii,files(s).name);
    
    dir_subj = fullfile(dir_final,['S' num2str(s)]);    
    mkdir(dir_subj);

    dir_scans = fullfile(dir_subj,'scans');    
    mkdir(dir_scans);
    
    dir_modal = fullfile(dir_scans,modality);    
    mkdir(dir_modal);
    
    copyfile(fname,dir_modal);
end
fprintf('\n')