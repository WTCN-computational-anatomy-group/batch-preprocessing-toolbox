clc; clear;

%% addpath
dir_auxiliary_toolbox   = '../../auxiliary-functions';
dir_mtv_preproc         = '../../../MTV-preproc';
dir_core_functions      = '../code';

addpath(genpath(dir_core_functions))
addpath(dir_auxiliary_toolbox)
addpath(dir_mtv_preproc)

%% Path to images
DirImages = '../data/traning_001';
PthT1     = fullfile(DirImages,'T1.nii');
PthT2     = fullfile(DirImages,'T2.nii');
PthPD     = fullfile(DirImages,'PD.nii');
PthCT     = fullfile(DirImages,'CT.nii');

%% Create dat struct (given as input to the preprocessing code)
dat = struct;
dat.modality{1}.name            = 'MRI';
dat.modality{1}.channel{1}.name = 'T1';
dat.modality{1}.channel{1}.nii  = nifti(PthT1);
dat.modality{1}.channel{2}.name = 'T2';
dat.modality{1}.channel{2}.nii  = nifti(PthT2);
dat.modality{1}.channel{3}.name = 'PD';
dat.modality{1}.channel{3}.nii  = nifti(PthPD);
dat.modality{2}.name            = 'CT';
dat.modality{2}.nii             = nifti(PthCT);

%% Define options
opt             = struct;
opt.dir_preproc = './temp';
opt             = preproc_default(opt);

opt.preproc.reset_origin   = true;
opt.preproc.do_realign2mni = true;
opt.preproc.do_coreg       = true;
opt.preproc.do_crop        = true;
opt.preproc.do_rem_neck    = true;
opt.preproc.do_denoise     = true;
opt.preproc.mc_denoise     = false;
opt.preproc.do_reslice     = true;
opt.preproc.vx             = 2;

dir_preproc = build_folder_structure(opt);    
opt         = init_opt(opt,dir_preproc);
    
%% Run preprocessing
dat = process_subject(dat,opt);