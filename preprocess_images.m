function preprocess_images

%--------------------------------------------------------------------------
% Give path to job definition JSON and required toolboxes
%--------------------------------------------------------------------------

job        = '/data/mbrud/jobs/batch-preprocessing-toolbox/ROB.json';
test_level = 1; % 0: no testing | 1: 1 subject | 2: 8 subjects (parfor) | 3: 16 subjects (holly)

pth_distributed_toolbox = '/data/mbrud/dev/distributed-computing';
pth_auxiliary_functions = '/data/mbrud/dev/auxiliary-functions';

% addpath
%--------------------------------------------------------------------------
addpath(genpath('./code'))
addpath(pth_distributed_toolbox)
addpath(pth_auxiliary_functions)
 
%--------------------------------------------------------------------------
% Init algorithm 
%--------------------------------------------------------------------------

% Get job parameters
[job,holly] = preproc_default(job,test_level);

% Create data object
dat = spm_json_manager('init_dat',job.dir_population);
dat = dat(1:min(job.S,numel(dat)));

% Build directory structure
[dir_preproc,dir_2d] = build_folder_structure(job);

% Create options struct
opt = init_opt(job,dir_preproc,dir_2d);

%--------------------------------------------------------------------------
% Start processing images
%--------------------------------------------------------------------------

print_progress('Started');
[~,~] = distribute(holly,'process_subject','inplace',dat,opt);
print_progress('Finished');
%==========================================================================

%==========================================================================
function print_progress(status)
date = datestr(now,'mmmm dd, yyyy HH:MM:SS');
fprintf('=======================================================\n')  
fprintf('| %s | %s processing images\n',date,status)
fprintf('=======================================================\n\n')   
%==========================================================================