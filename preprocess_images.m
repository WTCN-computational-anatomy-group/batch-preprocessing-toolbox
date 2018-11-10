function preprocess_images

%--------------------------------------------------------------------------
% Give path to job definition JSON and required toolboxes
%--------------------------------------------------------------------------

dir_distributed_toolbox = '../distributed-computing';
dir_auxiliary_toolbox   = '../auxiliary-functions';
dir_mtv_preproc         = '../../MTV-preproc';
dir_core_functions      = './code';

%--------------------------------------------------------------------------
% addpath
%--------------------------------------------------------------------------

addpath(genpath(dir_core_functions))
addpath(dir_distributed_toolbox)
addpath(dir_auxiliary_toolbox)
addpath(dir_mtv_preproc)

%--------------------------------------------------------------------------
% pars
%--------------------------------------------------------------------------

TEST_LEVEL = 0; % 0: no testing | 1: 1 subject | 2: 8 subjects (parfor) | 3: 16 subjects (holly)
VX         = [];
RAM        = '8G';
S          = Inf;

%--------------------------------------------------------------------------
% Define list of jobs
%--------------------------------------------------------------------------

jobs          = {};

% Haem template
% jobs{end + 1} = './jobs/default/MRBrainS18-3class.json';
% jobs{end + 1} = './jobs/default/ATLAS.json';
% jobs{end + 1} = './jobs/default/CROMIS.json';
% jobs{end + 1} = './jobs/default/CROMIS-LABELS.json';
% jobs{end + 1} = './jobs/default/IXI.json';
% jobs{end + 1} = './jobs/default/ROB.json';

% MRF-Net
% jobs{end + 1} = './jobs/MRF-Net/MRBrainS18-3class-T1.json';
jobs{end + 1} = './jobs/MRF-Net/OASIS-MICCAI-3class.json'; % PROBLEM HERE!

if TEST_LEVEL == 0 || TEST_LEVEL == 3
    
    %--------------------------------------------------------------------------
    % Prepare code for running on the Holly cluster   
    %--------------------------------------------------------------------------        
    
    % Copy most recent code to holly
    dir_holly = '/data/mbrud/Holly/code/batch-preprocessing-toolbox';
    if exist(dir_holly,'dir'), rmdir(dir_holly,'s'); end; mkdir(dir_holly);   
    copyfile(dir_core_functions,dir_holly);
    
    dir_holly = '/data/mbrud/Holly/code/auxiliary-functions';
    if exist(dir_holly,'dir'), rmdir(dir_holly,'s'); end; mkdir(dir_holly);   
    copyfile(fullfile(dir_auxiliary_toolbox,'*.m'),dir_holly);
    
    dir_holly = '/data/mbrud/Holly/code/mtv-preproc';
    if exist(dir_holly,'dir'), rmdir(dir_holly,'s'); end; mkdir(dir_holly);   
    copyfile(dir_mtv_preproc,dir_holly);
end

%--------------------------------------------------------------------------
% Begin processing
%--------------------------------------------------------------------------

for j=1:numel(jobs)
    
    %----------------------------------------------------------------------
    % Get a job
    %----------------------------------------------------------------------
    
    job = jobs{j};
    
    %----------------------------------------------------------------------
    % Init algorithm 
    %----------------------------------------------------------------------

    % Get job parameters
    [job,holly] = preproc_default(job,TEST_LEVEL);

    % Change some parameters
    if ~isempty(VX),  job.preproc.vx = VX;  end
    if ~isempty(RAM), holly.job.mem  = RAM; end
    if ~isempty(S),   job.S          = S;   end
    
    % Create data object
    dat = spm_json_manager('init_dat',job.dir_population);
    dat = dat(1:min(job.S,numel(dat)));
        
    % Build directory structure
    [dir_preproc,dir_2d] = build_folder_structure(job);

    % Create options struct
    opt = init_opt(job,dir_preproc,dir_2d);

    %----------------------------------------------------------------------
    % Start processing images
    %----------------------------------------------------------------------

    print_progress('Started');
    [~,~] = distribute(holly,'process_subject','inplace',dat,opt);
    
    spm_json_manager('make_pth_relative',dir_preproc,false);

    % Create dat.mat objects (for faster loading of population)
    fname = 'dat.mat';
    if job.write_3d
        cd(dir_preproc);
        spm_json_manager('init_dat','.',fullfile('.',fname));
        
        here = fileparts(mfilename('fullpath'));
        cd(here);
    else
        if exist(dir_preproc,'dir')
            rmdir(dir_preproc,'s');  
        end
    end
    if ~isempty(dir_2d)
        cd(dir_2d);
        spm_json_manager('make_pth_relative','.',false);
        spm_json_manager('init_dat','.',fullfile('.',fname)); 
        
        here = fileparts(mfilename('fullpath'));
        cd(here);
    end        
    
    print_progress('Finished');
end
%==========================================================================

%==========================================================================
function print_progress(status)
date = datestr(now,'mmmm dd, yyyy HH:MM:SS');
fprintf('=======================================================\n')  
fprintf('| %s | %s processing images\n',date,status)
fprintf('=======================================================\n\n')   
%==========================================================================