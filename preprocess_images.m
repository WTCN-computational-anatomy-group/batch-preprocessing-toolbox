function dat = preprocess_images

%--------------------------------------------------------------------------
% Give path to job definition JSON and required toolboxes
%--------------------------------------------------------------------------

dir_distributed_toolbox = '../distributed-computing';
dir_auxiliary_toolbox   = '../auxiliary-functions';
dir_mtv_preproc         = '../MTV-preproc';
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

S   = 50;
VX  = [];
RAM = '';

%--------------------------------------------------------------------------
% Define list of jobs
%--------------------------------------------------------------------------

jobs = {};

% T1 7 labels (spine)
% jobs{end + 1} = './jobs/t1-lab7-spine/MRBrainS18.json';
% jobs{end + 1} = './jobs/t1-lab7-spine/OASIS-MICCAI.json';
% jobs{end + 1} = './jobs/t1-lab7-spine/IXI.json';
% jobs{end + 1} = './jobs/t1-lab7-spine/ATLAS-LABELS.json';
% jobs{end + 1} = './jobs/t1-lab7-spine/BALGRIST-T1.json';

% % T1 7 labels
% jobs{end + 1} = './jobs/t1-lab7/MRBrainS18.json';
% jobs{end + 1} = './jobs/t1-lab7/OASIS-MICCAI.json';
% jobs{end + 1} = './jobs/t1-lab7/IXI.json';
% jobs{end + 1} = './jobs/t1-lab7/ATLAS-LABELS.json';
% jobs{end + 1} = './jobs/t1-lab7/BALGRIST-T1.json';

% % T1 3 labels
% jobs{end + 1} = './jobs/t1-lab3/MRBrainS18.json';
% jobs{end + 1} = './jobs/t1-lab3/OASIS-MICCAI.json';
% jobs{end + 1} = './jobs/t1-lab3/IXI.json';
% jobs{end + 1} = './jobs/t1-lab3/ATLAS-LABELS.json';

% MO
% jobs{end + 1} = './jobs/MO/MO.json';
% jobs{end + 1} = './jobs/MO/MO-OASIS.json';

% ABCD
% jobs{end + 1} = './jobs/abcd/ABCD.json';

% Spine-11
% jobs{end + 1} = './jobs/spine-11/OASIS-MICCAI.json';
% jobs{end + 1} = './jobs/spine-11/MRBrainS18.json';

% Spine-3
% jobs{end + 1} = './jobs/spine-3/IXI.json';
jobs{end + 1} = './jobs/spine-3/MRBrainS18.json';
% jobs{end + 1} = './jobs/spine-3/OASIS-MICCAI.json';
% jobs{end + 1} = './jobs/spine-3/BALGRIST-T1.json';
% jobs{end + 1} = './jobs/spine-3/DELIRIUM.json';
 
% Haem template
% jobs{end + 1} = './jobs/haem/ATLAS-LABELS.json';
% jobs{end + 1} = './jobs/haem/CROMIS-LABELS.json';
% jobs{end + 1} = './jobs/haem/MRBrainS18.json';
% jobs{end + 1} = './jobs/haem/IXI.json';
% jobs{end + 1} = './jobs/haem/OASIS-MICCAI.json';
% jobs{end + 1} = './jobs/haem/DELIRIUM.json';
% jobs{end + 1} = './jobs/haem/BALGRIST-T1.json';

% MTV
% jobs{end + 1} = './jobs/MTV/MatchingCases_T1T2DWIFlair.json';
% jobs{end + 1} = './jobs/MTV/DELIRIUM.json';
% jobs{end + 1} = './jobs/MTV/IXI.json';

% MRF-Net
% jobs{end + 1} = './jobs/MRF-Net/MRBrainS18-3class-T1.json';
% jobs{end + 1} = './jobs/MRF-Net/OASIS-MICCAI-3class.json';
% jobs{end + 1} = './jobs/MRF-Net/ISBR20.json';

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
    
    % Build directory structure
    [dir_preproc,dir_2d] = build_folder_structure(job);

    % Map directories to options struct
    opt = init_opt(job,dir_preproc,dir_2d);
    
    % Create data object
    dat = spm_json_manager('init_dat',job.dir_population,{},'',job.channel');
    dat = dat(1:min(job.S,numel(dat)));            

    %----------------------------------------------------------------------
    % Start processing images
    %----------------------------------------------------------------------

    print_progress('Started');
    [~,dat] = distribute(holly,'process_subject','inplace',dat,opt);
    
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

if strcmpi(job.holly.mode,'qsub')
    % Remove deleted files ending up in holly's hidden trash folders
    spm_misc('clean_holly_mbrud');
end
%==========================================================================

%==========================================================================
function print_progress(status)
date = datestr(now,'mmmm dd, yyyy HH:MM:SS');
fprintf('=======================================================\n')  
fprintf('| %s | %s processing images\n',date,status)
fprintf('=======================================================\n\n')   
%==========================================================================