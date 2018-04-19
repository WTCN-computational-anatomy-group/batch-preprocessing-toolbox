function preprocess_images

test_level = 3; % 0: no testing | 1: 1 subject | 2: 8 subjects (parfor) | 3: 8 subjects (holly)

%--------------------------------------------------------------------------
% OBS! Below parameters need to be set (for FIL users)
%--------------------------------------------------------------------------
% pth_distributed_toolbox = '/home/mbrud/dev/distributed-computing';
% pth_auxiliary_functions = '/home/mbrud/dev/auxiliary-functions';
pth_distributed_toolbox = '/cherhome/mbrud/dev/distributed-computing';
pth_auxiliary_functions = '/cherhome/mbrud/dev/auxiliary-functions';

holly_server_login   = 'mbrud';
holly_client_folder  = '/data/mbrud/TEMP/batch-preprocessing-toolbox/';
holly_matlab_add_src = '/home/mbrud/dev/batch-preprocessing-toolbox';
holly_matlab_add_aux = '/home/mbrud/dev/auxiliary-functions';

% addpath
%--------------------------------------------------------------------------
addpath(genpath('./code'))
addpath(pth_distributed_toolbox)
addpath(pth_auxiliary_functions)

%--------------------------------------------------------------------------
% Set distribute package parameters
%--------------------------------------------------------------------------

holly               = struct;
holly.server.ip     = 'holly';
holly.server.login  = holly_server_login;
holly.client.folder = holly_client_folder;
holly.server.folder = holly.client.folder;
holly.matlab.bin    = '/share/apps/matlab';
holly.matlab.addsub = holly_matlab_add_src;
holly.matlab.add    = holly_matlab_add_aux;
holly.restrict      = 'char';
holly.clean         = false;
holly.clean_init    = true;
holly.verbose       = false;
holly.job.mem       = '8G';
holly.job.use_dummy = true;

if     test_level==1, holly.server.ip  = ''; holly.client.workers = 0;
elseif test_level==2  holly.server.ip  = ''; holly.client.workers = Inf;
end
% holly.server.ip  = ''; holly.client.workers = 0;
% holly.server.ip  = ''; holly.client.workers = Inf;
    
holly = distribute_default(holly);

%--------------------------------------------------------------------------
% Set algorithm parameters
%--------------------------------------------------------------------------
pars = '/home/mbrud/Dropbox/PhD/Data/pars/batch-preprocessing-toolbox/denoise-CT-3d.json';

pars = preproc_default(pars,test_level);

%--------------------------------------------------------------------------
% Start processing images
%--------------------------------------------------------------------------

print_progress('Started');

pars = read_images(pars); 
obj  = init_obj(pars);

[obj,ix] = unfold_cell(obj,2);
[~,~]    = distribute(holly,'process_subject','inplace',obj);
obj      = fold_cell(obj,ix);

print_progress('Finished');

% Show the preprocessing results
m = 1; 
browse_subjects(obj{m}{1}.dir_preproc_2d,obj{m}{1}.modality);
%==========================================================================

%==========================================================================
function print_progress(status)
date = datestr(now,'mmmm dd, yyyy HH:MM:SS');
fprintf('=======================================================\n')  
fprintf('| %s | %s processing images\n',date,status)
fprintf('=======================================================\n\n')   
%==========================================================================