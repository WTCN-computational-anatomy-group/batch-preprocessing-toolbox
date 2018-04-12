function preprocess_images

%--------------------------------------------------------------------------
% OBS! Below parameters need to be set (for FIL users)
%--------------------------------------------------------------------------
pth2_distributed_toolbox = '/cherhome/mbrud/dev/distributed-computing';
pth2_auxiliary_functions = '/cherhome/mbrud/dev/auxiliary-functions';
holly_server_login       = 'mbrud';
holly_client_folder      = '/data/mbrud/tmp-preproc/';
holly_matlab_add_src     = '/home/mbrud/dev/batch-preprocessing-toolbox';
holly_matlab_add_aux     = '/home/mbrud/dev/auxiliary-functions';

% addpath
%--------------------------------------------------------------------------
addpath(genpath('./code'))
addpath(pth2_distributed_toolbox)
addpath(pth2_auxiliary_functions)

%--------------------------------------------------------------------------
% Set distribute package parameters
%--------------------------------------------------------------------------

test_level = 0; % 0: no testing | 1: 1 subject | 2: 8 subjects (parfor) | 3: 8 subjects (holly)

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
holly.job.mem       = '6G';
holly.job.use_dummy = true;

if     test_level==1, holly.server.ip  = ''; holly.client.workers = 0;
elseif test_level==2  holly.server.ip  = ''; holly.client.workers = Inf;
end

holly = distribute_default(holly);

%--------------------------------------------------------------------------
% Set algorithm parameters
%--------------------------------------------------------------------------
m = 0;

% CT data
%-------------
m = m + 1;

pars.dat{m}.dir_data = '/data/mbrud/images/CT/CHROMIS';
pars.dat{m}.modality = 'CT';

pars.dat{m}.preproc.do_realign2mni = true;
pars.dat{m}.preproc.do_crop = true;
pars.dat{m}.preproc.do_rem_neck = true;
pars.dat{m}.preproc.write_2d = true;

pars = pars_default(pars,test_level);

%--------------------------------------------------------------------------
% Start processing images
%--------------------------------------------------------------------------

print_progress('Started');

pars = read_images(pars); 
obj  = init_obj(pars);

obj   = unfold_cell(obj,2);
[~,~] = distribute(holly,'process_subject','inplace',obj);

print_progress('Finished');

m = 1; browse_subjects(obj{m}.dir_preproc_2d,obj{m}.modality);
%==========================================================================

%==========================================================================
function print_progress(status)
date = datestr(now,'mmmm dd, yyyy HH:MM:SS');
fprintf('=======================================================\n')  
fprintf('| %s | %s processing images.\n',date,status)
fprintf('=======================================================\n\n')   
%==========================================================================