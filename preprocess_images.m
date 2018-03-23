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

test_level = 2; % 0: no testing | 1: 1 subject | 2: 8 subjects (parfor) | 3: 8 subjects (holly)

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

holly = distribute_default(holly);

%--------------------------------------------------------------------------
% Set algorithm parameters
%--------------------------------------------------------------------------

m = 0; % data-set counter

% Example 1 (T1-weighted data)
% -Rigidly realign to MNI space
% -Remove data outside of head (+neck)
% -Bias-field correct
% -Skull strip
% -Make ML-labels
% -Normalise intensities
% -Write 2D versions
% --------------------------
m = m + 1;

pars.dat{m}.dir_data = '/mnt/cifs_share/share_data/Testing/batch-preprocessing-toolbox/MRI-T1'; % In Ashburner_group shared
pars.dat{m}.dir_preproc = '/home/mbrud/Desktop/TEMP/T1'; % If not specified, created automatically

pars.dat{m}.preproc.do_realign2mni = true;
pars.dat{m}.preproc.do_crop = true;
pars.dat{m}.preproc.do_rem_neck = true;
pars.dat{m}.preproc.do_bf_correct = true;      
pars.dat{m}.preproc.do_skull_strip = true;
pars.dat{m}.preproc.make_ml_labels = true;
pars.dat{m}.preproc.normalise_intensities = true;
pars.dat{m}.preproc.write_2d = true;

% Example 2 (T1, T2-, PD-weighted data)
% -Rigidly realign to MNI space
% -Remove data outside of head
% -Co-register
% -Reslice to image with largest FOV
% -Make 1 mm isotropic voxels
% -Write 2D versions
% --------------------------
m = m + 1;

pars.dat{m}.dir_data = '/mnt/cifs_share/share_data/Testing/batch-preprocessing-toolbox/MRI-T1T2PD'; % In Ashburner_group shared
pars.dat{m}.dir_preproc = '/home/mbrud/Desktop/TEMP/T1T2PD'; % If not specified, created automatically

pars.dat{m}.preproc.do_realign2mni = true;
pars.dat{m}.preproc.do_crop = true;
pars.dat{m}.preproc.do_coreg = true;
pars.dat{m}.preproc.do_reslice = true;
pars.dat{m}.preproc.vx = [1 1 1];
pars.dat{m}.preproc.write_2d = true;

% Example 3 (T1-weighted data)
% -Rigidly realign to MNI space
% -Remove data outside of head (+neck)
% -Segment
% -Write 2D versions
% --------------------------
m = m + 1;

pars.dat{m}.dir_data = '/mnt/cifs_share/share_data/Testing/batch-preprocessing-toolbox/MRI-T1'; % In Ashburner_group shared
pars.dat{m}.dir_preproc = '/home/mbrud/Desktop/TEMP/T1-seg'; % If not specified, created automatically

pars.dat{m}.preproc.do_realign2mni = true;
pars.dat{m}.preproc.do_crop = true;
pars.dat{m}.preproc.do_rem_neck = true;
pars.dat{m}.preproc.write_tc = [true(6,1) false(6,1) false(6,1) false(6,1)]; 
pars.dat{m}.preproc.write_2d = true;

% Example 4 (CT data)
% -Rigidly realign to MNI space
% -Remove data outside of head (+neck)
% -Skull strip
% -Write 2D versions
% --------------------------
m = m + 1;

pars.dat{m}.dir_data = '/mnt/cifs_share/share_data/Testing/batch-preprocessing-toolbox/CT'; % In Ashburner_group shared
pars.dat{m}.modality = 'CT';
pars.dat{m}.dir_preproc = '/home/mbrud/Desktop/TEMP/CT'; % If not specified, created automatically

pars.dat{m}.preproc.do_realign2mni = true;
pars.dat{m}.preproc.do_crop = true;
pars.dat{m}.preproc.do_rem_neck = true;
pars.dat{m}.preproc.do_skull_strip = true;
pars.dat{m}.preproc.write_2d = true;

% Finalise parameter struct
%--------------------------
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
%==========================================================================

%==========================================================================
function print_progress(status)
date = datestr(now,'mmmm dd, yyyy HH:MM:SS');
fprintf('=======================================================\n')  
fprintf('| %s | %s processing images.\n',date,status)
fprintf('=======================================================\n\n')   
%==========================================================================