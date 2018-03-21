function preprocess_images

addpath(genpath('./code'))
addpath('/cherhome/mbrud/dev/distributed-computing')
addpath('/cherhome/mbrud/dev/auxiliary-functions')

test_level = 0; % 0: no testing | 1: 1 subject | 2: 8 subjects (parfor) | 3: 8 subjects (holly)
dir_tmp    = '/data/mbrud/data-seg';
m          = 0;

%--------------------------------------------------------------------------
% Set distribute package parameters
%--------------------------------------------------------------------------

holly               = struct;
holly.server.ip     = 'holly';
holly.server.login  = 'mbrud';
holly.client.folder = fullfile(dir_tmp,'cluster');
holly.server.folder = holly.client.folder;
holly.matlab.bin    = '/share/apps/matlab';
holly.matlab.addsub = '/home/mbrud/dev/batch-preprocessing-toolbox';
holly.matlab.add    = '/home/mbrud/dev/auxiliary-functions';
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

S = Inf;

% New data-set (healthy)
%--------------------------
m = m + 1;

pars.dat{m}.dir_data = '/data/mbrud/images/CT/healthy/';
pars.dat{m}.modality = 'CT';

pars.dat{m}.S = S;

pars.dat{m}.preproc.do_crop = true;
pars.dat{m}.preproc.do_rem_neck = true;
pars.dat{m}.preproc.do_realign2mni = true;

pars.dat{m}.preproc.do_superres = true;
pars.dat{m}.preproc.superres.trunc_ct = [0 300];
pars.dat{m}.preproc.superres.proj_mat = 'smo';

pars.dat{m}.preproc.do_skull_strip = true; 

pars.dat{m}.preproc.write_2d = true;

% New data-set (CHROMIS)
%--------------------------
m = m + 1;

pars.dat{m}.dir_data = '/data/mbrud/images/CT/CHROMIS/';
pars.dat{m}.modality = 'CT';

pars.dat{m}.S = S;

pars.dat{m}.preproc.do_crop = true;
pars.dat{m}.preproc.do_rem_neck = true;
pars.dat{m}.preproc.do_realign2mni = true;

pars.dat{m}.preproc.do_superres = true;
pars.dat{m}.preproc.superres.trunc_ct = [0 300];
pars.dat{m}.preproc.superres.proj_mat = 'smo';

pars.dat{m}.preproc.do_skull_strip = true; 

pars.dat{m}.preproc.write_2d = true;

% New data-set
% --------------------------
m = m + 1;

pars.dat{m}.dir_data = '/data/mbrud/images/MRI/OASIS-longitudinal/';

pars.dat{m}.S = S;

pars.dat{m}.preproc.do_crop = true;
pars.dat{m}.preproc.do_rem_neck = true;
pars.dat{m}.preproc.do_realign2mni = true;

pars.dat{m}.preproc.write_2d = true;

% New data-set
% --------------------------
m = m + 1;

pars.dat{m}.dir_data = '/data/mbrud/images/MRI/IXI-T1T2PD/';

pars.dat{m}.S = S;

pars.dat{m}.preproc.do_crop = true;
pars.dat{m}.preproc.do_rem_neck = true;
pars.dat{m}.preproc.do_realign2mni = true;

pars.dat{m}.preproc.do_coreg = true;    
pars.dat{m}.preproc.do_reslice = true;

pars.dat{m}.preproc.write_2d = true;

% New data-set
% --------------------------
m = m + 1;

pars.dat{m}.dir_data = '/data/mbrud/images/MRI/IXI-T1/';

pars.dat{m}.S = S;

pars.dat{m}.preproc.do_crop = true;
pars.dat{m}.preproc.do_rem_neck = true;
pars.dat{m}.preproc.do_realign2mni = true;

pars.dat{m}.preproc.write_2d = true;

% New data-set
% --------------------------
m = m + 1;

pars.dat{m}.dir_data = '/data/mbrud/images/MRI/IXI-T1T2PD/';

pars.dat{m}.S = S;

pars.dat{m}.preproc.do_crop = true;
pars.dat{m}.preproc.do_realign2mni = true;

pars.dat{m}.preproc.do_coreg = true;    
pars.dat{m}.preproc.do_reslice = true;

pars.dat{m}.preproc.write_2d = true;

% New data-set
% --------------------------
m = m + 1;

pars.dat{m}.dir_data = '/data/mbrud/images/MRI/IXI-T1/';

pars.dat{m}.S = S;

pars.dat{m}.preproc.do_crop = true;
pars.dat{m}.preproc.do_realign2mni = true;

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

obj     = unfold_cell(obj,2);
[~,obj] = distribute(holly,'process_subject','inplace',obj);

print_progress('Finished');
%==========================================================================

%==========================================================================
function print_progress(status)
date = datestr(now,'mmmm dd, yyyy HH:MM:SS');
fprintf('=================================================\n')  
fprintf('| %s | %s processing images.\n',date,status)
fprintf('=================================================\n\n')   
%==========================================================================