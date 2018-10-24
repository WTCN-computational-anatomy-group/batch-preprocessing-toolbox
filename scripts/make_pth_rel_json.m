pth_distributed_toolbox = '/data/mbrud/dev/distributed-computing';
pth_auxiliary_functions = '/data/mbrud/dev/auxiliary-functions';

addpath(pth_distributed_toolbox)
addpath(pth_auxiliary_functions)

f          = {};
f{end + 1} = '/data/mbrud/populations/original/ATLAS';
f{end + 1} = '/data/mbrud/populations/original/CROMIS';
f{end + 1} = '/data/mbrud/populations/original/CROMIS-LABELS';
f{end + 1} = '/data/mbrud/populations/original/DELIRIUM';
f{end + 1} = '/data/mbrud/populations/original/IXI';
f{end + 1} = '/data/mbrud/populations/original/MatchingCases_T1T2DWIFlair';
f{end + 1} = '/data/mbrud/populations/original/MRBrainS18';
f{end + 1} = '/data/mbrud/populations/original/OASIS-LONG';
f{end + 1} = '/data/mbrud/populations/original/OASIS-MICCAI';
f{end + 1} = '/data/mbrud/populations/original/ROB';

for i=1:numel(f)
    spm_json_manager('make_pth_relative',f{i});
    dat = spm_json_manager('init_dat',f{i},fullfile(f{i},'dat.mat'));
end