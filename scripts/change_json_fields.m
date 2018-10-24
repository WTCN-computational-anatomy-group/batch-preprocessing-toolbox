pth_distributed_toolbox = '/data/mbrud/dev/distributed-computing';
pth_auxiliary_functions = '/data/mbrud/dev/auxiliary-functions';

addpath(pth_distributed_toolbox)
addpath(pth_auxiliary_functions)

f = '/home/mbrud/dev/SegModel/populations/2d/MRBrainS18-ra-cr-rn-res';

files = spm_select('FPList',f,'^.*\.json$');

%%
field = 'population';
nval  = 'MRBrainS18';

S0 = size(files,1);
for s=1:S0
    fprintf('.');
    
    fname = strtrim(files(s,:));
    
    spm_json_manager('modify_json_field',fname,field,nval)
end
fprintf('\n');
fprintf('Done!\n');
