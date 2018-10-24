pth_distributed_toolbox = '/data/mbrud/dev/distributed-computing';
pth_auxiliary_functions = '/data/mbrud/dev/auxiliary-functions';

addpath(pth_distributed_toolbox)
addpath(pth_auxiliary_functions)

f = '/data/mbrud/populations/original/MRBrainS18';

files = spm_select('FPList',f,'^.*\.json$');

%%
ofield = 'healthy';
nfield = 'lesion';
nval   = true;

S0 = size(files,1);
for s=1:S0
    fprintf('.');
    
    fname = strtrim(files(s,:));
    
    spm_json_manager('replace_json_field',fname,ofield,nfield,nval)
end
fprintf('\n');
fprintf('Done!\n');
