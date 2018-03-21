function nii2subject(dir_in,dir_out)
files = dir(fullfile(dir_in,'*.nii'));
S     = numel(files);

if exist(dir_out,'dir'), rmdir(dir_out,'s'); end; mkdir(dir_out);

for s=1:S
    dir_s = fullfile(dir_out,['S' num2str(s)]);
    mkdir(dir_s);
    
    fname = fullfile(dir_in,files(s).name);
    
    copyfile(fname,dir_s);
end

