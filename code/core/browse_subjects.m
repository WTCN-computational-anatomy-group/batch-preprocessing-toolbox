function browse_subjects(pth,modality)
folder    = dir(pth);
folder    = folder(3:end);
dirflag   = [folder.isdir];
subfolder = folder(dirflag);
S         = numel(subfolder);

fname = {};
cnt   = 1;
for s=1:S
    dir_scans = fullfile(pth,subfolder(s).name,'scans',modality);        

    scans = dir(fullfile(dir_scans,'*.nii'));
    if isempty(scans)
        scans = dir(fullfile(dir_scans,'*.img'));
    end

    if ~isempty(scans)
        fname{cnt} = fullfile(dir_scans,scans(1).name);       
        cnt        = cnt + 1;
    end
end 

spm_check_registration(char(fname));
%==========================================================================