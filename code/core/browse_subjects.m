function browse_subjects(pth)
folder    = dir(pth);
folder    = folder(3:end);
dirflag   = [folder.isdir];
subfolder = folder(dirflag);
S         = numel(subfolder);

dir_scans = fullfile(pth,subfolder(1).name);     

scans = dir(fullfile(dir_scans,'*.nii'));
if isempty(scans)
    scans = dir(fullfile(dir_scans,'*.img'));
end

N = numel(scans);
V = cell(1,S);
for s=1:S
    dir_scans = fullfile(pth,subfolder(s).name);        

    scans = dir(fullfile(dir_scans,'*.nii'));
    if isempty(scans)
        scans = dir(fullfile(dir_scans,'*.img'));
    end

    for n=1:N
        V{s}(n) = spm_vol(fullfile(dir_scans,scans(n).name));
    end              
end 

fname = cell(1,N*S);
cnt   = 1;
for s=1:S
    for n=1:N
        fname{cnt} = V{s}(n).fname;
        cnt        = cnt + 1;
    end
end

spm_check_registration(char(fname));
%==========================================================================