function pth_nii = search_and_convert_dcm(dir_dcm,dir_nii)
% Searches folder structure for sets of DICOM files and converts these to NIFTI format.
% FORMAT pth_nii = search_and_convert_dcm(dir_dcm,dir_nii)

seen = {}; % cell array to store seen edges
dcm_DFS(dir_dcm,[],seen,dir_nii);

pth_nii = dir(fullfile(dir_nii,'*.nii'));
pth_nii = {pth_nii(:).name};
pth_nii = cellfun(@(c) fullfile(dir_nii,c),pth_nii,'uni',false);
pth_nii = pth_nii';
%==========================================================================

%==========================================================================
function seen = dcm_DFS(R,v,seen,niidir)
% Recursive depth-first search (https://en.wikipedia.org/wiki/Depth-first_search)

if isempty(v)
    v = R;
end

% Look for DICOMs
P = [dir(fullfile(v,'*.dcm')); dir(fullfile(v,'*.DCM'))];
if ~isempty(P)
    % There are DICOM files in the folder
    P = {P(:).name};
    P = cellfun(@(c) fullfile(v,c),P,'uni',false);

    % Read DICOMs
    H   = spm_dicom_headers(char(P));
    my_spm_dicom_convert(H,'all','flat','nii',niidir);
end 

seen{end + 1,1} = v;                        % label v as discovered

w = dirs_in_dir(v);                         % get adjacents edges to v

w(ismember(w,v)) = [];                      % remove edge v

for i=1:numel(w)
    if ~ismember(seen,w{i})                 % if w not labeled as discovered
        seen = dcm_DFS(R,w{i},seen,niidir); % call dir_DFS recursively
    end
end
%==========================================================================

%==========================================================================
function dirs = dirs_in_dir(P)
% Get all subfolders of folder

d    = dir(P);
isub = [d(:).isdir];
dirs = {d(isub).name}';

dirs(ismember(dirs,{'.','..'})) = [];

for i=1:numel(dirs)
   dirs{i} = fullfile(P,dirs{i});
end
%==========================================================================