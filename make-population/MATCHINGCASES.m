clear; clc;

dir = '/data/mbrud/images-to-be-organised/MatchingCases_T1T2DWIFlair_anon/';
files = spm_select('FPListRec',dir,'^.*\.nii$');    

ndir = '/data/mbrud/populations/original/MatchingCases_T1T2DWIFlair/';
    
%%
S = size(files,1);
for s=1:S
    fprintf('%i ',s)

    fname       = strtrim(files(s,:));
    [pth,nam,ext] = fileparts(fname);

    pth = strsplit(pth,filesep);
    channel = pth{end};
    name    = pth{end - 1};
    
    if strcmp(channel,'DWI'), continue; end
    
    nfname = fullfile(ndir,[name '_' channel ext]);
    copyfile(fname,nfname);

    ix = strfind(nam,'_');
    
    age  = str2double(nam(ix(1) + 1:ix(2) - 1));
    sex  = nam(ix(2) + 1:ix(3) - 1);
    
    a            = struct;
    a.pth        = nfname;
    a.name       = name;
    a.modality   = 'MRI';
    a.channel    = channel;
    a.population = 'MatchingCases_T1T2DWIFlair';
    a.healthy    = true;
    a.age        = age;
    a.sex        = sex;
    
    a = orderfields(a);

    pth_json = fullfile(ndir,[name '_' channel '.json']);
    spm_jsonwrite(pth_json,a);
end
fprintf('\n')