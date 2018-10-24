clear; clc;

[files,dirs] = spm_select('FPListRec','/data/mbrud/to-be-deleted/CROMIS-old','^.*\.nii$');
ndir         = '/data/mbrud/populations/original/CROMIS';

%%
S = size(files,1);
for s=1:S 
    fprintf('%i ',s)
    
    fname0        = strtrim(files(s,:));
    [fol,nam,ext] = fileparts(fname0);   
    fname1        = fullfile(ndir,[nam ext]);
    
    copyfile(fname0,fname1);
end
fprintf('\n')

%%
files = spm_select('FPList',ndir,'^.*\.nii$');
S     = size(files,1);

for s=1:S
    fprintf('.')
    
    pth           = strtrim(files(s,:));    
    [fol,nam,ext] = fileparts(pth);        

    ix   = strfind(nam,'_');
    name = nam(ix + 1:end);
    
    a.pth        = pth;
    a.name       = name;
    a.modality   = 'CT';
    a.population = 'CROMIS';
    a.lesion     = true;
    
    a = orderfields(a);
    
    pth_json = fullfile(ndir,[nam '.json']);
    spm_jsonwrite(pth_json,a);
end
fprintf('\n')

%%
ndir  = '/data/mbrud/populations/original/CROMIS';
files = spm_select('FPList',ndir,'^.*\.nii$');

cnt = 0;
for s=1:size(files,1)   
    fname0 = strtrim(files(s,:));

    [fol0,nam0,ext] = fileparts(fname0);        
%     ix            = strfind(nam0,'_');
%     name0         = nam0(ix + 1:end);

    for s1=1:size(files,1)   
        if s==s1, continue; end

        fname1 = strtrim(files(s1,:));

        [fol1,nam1,ext] = fileparts(fname1);        
%         ix            = strfind(nam1,'_');
%         name1         = nam1(ix + 1:end);

        if strfind(nam0,nam1)
%             fname1
            delete(fname1);
            json = fullfile(fol1,[nam1 '.json']);
            delete(json);
            cnt = cnt + 1;
%             if strcmp(nam0(1:3),'iso')
%                 delete(fname1);
%                 delete(fullfile(fol1,[nam '.json']))
%             elseif strcmp(nam1(1:3),'iso')
%                 delete(fname0);
%                 delete(fullfile(fol0,[nam '.json']))
%             end
        end
    end
end

cnt