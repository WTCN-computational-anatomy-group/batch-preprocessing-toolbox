clear; clc;

dirs = {{'/home/mbrud/Data/images/IXI/IXI-T1','T1'},...
        {'/home/mbrud/Data/images/IXI/IXI-T2','T2'},...
        {'/home/mbrud/Data/images/IXI/IXI-PD','PD'}};%,...
%         {'/home/mbrud/Downloads/IXI-MRA','MRA'}};
ndir = '/data/mbrud/populations/original/IXI/';
    
metadata = xlsread('/home/mbrud/Data/images/IXI/IXI.xls');
col_sub   = 1;
col_age   = 12;
col_sex   = 2;
subjs     = metadata(:,1);

D = numel(dirs);
for d=1:D
    fprintf('%i ',d)
    
    files = spm_select('FPListRec',dirs{d},'^.*\.nii$');    
    
    S = size(files,1);
    for s=1:S
        fprintf('%i ',s)
        
        fname       = strtrim(files(s,:));
        [~,nam,ext] = fileparts(fname);
        nfname      = fullfile(ndir,[nam ext]);
                       
        sub = str2double(nam(4:6));
        ix  = find(subjs==sub);
        
        sex = metadata(ix,col_sex);
        
%         if isempty(sex) || numel(sex)~=1 || ~isfinite(sex)
%             continue        
%         end    
        
        if sex==1
            sex = 'M';
        else
            sex = 'F';
        end
        
        age = metadata(ix,col_age);
        
        if isempty(age) || numel(age)~=1 || ~isfinite(age)
            continue        
        end           
        
        a            = struct;
        a.pth        = nfname;
        a.name       = nam(1:6);
        a.modality   = 'MRI';
        a.channel    = dirs{d}{2};
        a.population = 'IXI';
        a.healthy    = true;
        a.age        = age;
        a.sex        = sex;

        a = orderfields(a);

        pth_json = fullfile(ndir,[nam '.json']);
        spm_jsonwrite(pth_json,a);
                
        copyfile(fname,nfname);        
    end
    fprintf('\n')
end   