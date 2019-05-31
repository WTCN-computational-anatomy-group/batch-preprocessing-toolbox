clear;

DirData = '/home/mbrud/Data/labels/Balgrist_controls/BALGRIST-T1';

Filenames = spm_select('FPList',DirData,'^.*\.img$');

%% Convert to NIfTI

for i=1:size(Filenames,1)
    Filename    = strtrim(Filenames(i,:));
    [~,nam,ext] = fileparts(Filename);
    
    V       = spm_vol(Filename);
    ima     = spm_read_vols(V);
    V.fname = fullfile(DirData,[nam '.nii']);
    spm_write_vol(V,ima);
end

%%

Subjects = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 21 22 23 24];

for s=1:numel(Subjects)
    if Subjects(s) < 10
        prefix = ['c0' num2str(Subjects(s)) '-'];
    else
        prefix = ['c' num2str(Subjects(s)) '-'];
    end
   
    SubjectFiles = spm_select('FPList',DirData,['^' prefix '.*\.nii$']);
    
    name = ['c0' num2str(Subjects(s))];
    
    for i=1:2
        Filename    = strtrim(SubjectFiles(i,:));
        [~,nam,ext] = fileparts(Filename);
        
        if contains(nam,'spine')
            % Labels
            a              = struct;
            a.name         = name;     
            a.rater        = 'unknown';
            a.population   = 'BALGRIST-T1';
            a.nam_modality = 'MRI';
            a.nam_channel  = 'T1';
            a.ix_img       = 1;
            a.pth          = [nam '.nii'];
            
            pth_json = fullfile(DirData,[nam '.json']);
            spm_jsonwrite(pth_json,a);
        else
            % T1
            a            = struct;
            a.name       = name;
            a.population = 'BALGRIST-T1';
            a.modality   = 'MRI';
            a.channel    = 'T1';
            a.healthy    = true;
            a.pth        = [nam '.nii'];

            a = orderfields(a);

            pth_json = fullfile(DirData,[nam '.json']);
            spm_jsonwrite(pth_json,a);
        end        
    end
end