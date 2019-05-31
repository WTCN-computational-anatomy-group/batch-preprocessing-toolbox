clear;

f_img = '/home/mbrud/Data/labels/IBSR20/20Normals_T1_8bit';
f_seg = '/home/mbrud/Data/labels/IBSR20/20Normals_T1_seg';
nf    = '/data/mbrud/populations/original/ISBR20';
mkdir(nf);

files_img = spm_select('FPList',f_img,'^.*\.img$');
files_seg = spm_select('FPList',f_seg,'^.*\.img$');
S     = size(files_img,1);

%%
for s=1:S
    fname       = strtrim(files_img(s,:));        
    [~,nam,ext] = fileparts(fname);
    
    name = strsplit(nam,'_');
    name = name{1};
    
    V=spm_vol(fname);
    ima=spm_read_vols(V);
    V.fname=fullfile(nf,[nam '.nii']);
    spm_write_vol(V,ima);

%     Nii=nifti(fname);
%     Nio=Nii;
%     Nio.dat.fname = fullfile(nf,[nam '.nii']);
%     Nio.dat.offset = 352;
%     Nio.dat(:,:,:,:,:) = Nii.dat(:,:,:,:,:);
        
    a            = struct;
    a.name       = name;
    a.population = 'ISBR20';
    a.modality   = 'MRI';
    a.channel    = 'T1';
    a.healthy    = true;
    a.pth        = [nam '.nii'];

    a = orderfields(a);

    pth_json = fullfile(nf,[name '.json']);
    spm_jsonwrite(pth_json,a);
    
    % Labels
    fname       = strtrim(files_seg(s,:));        
    [~,nam,ext] = fileparts(fname);
    
    name = strsplit(nam,'_');
    name = name{1};
    
    V=spm_vol(fname);
    ima=spm_read_vols(V);
    V.fname=fullfile(nf,[nam '_seg.nii']);
    spm_write_vol(V,ima);
    
%     Nii=nifti(fname);
%     Nio=Nii;
%     Nio.dat.fname = fullfile(nf,[nam '_seg.nii']);
%     Nio.dat.offset = 352;
%     Nio.dat(:,:,:,:,:) = Nii.dat(:,:,:,:,:);
    
    a              = struct;
    a.name         = name;
    a.rater        = 'unknown';
    a.population   = 'ISBR20';
    a.nam_modality = 'MRI';
    a.nam_channel  = 'T1';
    a.ix_img       = 1;
    a.pth          = [nam '_seg.nii'];

    a = orderfields(a);

    pth_json = fullfile(nf,[name '-seg' '.json']);
    spm_jsonwrite(pth_json,a);
end