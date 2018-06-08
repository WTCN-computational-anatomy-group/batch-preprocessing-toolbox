function dat = add_segmentations2dat(dat)
% For segmenting we assume... 
m = 1; % ...one modality...
n = 1; % ...and one image per channel       
    
type   = {'c','rc','wc','mwc'};
tissue = {'GM','WM','CSF','BONE','ST','BG'};

T = numel(type);
K = numel(tissue);

dat.segmentation     = {};
dat.segmentation_map = containers.Map;
for t=1:T
    if isfield(dat.modality{m},'channel')    
        C = numel(dat.modality{m}.channel);
        for c=1:C            
            fname = dat.modality{m}.channel{c}.nii(n).dat.fname;        
            
            [pth,nam,ext] = fileparts(fname);
            fname1        = fullfile(pth,[type{t} '1' nam ext]);
            
            if exist(fname1,'file')==2
                break
            end
        end
    else        
        fname = dat.modality{m}.nii(n).dat.fname;   
    end          

    dat.segmentation_map(type{t}) = t;
    dat.segmentation{t}.class     = {};
    dat.segmentation{t}.name      = type{t};
    dat.segmentation{t}.class_map = containers.Map;
    
    for k=1:K        
        [pth,nam,ext] = fileparts(fname);
        fname1        = fullfile(pth,[type{t} num2str(k) nam ext]);
        Nii           = nifti(fname1);
        
        dat.segmentation{t}.class_map(tissue{k}) = k;
        dat.segmentation{t}.class{k}.name        = tissue{k};
        dat.segmentation{t}.class{k}.nii         = Nii;        
                               
        a            = struct;
        a.name       = dat.name;
        a.population = dat.population;
        a.pth        = fname1;
        a.type       = type{t};
        a.tissue     = tissue{k};

        a = orderfields(a);
        
        pth_json = fullfile(pth,[type{t} num2str(k) nam '.json']);
        spm_jsonwrite(pth_json,a);        
        dat.segmentation{t}.class{k}.json.pth = pth_json;
    end
end
%========================================================================== 