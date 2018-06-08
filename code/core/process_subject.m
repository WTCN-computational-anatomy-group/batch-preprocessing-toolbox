function dat = process_subject(dat,opt) 

% Get options
%--------------------------------------------------------------------------
dir_preproc = opt.dir_preproc;
dir_2d      = opt.dir_2d;
axis_2d     = opt.axis_2d;
write_2d    = opt.write_2d;
preproc     = opt.preproc;
denoise     = preproc.denoise;
superres    = preproc.superres;
M           = numel(dat.modality);

% Make copies into dir_preproc (in order to not modify the original data)
%--------------------------------------------------------------------------
for m=1:M
    if isfield(dat.modality{m},'channel')
        C = numel(dat.modality{m}.channel);
        for c=1:C
            N = numel(dat.modality{m}.channel{c}.nii);
            for n=1:N
                fname       = dat.modality{m}.channel{c}.nii(n).dat.fname;
                [~,nam,ext] = fileparts(fname);
                nfname      = fullfile(dir_preproc,[nam ext]);                
                copyfile(fname,nfname);
                
                dat.modality{m}.channel{c}.nii(n) = nifti(nfname);
                
                fname       = dat.modality{m}.channel{c}.json(n).pth;
                [~,nam,ext] = fileparts(fname);
                nfname1     = fullfile(dir_preproc,[nam ext]);                
                copyfile(fname,nfname1);   
                
                dat.modality{m}.channel{c}.json(n).pth = nfname1;
                
                spm_json_manager('modify_json_field',dat.modality{m}.channel{c}.json(n).pth,'pth',nfname);
            end
        end
    else
        N = numel(dat.modality{m}.nii);
        for n=1:N
            fname       = dat.modality{m}.nii(n).dat.fname;
            [~,nam,ext] = fileparts(fname);
            nfname      = fullfile(dir_preproc,[nam ext]);            
            copyfile(fname,nfname);
            
            dat.modality{m}.nii(n) = nifti(nfname);
            
            fname       = dat.modality{m}.json(n).pth;
            [~,nam,ext] = fileparts(fname);
            nfname1     = fullfile(dir_preproc,[nam ext]);            
            copyfile(fname,nfname1);
            
            dat.modality{m}.json(n).pth = nfname1;
            
            spm_json_manager('modify_json_field',dat.modality{m}.json(n).pth,'pth',nfname);
        end
    end
end

if isfield(dat,'label')
    % Also copy labels, if exists
    %----------------------------------------------------------------------
    R = numel(dat.label);
    for r=1:R
        fname       = dat.label{r}.nii.dat.fname;
        [~,nam,ext] = fileparts(fname);
        nfname      = fullfile(dir_preproc,[nam ext]);            
        copyfile(fname,nfname);

        dat.label{r}.nii = nifti(nfname);

        fname       = dat.label{r}.json.pth;
        [~,nam,ext] = fileparts(fname);
        nfname1     = fullfile(dir_preproc,[nam ext]);            
        copyfile(fname,nfname1);

        dat.label{r}.json.pth = nfname1;

        spm_json_manager('modify_json_field',dat.label{r}.json.pth,'pth',nfname);
    end
end

% Rigidly realign to MNI space
%--------------------------------------------------------------------------
if preproc.do_realign2mni   
    origin_reset = false;
    for m=1:M
        if ~isfield(dat.modality{m},'channel')
            % Reset the origin (not compatible with multi-channel data)
            N = numel(dat.modality{m}.nii);
            for n=1:N                
                fname = dat.modality{m}.nii(n).dat.fname;
                mat0  = dat.modality{m}.nii(n).mat;
                vx    = spm_misc('vxsize',mat0);
                
                spm_impreproc('nm_reorient',fname,vx,1,'ro_');  
                
                [dat.modality{m}.nii(n),nfname] = update_nii(fname,'ro_');
                
                spm_json_manager('modify_json_field',dat.modality{m}.json(n).pth,'pth',nfname);
                
                fname = dat.modality{m}.nii(n).dat.fname;
                spm_impreproc('reset_origin',fname);
                
                dat.modality{m}.nii(n) = nifti(fname);
            end            
            
            origin_reset = true;
        end
        
        if isfield(dat.modality{m},'channel')
            fname                             = dat.modality{m}.channel{1}.nii(1).dat.fname;  
            mat1                              = spm_impreproc('rigid_align',fname);                     
            dat.modality{m}.channel{1}.nii(1) = nifti(fname);   
                
            C = numel(dat.modality{m}.channel);
            for c=1:C
                N = numel(dat.modality{m}.channel{c}.nii);
                for n=1:N
                    if c==1 && n==1, continue; end
                    
                    fname                             = dat.modality{m}.channel{c}.nii(n).dat.fname; 
                    mat0                              = dat.modality{m}.channel{c}.nii(n).mat; 
                    spm_get_space(fname,mat1\mat0);  
                    dat.modality{m}.channel{c}.nii(n) = nifti(fname);  
                end
            end            
        else
            N = numel(dat.modality{m}.nii);
            for n=1:N  
                fname                  = dat.modality{m}.nii(n).dat.fname;  
                mat1                   = spm_impreproc('rigid_align',fname);                     
                dat.modality{m}.nii(n) = nifti(fname);   
            end
        end
    end
    
    if isfield(dat,'label')
        R = numel(dat.label);
        for r=1:R
            fname = dat.label{r}.nii.dat.fname;

            if origin_reset
                mat0 = dat.label{r}.nii.mat;
                vx   = spm_misc('vxsize',mat0);

                spm_impreproc('nm_reorient',fname,vx,1,'ro_'); 


                [dat.label{r}.nii,nfname] = update_nii(fname,'ro_');

                spm_json_manager('modify_json_field',dat.label{r}.json.pth,'pth',nfname);

                fname = dat.label{r}.nii.dat.fname;
                spm_impreproc('reset_origin',fname);

                dat.label{r}.nii = nifti(fname);
            end

            mat0             = dat.label{r}.nii.mat;
            spm_get_space(fname,mat1\mat0); 
            dat.label{r}.nii = nifti(fname); 
        end
    end
end               

% Remove image data outside of the head (air..)
%--------------------------------------------------------------------------
if preproc.do_crop    
    for m=1:M
        if isfield(dat.modality{m},'channel')
            C = numel(dat.modality{m}.channel);
            for c=1:C
                N = numel(dat.modality{m}.channel{c}.nii);
                for n=1:N
                    fname = dat.modality{m}.channel{c}.nii(n).dat.fname;

                    [~,bb] = spm_impreproc('atlas_crop',fname,'cr_',preproc.do_rem_neck); 

                    [dat.modality{m}.channel{c}.nii(n),nfname] = update_nii(fname,'cr_');

                    spm_json_manager('modify_json_field',dat.modality{m}.channel{c}.json(n).pth,'pth',nfname);
                end
            end
        else
            N = numel(dat.modality{m}.nii);
            for n=1:N
                fname = dat.modality{m}.nii(n).dat.fname;

                [~,bb] = spm_impreproc('atlas_crop',fname,'cr_',preproc.do_rem_neck); 

                [dat.modality{m}.nii(n),nfname] = update_nii(fname,'cr_');

                spm_json_manager('modify_json_field',dat.modality{m}.json(n).pth,'pth',nfname);
            end
        end  
    end
    
    if isfield(dat,'label')
        R = numel(dat.label);
        for r=1:R
            fname = dat.label{r}.nii.dat.fname;
            
            spm_impreproc('subvol',spm_vol(fname),bb,'cr_'); 
                        
            [dat.label{r}.nii,nfname] = update_nii(fname,'cr_');
            
            spm_json_manager('modify_json_field',dat.label{r}.json.pth,'pth',nfname);                                            
        end
    end
end  

% Co-register images
%--------------------------------------------------------------------------
if preproc.do_coreg
    % First perform within-channel co-registration    
    for m=1:M
        if isfield(dat.modality{m},'channel')            
            V   = spm_vol;
            cnt = 0;            
            C   = numel(dat.modality{m}.channel);
            for c=1:C
                N = numel(dat.modality{m}.channel{c}.nii);
                for n=1:N
                    fname  = dat.modality{m}.channel{c}.nii(n).dat.fname;
                    cnt    = cnt + 1;
                    V(cnt) = spm_vol(fname);                    
                end
            end
            
            if cnt>1
                V = spm_impreproc('coreg',V);
                
                cnt = 0;
                for c=1:C
                    N = numel(dat.modality{m}.channel{c}.nii);
                    for n=1:N
                        cnt                               = cnt + 1;
                        dat.modality{m}.channel{c}.nii(n) = nifti(V(cnt).fname);                                                
                    end
                end
            end
        end
    end   

    % TODO: Then perform between-channel co-registration 
    %...
    
    clear V
end

% NN down-sampling in-plane
%--------------------------------------------------------------------------
if preproc.do_ds_inplane
    for m=1:M
        if isfield(dat.modality{m},'channel')
            C = numel(dat.modality{m}.channel);
            for c=1:C
                N = numel(dat.modality{m}.channel{c}.nii);
                for n=1:N
                    fname = dat.modality{m}.channel{c}.nii(n).dat.fname;

                    spm_impreproc('downsample_inplane',fname,1);

                    [dat.modality{m}.channel{c}.nii(n),nfname] = update_nii(fname,'ds_');

                    spm_json_manager('modify_json_field',dat.modality{m}.channel{c}.json(n).pth,'pth',nfname);
                end
            end
        else
            N = numel(dat.modality{m}.nii);
            for n=1:N
                fname = dat.modality{m}.nii(n).dat.fname;

                spm_impreproc('downsample_inplane',fname,1);

                [dat.modality{m}.nii(n),nfname] = update_nii(fname,'ds_');

                spm_json_manager('modify_json_field',dat.modality{m}.json(n).pth,'pth',nfname);
            end
        end  
    end
    
    if isfield(dat,'label')
        R = numel(dat.label);
        for r=1:R
            fname = dat.label{r}.nii.dat.fname;
            
            spm_impreproc('downsample_inplane',fname,1);
                        
            [dat.label{r}.nii,nfname] = update_nii(fname,'ds_');
            
            spm_json_manager('modify_json_field',dat.label{r}.json.pth,'pth',nfname);                                            
        end
    end    
end

% Create equally sized images by super-resolution
%--------------------------------------------------------------------------
if preproc.do_superres    
    for m=1:M    
        modality = dat.modality{m}.name;
        if ~strcmpi(modality,'MRI'), continue; end
                 
        if isfield(dat.modality{m},'channel')
            C   = numel(dat.modality{m}.channel);
            img = cell(1,C); % Input object to super-resolution algorithm 
            for c=1:C
                img{c} = nifti;
                N = numel(dat.modality{m}.channel{c}.nii);
                for n=1:N
                    fname     = dat.modality{m}.channel{c}.nii(n).dat.fname;
                    img{c}(n) = nifti(fname);
                end
            end
        else
            img{1} = nifti;
            N      = numel(dat.modality{m}.nii);            
            for n=1:N
                fname        = dat.modality{m}.nii(n).dat.fname;
                img{1}(n) = nifti(fname);
            end
        end
        
        % Do super-resolution
        nii = spm_superres(img,modality,superres);
        clear img   
        
        if isfield(dat.modality{m},'channel')
            C = numel(dat.modality{m}.channel);
            for c=1:C
                nfname                            = nii(c).dat.fname;
                dat.modality{m}.channel{c}.nii(1) = nii(c);
                spm_json_manager('modify_json_field',dat.modality{m}.channel{c}.json(1).pth,'pth',nfname);
                
                N = numel(dat.modality{m}.channel{c}.nii);
                for n=2:N                      
                    delete(dat.modality{m}.channel{c}.json(n).pth);
                    dat.modality{m}.channel{c}.nii(n)  = [];
                    dat.modality{m}.channel{c}.json(n) = [];
                end
            end
        else
            nfname                 = nii(c).dat.fname;
            dat.modality{m}.nii(1) = nii(c);
            spm_json_manager('modify_json_field',dat.modality{m}.json(1).pth,'pth',nfname);

            N = numel(dat.modality{m}.nii);
            for n=2:N                      
                delete(dat.modality{m}.json(n).pth);
                dat.modality{m}.nii(n)  = [];
                dat.modality{m}.json(n) = [];
            end
        end        
    end            
end                  

% Change voxel size of image(s)
%--------------------------------------------------------------------------
if ~isempty(preproc.vx) && ~preproc.do_superres
    for m=1:M
        if isfield(dat.modality{m},'channel')
            C = numel(dat.modality{m}.channel);
            for c=1:C
                N = numel(dat.modality{m}.channel{c}.nii);
                for n=1:N
                    fname = dat.modality{m}.channel{c}.nii(n).dat.fname;

                    spm_impreproc('nm_reorient',fname,preproc.vx,1,'vx_');  

                    [dat.modality{m}.channel{c}.nii(n),nfname] = update_nii(fname,'vx_');

                    spm_json_manager('modify_json_field',dat.modality{m}.channel{c}.json(n).pth,'pth',nfname);
                end
            end
        else
            N = numel(dat.modality{m}.nii);
            for n=1:N
                fname = dat.modality{m}.nii(n).dat.fname;

                spm_impreproc('nm_reorient',fname,preproc.vx,1,'vx_');  

                [dat.modality{m}.nii(n),nfname] = update_nii(fname,'vx_');

                spm_json_manager('modify_json_field',dat.modality{m}.json(n).pth,'pth',nfname);
            end
        end       
    end
    
    if isfield(dat,'label')
        R = numel(dat.label);
        for r=1:R
            fname = dat.label{r}.nii.dat.fname;
            
            spm_impreproc('nm_reorient',fname,preproc.vx,0,'vx_');  
                        
            [dat.label{r}.nii,nfname] = update_nii(fname,'vx_');
            
            spm_json_manager('modify_json_field',dat.label{r}.json.pth,'pth',nfname);                                            
        end
    end    
end

% Reslice to size of image with largest FOV
%--------------------------------------------------------------------------
if preproc.do_reslice && ~preproc.do_superres       
    V   = spm_vol;
    cnt = 0;
    for m=1:M
        if isfield(dat.modality{m},'channel')
            C = numel(dat.modality{m}.channel);
            for c=1:C
                N = numel(dat.modality{m}.channel{c}.nii);
                for n=1:N
                    fname  = dat.modality{m}.channel{c}.nii(n).dat.fname;                    
                    cnt    = cnt + 1;
                    V(cnt) = spm_vol(fname);                    
                end
            end
        else
            N = numel(dat.modality{m}.nii);
            for n=1:N
                fname  = dat.modality{m}.nii(n).dat.fname;
                cnt    = cnt + 1;
                V(cnt) = spm_vol(fname);                
            end
        end  
    end
    
    if cnt>1
        V = spm_impreproc('reslice',V); 
        
        cnt = 0;
        for m=1:M
            if isfield(dat.modality{m},'channel')
                C = numel(dat.modality{m}.channel);
                for c=1:C
                    N = numel(dat.modality{m}.channel{c}.nii);
                    for n=1:N
                        cnt                               = cnt + 1;                                           
                        dat.modality{m}.channel{c}.nii(n) = nifti(V(cnt).fname);  
                        
                        spm_json_manager('modify_json_field',dat.modality{m}.channel{c}.json(n).pth,'pth',V(cnt).fname);
                    end
                end
            else
                N = numel(dat.modality{m}.nii);
                for n=1:N
                    cnt                    = cnt + 1;                                           
                    dat.modality{m}.nii(n) = nifti(V(cnt).fname);
                    
                    spm_json_manager('modify_json_field',dat.modality{m}.json(n).pth,'pth',V(cnt).fname);
                end
            end  
        end
    end
    
    clear V
end

% Denoise images
%--------------------------------------------------------------------------
if preproc.do_denoise
    for m=1:M    
        nii      = nifti;
        modality = dat.modality{m}.name;
        cnt      = 0;
        if isfield(dat.modality{m},'channel')
            C = numel(dat.modality{m}.channel);
            for c=1:C
                N = numel(dat.modality{m}.channel{c}.nii);
                for n=1:N
                    cnt      = cnt + 1;
                    fname    = dat.modality{m}.channel{c}.nii(n).dat.fname;
                    nii(cnt) = nifti(fname);
                end
            end
        else
            N = numel(dat.modality{m}.nii);
            for n=1:N
                cnt      = cnt + 1;
                fname    = dat.modality{m}.nii(n).dat.fname;
                nii(cnt) = nifti(fname);
            end
        end
        
        % Do denoising
        nii = spm_denoise(nii,modality,denoise);  
        
        cnt = 0;
        if isfield(dat.modality{m},'channel')
            C = numel(dat.modality{m}.channel);
            for c=1:C
                N = numel(dat.modality{m}.channel{c}.nii);
                for n=1:N
                    cnt                               = cnt + 1;
                    dat.modality{m}.channel{c}.nii(n) = nii(cnt);   
                    
                    nfname = nii(cnt).dat.fname;
                    spm_json_manager('modify_json_field',dat.modality{m}.channel{c}.json(n).pth,'pth',nfname);
                end
            end
        else
            N = numel(dat.modality{m}.nii);
            for n=1:N
                cnt                    = cnt + 1;
                dat.modality{m}.nii(n) = nii(cnt);     
                
                nfname = nii(cnt).dat.fname;
                spm_json_manager('modify_json_field',dat.modality{m}.json(n).pth,'pth',nfname);
            end
        end        
    end    
    
    clear nii
end

% Simple normalisation of image intensities (make mean(img)=100)
%--------------------------------------------------------------------------
if preproc.do_normalise_intensities    
    for m=1:M
        if isfield(dat.modality{m},'channel')
            C = numel(dat.modality{m}.channel);
            for c=1:C
                N = numel(dat.modality{m}.channel{c}.nii);
                for n=1:N                       
                    img = dat.modality{m}.channel{c}.nii(n).dat(:,:,:);                    
                    img = normalise_intensities(img,dat.modality{m}.name);                    

                    dat.modality{m}.channel{c}.nii(n).dat(:,:,:) = img;
                end
            end
        else
            N = numel(dat.modality{m}.nii);
            for n=1:N
                img = dat.modality{m}.nii(n).dat(:,:,:);                    
                img = normalise_intensities(img,dat.modality{m}.name);                    

                dat.modality{m}.nii(n).dat(:,:,:) = img;
            end
        end  
    end
    clear img
end                  

% Segment using spm_preproc8
%--------------------------------------------------------------------------
if preproc.do_segment || preproc.do_skull_strip || preproc.do_bf_correct              

    % For segmenting we assume... 
    m = 1; % ...one modality...
    n = 1; % ...and one image per channel        
    
    % Set-up output
    write_tc = true(6,4); % c, rc, wc, mwc
    write_bf = [false true]; 
    write_df = false(3,1);                       
    
    % Build spm_vol input-object to segmentation function
    if isfield(dat.modality{m},'channel')
        V = spm_vol;
        C = numel(dat.modality{m}.channel);
        for c=1:C            
            fname = dat.modality{m}.channel{c}.nii(n).dat.fname;
            V(c)  = spm_vol(fname);
        end
    else        
        fname = dat.modality{m}.nii(n).dat.fname;
        V     = spm_vol(fname);
    end  
    
    % Run segmentation
    segment_subject(V,write_tc,write_bf,write_df,dir_preproc,dat.modality{m}.name);
     
    % Get output files
    dat = add_segmentations2dat(dat);
end         

% Overwrite image data with bias-corrected version from spm_preproc8
%--------------------------------------------------------------------------
if preproc.do_bf_correct
    if isfield(dat.modality{m},'channel')            
        C = numel(dat.modality{m}.channel);
        for c=1:C            
            fname = dat.modality{m}.channel{c}.nii(n).dat.fname;

            [dat.modality{m}.channel{c}.nii(n),nfname] = update_nii(fname,'m');

            spm_json_manager('modify_json_field',dat.modality{m}.channel{c}.json(n).pth,'pth',nfname);
        end
    else        
        fname = dat.modality{m}.nii(n).dat.fname;

        [dat.modality{m}.nii(n),nfname] = update_nii(fname,'m');

        spm_json_manager('modify_json_field',dat.modality{m}.json(n).pth,'pth',nfname);
    end              
end

% Skull-strip image data using mask built from GM, WM and CSF classes from
% spm_preproc8 segmentation
%--------------------------------------------------------------------------
if preproc.do_skull_strip    
    dat = skull_strip(dat);                 
end  

% Clean-up
%--------------------------------------------------------------------------
if ~preproc.do_segment && (preproc.do_skull_strip || preproc.do_bf_correct)
    for t=1:numel(dat.segmentation)
        for k=1:numel(dat.segmentation{t}.class)
            delete(dat.segmentation{t}.class{k}.nii.dat.fname);
            delete(dat.segmentation{t}.class{k}.json.pth);
        end
    end
    
    dat = rmfield(dat,'segmentation');
    dat = rmfield(dat,'segmentation_map');
end

% Create 2D versions
%--------------------------------------------------------------------------
if write_2d
    % Make copies
    for m=1:M
        if isfield(dat.modality{m},'channel')
            C = numel(dat.modality{m}.channel);
            for c=1:C
                N = numel(dat.modality{m}.channel{c}.nii);
                for n=1:N
                    fname       = dat.modality{m}.channel{c}.nii(n).dat.fname;
                    [~,nam,ext] = fileparts(fname);
                    nfname      = fullfile(dir_2d,[nam ext]);                
                    copyfile(fname,nfname);

                    dat.modality{m}.channel{c}.nii(n) = nifti(nfname);

                    fname       = dat.modality{m}.channel{c}.json(n).pth;
                    [~,nam,ext] = fileparts(fname);
                    nfname      = fullfile(dir_2d,[nam ext]);                
                    copyfile(fname,nfname);   

                    dat.modality{m}.channel{c}.json(n).pth = nfname;
                    
                    spm_json_manager('modify_json_field',dat.modality{m}.channel{c}.json(n).pth,'pth',nfname);
                end
            end
        else
            N = numel(dat.modality{m}.nii);
            for n=1:N
                fname       = dat.modality{m}.nii(n).dat.fname;
                [~,nam,ext] = fileparts(fname);
                nfname      = fullfile(dir_2d,[nam ext]);            
                copyfile(fname,nfname);

                dat.modality{m}.nii(n) = nifti(nfname);

                fname       = dat.modality{m}.json(n).pth;
                [~,nam,ext] = fileparts(fname);
                nfname      = fullfile(dir_2d,[nam ext]);            
                copyfile(fname,nfname);

                dat.modality{m}.json(n).pth = nfname;
                
                spm_json_manager('modify_json_field',dat.modality{m}.json(n).pth,'pth',nfname);
            end
        end
    end  
    
    if isfield(dat,'label')
        % Also copy labels, if exists
        %------------------------------------------------------------------
        R = numel(dat.label);
        for r=1:R
            fname       = dat.label{r}.nii.dat.fname;
            [~,nam,ext] = fileparts(fname);
            nfname      = fullfile(dir_2d,[nam ext]);            
            copyfile(fname,nfname);

            dat.label{r}.nii = nifti(nfname);

            fname       = dat.label{r}.json.pth;
            [~,nam,ext] = fileparts(fname);
            nfname1     = fullfile(dir_2d,[nam ext]);            
            copyfile(fname,nfname1);

            dat.label{r}.json.pth = nfname1;

            spm_json_manager('modify_json_field',dat.label{r}.json.pth,'pth',nfname);
        end
    end

    if isfield(dat,'segmentation')
        % Also copy segmentations, if exists
        %------------------------------------------------------------------
        T = numel(dat.segmentation);
        for t=1:T
            K = numel(dat.segmentation{t}.class);
            for k=1:K
                fname       = dat.segmentation{t}.class{k}.nii.dat.fname;
                [~,nam,ext] = fileparts(fname);
                nfname      = fullfile(dir_2d,[nam ext]);            
                copyfile(fname,nfname);

                dat.segmentation{t}.class{k}.nii = nifti(nfname);
               
                fname       = dat.segmentation{t}.class{k}.json.pth;
                [~,nam,ext] = fileparts(fname);
                nfname1     = fullfile(dir_2d,[nam ext]);            
                copyfile(fname,nfname1);
                
                dat.segmentation{t}.class{k}.json.pth = nfname1;
                
                spm_json_manager('modify_json_field',dat.segmentation{t}.class{k}.json.pth,'pth',nfname);
           end
        end
    end
    
    % Create 2d-slices
    %----------------------------------------------------------------------
    for m=1:M
        if isfield(dat.modality{m},'channel')
            C = numel(dat.modality{m}.channel);
            for c=1:C
                N = numel(dat.modality{m}.channel{c}.nii);
                for n=1:N
                    fname = dat.modality{m}.channel{c}.nii(n).dat.fname;

                    nfname = create_2d_slice(fname,axis_2d);

                    dat.modality{m}.channel{c}.nii(n) = nifti(nfname);

                    spm_json_manager('modify_json_field',dat.modality{m}.channel{c}.json(n).pth,'pth',nfname);
                end
            end
        else
            N = numel(dat.modality{m}.nii);
            for n=1:N
                fname = dat.modality{m}.nii(n).dat.fname;

                nfname = create_2d_slice(fname,axis_2d);

                dat.modality{m}.nii(n) = nifti(nfname);

                spm_json_manager('modify_json_field',dat.modality{m}.json(n).pth,'pth',nfname);
            end
        end   
    end
    
    if isfield(dat,'label')
        R = numel(dat.label);
        for r=1:R
            fname  = dat.label{r}.nii.dat.fname;
            
            nfname = create_2d_slice(fname,axis_2d);

            dat.label{r}.nii = nifti(nfname);

            spm_json_manager('modify_json_field',dat.label{r}.json.pth,'pth',nfname);
        end
    end    
    
    if isfield(dat,'segmentation')
        T = numel(dat.segmentation);
        for t=1:T
            K = numel(dat.segmentation{t}.class);
            for k=1:K
                fname = dat.segmentation{t}.class{k}.nii.dat.fname;
                
                nfname = create_2d_slice(fname,axis_2d);
                
                dat.segmentation{t}.class{k}.nii = nifti(nfname);

                spm_json_manager('modify_json_field',dat.segmentation{t}.class{k}.json.pth,'pth',nfname);
           end
        end
    end    
end      
%==========================================================================

%==========================================================================
function [nii,nfname] = update_nii(fname,prefix)  
[pth,nam,ext] = fileparts(fname);
delete(fname);
nfname        = fullfile(pth,[prefix nam ext]);        
nii           = nifti(nfname);
%==========================================================================

%==========================================================================
function img = normalise_intensities(img,modality,val)               
if nargin<3, val = 100; end

msk  = spm_misc('msk_modality',img,modality);     
sint = sum(reshape(img(msk),[],1));
nm   = nnz(msk);        
scl  = (val/(sint/nm));
img  = scl*img;
%==========================================================================   