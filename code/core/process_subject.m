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
segment     = opt.segment;
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

% Down-sample in-plane
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

% Simple normalisation of image intensities
%--------------------------------------------------------------------------
if preproc.normalise_intensities    
    for m=1:M
        if isfield(dat.modality{m},'channel')
            C = numel(dat.modality{m}.channel);
            for c=1:C
                N = numel(dat.modality{m}.channel{c}.nii);
                for n=1:N                       
                    img = dat.modality{m}.channel{c}.nii(n).dat(:,:,:);                    
                    img = normalise_intensities(img,dat.modality{m}.name);                    

                    fname = dat.modality{m}.channel{c}.nii(n).dat.fname;
                    mat   = dat.modality{m}.channel{c}.nii(n).dat.mat;
                    dtype = dat.modality{m}.channel{c}.nii(n).dat.dtype;
                    
                    [pth,nam,ext] = fileparts(fname);
                    nfname        = fullfile(pth,['ni_' nam ext]);
                    spm_misc('create_nii',nfname,img,mat,dtype,'norm-int');
                    
                    dat.modality{m}.channel{c}.nii(n) = nifti(nfname);

                    spm_json_manager('modify_json_field',dat.modality{m}.channel{c}.json(n).pth,'pth',nfname);
                end
            end
        else
            N = numel(dat.modality{m}.nii);
            for n=1:N
                img = dat.modality{m}.nii(n).dat(:,:,:);                    
                img = normalise_intensities(img,dat.modality{m}.name);                    

                fname = dat.modality{m}.nii(n).dat.fname;
                mat   = dat.modality{m}.nii(n).dat.mat;
                dtype = dat.modality{m}.nii(n).dat.dtype;

                [pth,nam,ext] = fileparts(fname);
                nfname        = fullfile(pth,['ni_' nam ext]);
                spm_misc('create_nii',nfname,img,mat,dtype,'norm-int');

                dat.modality{m}.nii(n) = nifti(nfname);

                spm_json_manager('modify_json_field',dat.modality{m}.json(n).pth,'pth',nfname);
            end
        end  
    end
end                  

% Segment and/or bias-field correct and/or skull-strip
%--------------------------------------------------------------------------
if segment.write_tc || preproc.do_skull_strip || preproc.do_bf_correct              

    % For segmenting we assume one modality and one image per channel
    m = 1;
    n = 1;
    
    tissue = {'GM','WM','CSF','BONE','ST','BG'};
    
    % Set-up output
    write_tc = false(6,4);
    write_bf = false(1,2); 
    write_df = false(3,1);        
    if preproc.do_skull_strip
        write_tc(1:3,1) = true;
    end    
    if segment.write_tc
        write_tc(1:6,[1 2 3 4]) = true;
    end    
    if preproc.do_bf_correct
        write_bf(1,2) = true;
    end
       
    K_a = nnz(sum(write_tc,2));
    
    % Build spm_vol object    
    if isfield(dat.modality{m},'channel')
        V = spm_vol;
        C = numel(dat.modality{m}.channel);
        for c=1:C            
            fname = dat.modality{m}.channel{c}.nii(n).dat.fname;
            V(c)  = spm_vol(fname);
            
            pth_json = dat.modality{m}.channel{c}.json(n).pth;
        end
    else        
        fname = dat.modality{m}.nii(n).dat.fname;
        V     = spm_vol(fname);
        
        pth_json = dat.modality{m}.json(n).pth;
    end  
    
    % Run segmentation
    segment_subject(V,write_tc,write_bf,write_df,dir_preproc,dat.modality{m}.name);
     
    % Get output files
    a0 = spm_jsonread(pth_json);
    
    for c=1:numel(V)        
        fname0     = V(c).fname;
        [~,nam0,~] = fileparts(fname0);
        
        files = spm_select('FPList',dir_preproc,nam0);
        
        for i=1:size(files,1)
            fname1        = files(i,:);
            [pth1,nam1,~] = fileparts(fname1);
            
            
            for k=1:K_a
                type = 'wc';
                if strcmp(nam1,[type num2str(k) nam0])                                        
                    a            = struct;
                    a.name       = a0.name;
                    a.population = a0.population;
                    a.pth        = fname1;
                    a.type       = type;
                    a.tissue     = tissue{k};
                    
                    pth_json1 = fullfile(pth1,[nam1 '.json']);
                    
                    a = orderfields(a);
                    spm_jsonwrite(pth_json1,a);
                end      
                
                type = 'c';
                if strcmp(nam1,[type num2str(k) nam0])                                        
                    a            = struct;
                    a.name       = a0.name;
                    a.population = a0.population;
                    a.pth        = fname1;
                    a.type       = type;
                    a.tissue     = tissue{k};
                    
                    pth_json1 = fullfile(pth1,[nam1 '.json']);
                    
                    a = orderfields(a);
                    spm_jsonwrite(pth_json1,a);
                end   
                
                type = 'rc';
                if strcmp(nam1,[type num2str(k) nam0])                                        
                    a            = struct;
                    a.name       = a0.name;
                    a.population = a0.population;
                    a.pth        = fname1;
                    a.type       = type;
                    a.tissue     = tissue{k};
                    
                    pth_json1 = fullfile(pth1,[nam1 '.json']);
                    
                    a = orderfields(a);
                    spm_jsonwrite(pth_json1,a);
                end    
                
                type = 'mwc';
                if strcmp(nam1,[type num2str(k) nam0])                                        
                    a            = struct;
                    a.name       = a0.name;
                    a.population = a0.population;
                    a.pth        = fname1;
                    a.type       = type;
                    a.tissue     = tissue{k};
                    
                    pth_json1 = fullfile(pth1,[nam1 '.json']);
                    
                    a = orderfields(a);
                    spm_jsonwrite(pth_json1,a);
                end                  
            end
        end
    end
    
%     if obj.preproc.do_bf_correct
%         % Replace data with bias-corrected versions
%         files = spm_select('FPList',dir_seg,'^m.*\.nii$'); % Get bias-field corrected images
%         for n=1:C          
%             fname   = obj.scans{n}{1}.fname;                        
%             [~,nam] = fileparts(fname);
%             for n1=1:C
%                 [~,nam_bf] = fileparts(files(n1,:));
%                 if strcmp(nam,nam_bf(2:end))
%                     fname_bf = files(n1,:);
%                     V_bf     = spm_vol(fname_bf);
%                     
%                     [pth,nam,ext] = fileparts(fname);
%                     nfname        = fullfile(pth,['bf_' nam ext]);
% 
%                     Nii         = nifti;
%                     Nii.dat     = file_array(nfname,V_bf.dim,V_bf.dt,0,1,0);
%                     Nii.mat     = V_bf.mat;
%                     Nii.mat0    = V_bf.mat;
%                     Nii.descrip = 'bf-corrected';
%                     create(Nii);
% 
%                     Nii1           = nifti(fname_bf);
%                     img            = single(Nii1.dat(:,:,:));                    
%                     Nii.dat(:,:,:) = img;
%                     
%                     obj.scans{n}{1} = spm_vol(nfname);
%                     delete(fname)
%                 end
%             end  
%         end
%     end
%     
%     if obj.preproc.do_skull_strip
%         % Overwrite image data with skull-stripped versions
%         files = spm_select('FPList',dir_seg,'^c[1,2,3].*\.nii$');
%         V0    = spm_vol(files);
%         K     = numel(V0);
%         msk   = zeros(V0(1).dim,'single');
%         for k=1:K
%             Nii  = nifti(V0(k).fname);
%             resp = single(Nii.dat(:,:,:)); 
%             msk  = msk + resp;
%         end
%         clear resp
%                                                
%         for z=1:V0(1).dim(3) % loop over slices
%             msk(:,:,z) = imgaussfilt(msk(:,:,z),1);    % Smooth
%             msk(:,:,z) = msk(:,:,z)>0.5;               % Threshold
%             msk(:,:,z) = imfill(msk(:,:,z),4,'holes'); % Fill holes
%         end
% 
%         % Mask out voxels based on SPM TPM size
%         pth_tpm = fullfile(spm('dir'),'tpm','TPM.nii,');
%         V1      = spm_vol(pth_tpm);
% 
%         M0  = V0(1).mat;      
%         dm0 = V0(1).dim; 
%         M1  = V1(1).mat;  
%         dm1 = V1(1).dim;
% 
%         T = M1\M0; % Get the mapping from M0 to M1
% 
%         % Use ndgrid to give an array of voxel indices
%         [x0,y0,z0] = ndgrid(single(1:dm0(1)),...
%                             single(1:dm0(2)),...
%                             single(1:dm0(3)));
% 
%         % Transform these indices to the indices that they point to in the reference image
%         D = cat(4,T(1,1)*x0 + T(1,2)*y0 + T(1,3)*z0 + T(1,4), ...
%                   T(2,1)*x0 + T(2,2)*y0 + T(2,3)*z0 + T(2,4), ...
%                   T(3,1)*x0 + T(3,2)*y0 + T(3,3)*z0 + T(3,4));
%         clear x0 y0 z0
%         
%         % Mask according to whether these are < 1 or > than the dimensions of the reference image.        
%         msk1 = cell(1,3);
%         ix   = [1 1 20];
%         for i=1:3
%             msk1{i} = D(:,:,:,i) >= ix(i) & D(:,:,:,i) <= dm1(i);
%         end
%         clear D
%         
%         % Generate masked image
%         for i=1:3
%             msk = msk1{i}.*msk;
%         end
%         clear msk1
%         
%         if 0
%             % For testing skull-stripping
%             split = 4;
%             dm0   = V0(1).dim;
%             nfigs = floor(dm0(3)/split);
%             
%             F1 = floor(sqrt(nfigs));
%             F2 = ceil(nfigs/F1);      
%             figure(666); 
%             for f=1:nfigs
%                 subplot(F1,F2,f);            
%                 imagesc(msk(:,:,split*f)'); colormap(gray); axis off xy image;
%             end
%             
%             figure(667);            
%             msk1 = permute(msk,[2 3 1]);
%             slice = msk1(:,:,floor(size(msk1,3)/2) + 1);
%             subplot(121); imagesc(slice'); colormap(gray); axis off xy;
%             
%             msk1 = permute(obj.scans{1}{1}.private.dat(:,:,:),[2 3 1]);
%             slice = msk1(:,:,floor(size(msk1,3)/2) + 1);
%             subplot(122); imagesc(slice',[0 100]); colormap(gray); axis off xy;
%         end
%         
%         for n=1:C
%             fname         = obj.scans{n}{1}.fname;  
%             [pth,nam,ext] = fileparts(fname);
%             nfname        = fullfile(pth,['ss_' nam ext]);
% 
%             Nii         = nifti;
%             Nii.dat     = file_array(nfname,obj.scans{n}{1}.dim,[16 0],0,1,0);
%             Nii.mat     = obj.scans{n}{1}.mat;
%             Nii.mat0    = obj.scans{n}{1}.mat;
%             Nii.descrip = 'skull-stripped';
%             create(Nii);
% 
%             Nii1           = nifti(fname);
%             img            = single(Nii1.dat(:,:,:));
%             img(~msk)      = NaN;  
%             Nii.dat(:,:,:) = img;
%             
%             obj.scans{n}{1} = spm_vol(nfname);
%             delete(fname);
%         end
%     end  
%     
%     if obj.preproc.make_ml_labels && isempty(obj.labels)
%         % Write maximum-likelihoods labels (only if labels are not available)                
%         files = spm_select('FPList',dir_seg,'^c.*\.nii$');
%         
%         if ~isempty(files)
%             V0    = spm_vol(files);
%             K     = numel(V0);
%             img   = zeros([V0(1).dim K],'single');
%             for k=1:K
%                 Nii          = nifti(V0(k).fname);
%                 img(:,:,:,k) = single(Nii.dat(:,:,:));
%             end
% 
%             if K<6
%                 % Less than the default SPM template number of classes requested
%                 % -> correct ML labels
%                 img1 = ones(V0(1).dim,'single');
%                 img1 = img1 - sum(img,4);
%                 img  = cat(4,img,img1);
%                 clear img1
%             end
% 
%             [~,ml_labels] = max(img,[],4);        
%             clear img
% 
%             fname       = obj.scans{1}{1}.fname;  
%             [~,nam,ext] = fileparts(fname);
%             nfname      = fullfile(dir_labels,['ml_' nam ext]);
% 
%             Nii      = nifti;
%             Nii.dat  = file_array(nfname,size(ml_labels),'uint8',0,1/K,0);
%             Nii.mat  = V0(1).mat;
%             Nii.mat0 = V0(1).mat;
%             Nii.descrip = 'ML-labels';
%             create(Nii);
% 
%             Nii.dat(:,:,:) = ml_labels;
%             clear ml_labels
% 
%             obj.labels = spm_vol(nfname);
%         end
%     end
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
        %----------------------------------------------------------------------
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

    % Create 2d-slices
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