function dat = process_subject(dat,opt) 

%--------------------------------------------------------------------------
% 1.  Make copies
% 2.  Reset origin
% 3.  Rigidly realign to MNI space
% 4.  Remove image data outside of the head
% 5.  Coregister images
% 6.  NN down-sampling in-plane
% 7.  NN down-sampling through-plane
% 8.  Create equally sized images by super-resolution
% 9.  Reslice to size of image with largest FOV
% 10. Change voxel size of image(s)
% 11. Resize to SPM template dimensions
% 12. Denoise images
% 13. Simple normalisation of image intensities
% 14. Segment using spm_preproc8
% 15. Overwrite image data with bias-corrected version from spm_preproc8
% 16. Skull-strip image data
% 17. Write 2D versions
%--------------------------------------------------------------------------

% Get NN down-sampling through-planeoptions
%--------------------------------------------------------------------------
dir_preproc = opt.dir_preproc;
dir_2d      = opt.dir_2d;
axis_2d     = opt.axis_2d;
write_2d    = opt.write_2d;
preproc     = opt.preproc;
M           = numel(dat.modality);
deg         = preproc.deg;

% Get indices of highest resolution image
%--------------------------------------------------------------------------
if numel(dat.modality) > 1
    error('numel(dat.modality) > 1');
end
m              = 1; % We assume only one modality for now..
[ixc,ixn,ixhr] = get_ix_hr(dat.modality{m});

% Make copies into dir_preproc (in order to not modify the original data)
%--------------------------------------------------------------------------
fprintf('\n')
fprintf('Copying images... ')

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
                
                if isfield(dat.modality{m}.channel{c},'json')
                    fname       = dat.modality{m}.channel{c}.json(n).pth;
                    [~,nam,ext] = fileparts(fname);
                    nfname1     = fullfile(dir_preproc,[nam ext]);                
                    copyfile(fname,nfname1);   

                    dat.modality{m}.channel{c}.json(n).pth = nfname1;

                    spm_json_manager('modify_json_field',dat.modality{m}.channel{c}.json(n).pth,'pth',nfname);
                end
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
            
            if isfield(dat.modality{m},'json')
                fname       = dat.modality{m}.json(n).pth;
                [~,nam,ext] = fileparts(fname);
                nfname1     = fullfile(dir_preproc,[nam ext]);            
                copyfile(fname,nfname1);

                dat.modality{m}.json(n).pth = nfname1;

                spm_json_manager('modify_json_field',dat.modality{m}.json(n).pth,'pth',nfname);
            end
        end
    end
end
fprintf('done!\n')

if isfield(dat,'label')
    fprintf('Copying labels... ')
    
    % Also copy labels, if exists
    %----------------------------------------------------------------------
    R = numel(dat.label);
    for r=1:R
        fname       = dat.label{r}.nii.dat.fname;
        [~,nam,ext] = fileparts(fname);
        nfname      = fullfile(dir_preproc,[nam ext]);            
        copyfile(fname,nfname);

        dat.label{r}.nii = nifti(nfname);

        if isfield(dat.label{r},'json')
            fname       = dat.label{r}.json.pth;
            [~,nam,ext] = fileparts(fname);
            nfname1     = fullfile(dir_preproc,[nam ext]);            
            copyfile(fname,nfname1);

            dat.label{r}.json.pth = nfname1;

            spm_json_manager('modify_json_field',dat.label{r}.json.pth,'pth',nfname);
        end
    end
    
    fprintf('done!\n')
end

% Reset origin
%--------------------------------------------------------------------------
if preproc.reset_origin
    fprintf('Resetting origin... ')
            
    for m=1:M        
        n = 1;
        c = 1;
        if isfield(dat.modality{m},'channel')                        
            fname = dat.modality{m}.channel{c}.nii(n).dat.fname;
            mat0  = dat.modality{m}.channel{c}.nii(n).mat;
        else
            fname = dat.modality{m}.nii(n).dat.fname;
            mat0  = dat.modality{m}.nii(n).mat;
        end            
        vx = spm_misc('vxsize',mat0);
        
        spm_impreproc('nm_reorient',fname,vx,'ro_',deg);  
        if isfield(dat.modality{m},'channel')
            [dat.modality{m}.channel{c}.nii(n),nfname] = update_nii(fname,'ro_');
        else
            [dat.modality{m}.nii(n),nfname]            = update_nii(fname,'ro_');
        end
        
        if isfield(dat.modality{m},'channel') && isfield(dat.modality{m}.channel{c},'json')
            spm_json_manager('modify_json_field',dat.modality{m}.channel{c}.json(n).pth,'pth',nfname);
        elseif isfield(dat.modality{m},'json')
            spm_json_manager('modify_json_field',dat.modality{m}.json(n).pth,'pth',nfname);
        end

        spm_impreproc('reset_origin',nfname);

        if isfield(dat.modality{m},'channel')
            dat.modality{m}.channel{c}.nii(n) = nifti(nfname);
        else
            dat.modality{m}.nii(n)            = nifti(nfname);
        end
    end
    
    if isfield(dat,'label')
        % Do labels too
        R = numel(dat.label);
        for r=1:R          
            fname = dat.label{r}.nii.dat.fname;
            mat0  = dat.label{r}.nii.mat;
            vx    = spm_misc('vxsize',mat0);
            
            spm_impreproc('nm_reorient',fname,vx,'ro_',deg); 

            [dat.label{r}.nii,nfname] = update_nii(fname,'ro_');

            if isfield(dat.label{r},'json')
                spm_json_manager('modify_json_field',dat.label{r}.json.pth,'pth',nfname);
            end

            spm_impreproc('reset_origin',nfname);
            
            dat.label{r}.nii = nifti(nfname);
        end
    end
        
    fprintf('done!\n')
end

% Collapse labels
%--------------------------------------------------------------------------
if isfield(dat,'label') && ~isempty(preproc.part_labels)    
    R = numel(dat.label);
    for r=1:R
        collapse_labels(dat.label{r}.nii,preproc.part_labels);                 
    end    
end

% Rigidly realign to MNI space
%--------------------------------------------------------------------------
if preproc.do_realign2mni   
    fprintf('Realigning to MNI... ')
    
    for m=1:M % Loop over modalities        
        
        if isfield(dat.modality{m},'channel')                        
            fname                                 = dat.modality{m}.channel{ixc}.nii(ixn).dat.fname;  
            mat1                                  = spm_impreproc('rigid_align',fname);                     
            dat.modality{m}.channel{ixc}.nii(ixn) = nifti(fname);   
                
            C = numel(dat.modality{m}.channel);
            for c=1:C
                
                N = numel(dat.modality{m}.channel{c}.nii);
                for n=1:N

                    
                    if c==ixc && n==ixn
                        continue; 
                    end
                    
                    fname                             = dat.modality{m}.channel{c}.nii(n).dat.fname; 
                    mat0                              = dat.modality{m}.channel{c}.nii(n).mat; 
                    spm_get_space(fname,mat1\mat0);  
                    dat.modality{m}.channel{c}.nii(n) = nifti(fname);                      
                end                
            end            
        else                        
            fname                    = dat.modality{m}.nii(ixn).dat.fname;  
            mat1                     = spm_impreproc('rigid_align',fname);                     
            dat.modality{m}.nii(ixn) = nifti(fname);   

            N = numel(dat.modality{m}.nii);
            for n=1:N

                if n==ixn
                    continue; 
                end

                fname                  = dat.modality{m}.nii(n).dat.fname; 
                mat0                   = dat.modality{m}.nii(n).mat; 
                spm_get_space(fname,mat1\mat0);  
                dat.modality{m}.nii(n) = nifti(fname);                      
            end                           
        end
    end % Loop over modalities    
                    
    if isfield(dat,'label')
        % Do labels too
        R = numel(dat.label);
        for r=1:R          
            fname = dat.label{r}.nii.dat.fname;

            mat0             = dat.label{r}.nii.mat;
            spm_get_space(fname,mat1\mat0); 
            dat.label{r}.nii = nifti(fname); 
        end
    end
                    
    fprintf('done!\n')
end               

% Remove image data outside of the head
%--------------------------------------------------------------------------
if preproc.do_crop    
    fprintf('Removing air voxels... ')
    
    for m=1:M
        if isfield(dat.modality{m},'channel')
            C = numel(dat.modality{m}.channel);
            for c=1:C
                channel = dat.modality{m}.channel{c}.name;
                
                N = numel(dat.modality{m}.channel{c}.nii);
                for n=1:N
                    fname = dat.modality{m}.channel{c}.nii(n).dat.fname;

                    if strcmpi(channel,'IR')
                        spm_impreproc('atlas_crop',fname,'cr_',2); 
                    else
                        spm_impreproc('atlas_crop',fname,'cr_',preproc.do_rem_neck); 
                    end                                       

                    [dat.modality{m}.channel{c}.nii(n),nfname] = update_nii(fname,'cr_');

                    if isfield(dat.modality{m}.channel{c},'json')
                        spm_json_manager('modify_json_field',dat.modality{m}.channel{c}.json(n).pth,'pth',nfname);
                    end
                end                
            end            
        else
            N = numel(dat.modality{m}.nii);
            for n=1:N
                fname = dat.modality{m}.nii(n).dat.fname;

                spm_impreproc('atlas_crop',fname,'cr_',preproc.do_rem_neck); 

                [dat.modality{m}.nii(n),nfname] = update_nii(fname,'cr_');

                if isfield(dat.modality{m},'json')
                    spm_json_manager('modify_json_field',dat.modality{m}.json(n).pth,'pth',nfname);
                end
            end
        end  
    end    
    
    fprintf('done!\n')
end  

% In the case of IR data, remove bottom part of image
%--------------------------------------------------------------------------  
if ~preproc.do_crop  
    for m=1:M    
        if isfield(dat.modality{m},'channel')
            C = numel(dat.modality{m}.channel);
            for c=1:C
                channel = dat.modality{m}.channel{c}.name;

                N = numel(dat.modality{m}.channel{c}.nii);
                for n=1:N
                    fname = dat.modality{m}.channel{c}.nii(n).dat.fname;

                    if strcmpi(channel,'IR')
                        fprintf('Removing bottom of IR image... ')

                        spm_impreproc('atlas_crop',fname,'cr_',2); 

                        fprintf('done!\n')
                    else
                        continue;
                    end                                       

                    [dat.modality{m}.channel{c}.nii(n),nfname] = update_nii(fname,'cr_');

                    if isfield(dat.modality{m}.channel{c},'json')
                        spm_json_manager('modify_json_field',dat.modality{m}.channel{c}.json(n).pth,'pth',nfname);
                    end
                end                
            end 
        end
    end    
end

% Co-register images
%--------------------------------------------------------------------------
if preproc.do_coreg
    fprintf('Co-registering... ')
    
    % Perform co-registration    
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
                [V,res,source_ix] = spm_impreproc('coreg',V,ixhr);
                
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
        
    if isfield(dat,'label')
        % Do labels too
        R = numel(dat.label);
        for r=1:R          
            fname       = dat.label{r}.nii.dat.fname;
            nam_channel = dat.label{r}.nam_channel;
            mat0        = dat.label{r}.nii.mat;
            
            for i=source_ix
                if strcmpi(dat.modality{m}.channel{i}.name,nam_channel)
                    mat1 = res{i}{1}.M;
                    break
                end
            end
            
            spm_get_space(fname,mat1\mat0); 
            dat.label{r}.nii = nifti(fname); 
        end
    end    
    clear V
    
    fprintf('done!\n')
end

% NN down-sampling in-plane
%--------------------------------------------------------------------------
if preproc.do_ds_inplane
    fprintf('Down-sampling in-plane... ')
    
    for m=1:M
        if isfield(dat.modality{m},'channel')
            C = numel(dat.modality{m}.channel);
            for c=1:C
                N = numel(dat.modality{m}.channel{c}.nii);
                for n=1:N
                    fname = dat.modality{m}.channel{c}.nii(n).dat.fname;

                    spm_impreproc('downsample_inplane',fname);

                    [dat.modality{m}.channel{c}.nii(n),nfname] = update_nii(fname,'ds_');

                    if isfield(dat.modality{m}.channel{c},'json')
                        spm_json_manager('modify_json_field',dat.modality{m}.channel{c}.json(n).pth,'pth',nfname);
                    end
                end
            end
        else
            N = numel(dat.modality{m}.nii);
            for n=1:N
                fname = dat.modality{m}.nii(n).dat.fname;

                spm_impreproc('downsample_inplane',fname);

                [dat.modality{m}.nii(n),nfname] = update_nii(fname,'ds_');

                if isfield(dat.modality{m},'json')
                    spm_json_manager('modify_json_field',dat.modality{m}.json(n).pth,'pth',nfname);
                end
            end
        end  
    end     
    
    fprintf('done!\n')
end

% NN down-sampling through-plane
%--------------------------------------------------------------------------
if preproc.do_ds_throughplane
    fprintf('Down-sampling through-plane... ')
    
    for m=1:M
        if isfield(dat.modality{m},'channel')
            C = numel(dat.modality{m}.channel);
            for c=1:C
                N = numel(dat.modality{m}.channel{c}.nii);
                for n=1:N
                    fname = dat.modality{m}.channel{c}.nii(n).dat.fname;

                    spm_impreproc('downsample_throughplane',fname);

                    [dat.modality{m}.channel{c}.nii(n),nfname] = update_nii(fname,'dsz_');

                    if isfield(dat.modality{m}.channel{c},'json')
                        spm_json_manager('modify_json_field',dat.modality{m}.channel{c}.json(n).pth,'pth',nfname);
                    end
                end
            end
        else
            N = numel(dat.modality{m}.nii);
            for n=1:N
                fname = dat.modality{m}.nii(n).dat.fname;

                spm_impreproc('downsample_throughplane',fname);

                [dat.modality{m}.nii(n),nfname] = update_nii(fname,'dsz_');

                if isfield(dat.modality{m},'json')
                    spm_json_manager('modify_json_field',dat.modality{m}.json(n).pth,'pth',nfname);
                end
            end
        end  
    end
    
    fprintf('done!\n')
end

% Create equally sized images by super-resolution
%--------------------------------------------------------------------------
if preproc.do_superres    
    fprintf('Super-resolving... ')
    
    for m=1:M    
        modality = dat.modality{m}.name;
        if ~strcmpi(modality,'MRI'), continue; end
                 
        if isfield(dat.modality{m},'channel')                        
            C        = numel(dat.modality{m}.channel);
            channels = cell(1,C);
            img      = cell(1,C); % Input object to super-resolution algorithm 
            for c=1:C
                channels{c} = dat.modality{m}.channel{c}.name;
                img{c}      = nifti;
                N = numel(dat.modality{m}.channel{c}.nii);
                for n=1:N
                    fname     = dat.modality{m}.channel{c}.nii(n).dat.fname;
                    img{c}(n) = nifti(fname);
                end
            end
        else
            channels = 'NaN';
            
            img{1} = nifti;
            N      = numel(dat.modality{m}.nii);            
            for n=1:N
                fname     = dat.modality{m}.nii(n).dat.fname;
                img{1}(n) = nifti(fname);
            end
        end
        
        % Do super-resolution
        error('Add spm_preproc_mtv here!');          
        
        if isfield(dat.modality{m},'channel')
            C = numel(dat.modality{m}.channel);
            for c=1:C
                nfname                            = nii(c).dat.fname;
                dat.modality{m}.channel{c}.nii(1) = nii(c);
                
                if isfield(dat.modality{m}.channel{c},'json')
                    spm_json_manager('modify_json_field',dat.modality{m}.channel{c}.json(1).pth,'pth',nfname);
                end
                
                N = numel(dat.modality{m}.channel{c}.nii);
                for n=2:N                      
                    if isfield(dat.modality{m}.channel{c},'json')
                        delete(dat.modality{m}.channel{c}.json(n).pth);
                        dat.modality{m}.channel{c}.json(n) = [];
                    end
                    dat.modality{m}.channel{c}.nii(n)  = [];                    
                end
            end                        
        else
            nfname                 = nii(c).dat.fname;
            dat.modality{m}.nii(1) = nii(c);
            if isfield(dat.modality{m},'json')
                spm_json_manager('modify_json_field',dat.modality{m}.json(1).pth,'pth',nfname);
            end
            
            N = numel(dat.modality{m}.nii);
            for n=2:N                    
                if isfield(dat.modality{m},'json')
                    delete(dat.modality{m}.json(n).pth);
                    dat.modality{m}.json(n) = [];
                end
                dat.modality{m}.nii(n)  = [];                
            end
        end        
    end          
    
    fprintf('done!\n')
end                  

% Reslice
%--------------------------------------------------------------------------
if preproc.do_reslice && ~preproc.do_superres  
    fprintf('Reslicing... ')
    
    V      = spm_vol;
    cnt    = 0;
    for m=1:M
        if isfield(dat.modality{m},'channel')
            C        = numel(dat.modality{m}.channel);
            chn_nams = cell(1,C);
            for c=1:C
                N           = numel(dat.modality{m}.channel{c}.nii);
                chn_nams{c} = dat.modality{m}.channel{c}.name;
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
        V = spm_impreproc('reslice',V,deg,ixhr); 
        
        cnt = 0;
        for m=1:M
            if isfield(dat.modality{m},'channel')
                C = numel(dat.modality{m}.channel);
                for c=1:C
                    N = numel(dat.modality{m}.channel{c}.nii);
                    for n=1:N
                        cnt                               = cnt + 1;                                           
                        dat.modality{m}.channel{c}.nii(n) = nifti(V(cnt).fname);  
                        
                        if isfield(dat.modality{m}.channel{c},'json')
                            spm_json_manager('modify_json_field',dat.modality{m}.channel{c}.json(n).pth,'pth',V(cnt).fname);
                        end
                    end
                end
            else
                N = numel(dat.modality{m}.nii);
                for n=1:N
                    cnt                    = cnt + 1;                                           
                    dat.modality{m}.nii(n) = nifti(V(cnt).fname);
                    
                    if isfield(dat.modality{m},'json')
                        spm_json_manager('modify_json_field',dat.modality{m}.json(n).pth,'pth',V(cnt).fname);
                    end
                end
            end  
        end
        
    end
    clear V
        
    fprintf('done!\n')
end

% Change voxel size of image(s)
%--------------------------------------------------------------------------
if ~isempty(preproc.vx) && ~preproc.do_superres
    fprintf('Changing voxel sizes... ')
    
    for m=1:M        
        if isfield(dat.modality{m},'channel')
            C = numel(dat.modality{m}.channel);
            for c=1:C
                N = numel(dat.modality{m}.channel{c}.nii);
                for n=1:N
                    fname = dat.modality{m}.channel{c}.nii(n).dat.fname;

                    spm_impreproc('nm_reorient',fname,preproc.vx,'vx_',deg);  

                    [dat.modality{m}.channel{c}.nii(n),nfname] = update_nii(fname,'vx_');

                    if isfield(dat.modality{m}.channel{c},'json')
                        spm_json_manager('modify_json_field',dat.modality{m}.channel{c}.json(n).pth,'pth',nfname);
                    end
                end
            end
        else
            N = numel(dat.modality{m}.nii);
            for n=1:N
                fname = dat.modality{m}.nii(n).dat.fname;

                spm_impreproc('nm_reorient',fname,preproc.vx,'vx_',deg);  

                [dat.modality{m}.nii(n),nfname] = update_nii(fname,'vx_');

                if isfield(dat.modality{m},'json')
                    spm_json_manager('modify_json_field',dat.modality{m}.json(n).pth,'pth',nfname);
                end
            end
        end       
    end
    
    fprintf('done!\n')
end

% Resize
%--------------------------------------------------------------------------
if preproc.resize
    fprintf('Resizing... ')
    
    for m=1:M        
        if isfield(dat.modality{m},'channel')
            C = numel(dat.modality{m}.channel);
            for c=1:C
                N = numel(dat.modality{m}.channel{c}.nii);
                for n=1:N
                    fname = dat.modality{m}.channel{c}.nii(n).dat.fname;

                    resize_nii(dat.modality{m}.channel{c}.nii(n)); 

                    [dat.modality{m}.channel{c}.nii(n),nfname] = update_nii(fname,'r_');

                    if isfield(dat.modality{m}.channel{c},'json')
                        spm_json_manager('modify_json_field',dat.modality{m}.channel{c}.json(n).pth,'pth',nfname);
                    end
                end
            end
        else
            N = numel(dat.modality{m}.nii);
            for n=1:N
                fname = dat.modality{m}.nii(n).dat.fname;

                resize_nii(dat.modality{m}.nii(n)); 

                [dat.modality{m}.nii(n),nfname] = update_nii(fname,'r_');

                if isfield(dat.modality{m},'json')
                    spm_json_manager('modify_json_field',dat.modality{m}.json(n).pth,'pth',nfname);
                end
            end
        end       
    end   
    
    fprintf('done!\n')       
end

% Reslice labels
%--------------------------------------------------------------------------
if isfield(dat,'label')
    fprintf('Reslicing labels... ')

    % Create spm_vol of reference volume
    V   = spm_vol;
    cnt = 0;
    for m=1:M
        if isfield(dat.modality{m},'channel')
            C        = numel(dat.modality{m}.channel);
            chn_nams = cell(1,C);
            for c=1:C
                N           = numel(dat.modality{m}.channel{c}.nii);
                chn_nams{c} = dat.modality{m}.channel{c}.name;
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
   
    % Reslice labels
    R = numel(dat.label);
    for r=1:R          
        nii = reslice_labels(V(ixhr),dat.label{r}.nii);

        dat.label{r}.nii = nii;                 

        if isfield(dat.label{r},'json')
            spm_json_manager('modify_json_field',dat.label{r}.json.pth,'pth',nii.dat.fname);
        end
    end
        
    fprintf('done!\n')
end

% Denoise images
%--------------------------------------------------------------------------
if preproc.do_denoise
    fprintf('Denoising... ')
    
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
        nii_in = nii;
        nii    = spm_mtv_preproc('InputImages',nii_in,'Modality',modality,'Verbose',0,'OutputDirectory',fileparts(nii(1).dat.fname));
        for i=1:numel(nii)
            delete(nii_in(i).dat.fname);
        end
        
        cnt = 0;
        if isfield(dat.modality{m},'channel')
            C = numel(dat.modality{m}.channel);
            for c=1:C
                N = numel(dat.modality{m}.channel{c}.nii);
                for n=1:N
                    cnt                               = cnt + 1;
                    dat.modality{m}.channel{c}.nii(n) = nii(cnt);   
                    
                    nfname = nii(cnt).dat.fname;
                    
                    if isfield(dat.modality{m}.channel{c},'json')
                        spm_json_manager('modify_json_field',dat.modality{m}.channel{c}.json(n).pth,'pth',nfname);
                    end
                end
            end
        else
            N = numel(dat.modality{m}.nii);
            for n=1:N
                cnt                    = cnt + 1;
                dat.modality{m}.nii(n) = nii(cnt);     
                
                nfname = nii(cnt).dat.fname;
                
                if isfield(dat.modality{m},'json')
                    spm_json_manager('modify_json_field',dat.modality{m}.json(n).pth,'pth',nfname);
                end
            end
        end        
    end    
    clear nii
    
    fprintf('done!\n')
end

% Simple normalisation of image intensities (make mean(img)=100)
%--------------------------------------------------------------------------
if preproc.do_normalise_intensities    
    fprintf('Intensity normalising... ')
    
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
    
    fprintf('done!\n')
end   

% Segment using spm_preproc8
%--------------------------------------------------------------------------
if preproc.do_segment || preproc.do_skull_strip || preproc.do_bf_correct              
    fprintf('Segmenting... ')

    % For segmenting we assume... 
    m = 1; % ...one modality...
    n = 1; % ...and one image per channel        
    
    % Set-up output
    write_tc = false(6,4); % c, rc, wc, mwc    
    write_bf = [false false]; 
    write_df = false(3,1);                       
    
    write_tc([1 2 3],1) = true;
    
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
    segment_preproc8(V,write_tc,write_bf,write_df,dir_preproc,dat.modality{m}.name);
     
    % Get output files
    dat = add_segmentations2dat(dat);
    
    fprintf('done!\n')
end         

% Overwrite image data with bias-corrected version from spm_preproc8
%--------------------------------------------------------------------------
if preproc.do_bf_correct
    fprintf('Bias-field correcting... ')
    
    if isfield(dat.modality{m},'channel')            
        C = numel(dat.modality{m}.channel);
        for c=1:C            
            fname = dat.modality{m}.channel{c}.nii(n).dat.fname;

            [dat.modality{m}.channel{c}.nii(n),nfname] = update_nii(fname,'m');

            if isfield(dat.modality{m}.channel{c},'json')
                spm_json_manager('modify_json_field',dat.modality{m}.channel{c}.json(n).pth,'pth',nfname);
            end
        end
    else        
        fname = dat.modality{m}.nii(n).dat.fname;

        [dat.modality{m}.nii(n),nfname] = update_nii(fname,'m');

        if isfield(dat.modality{m},'json')
            spm_json_manager('modify_json_field',dat.modality{m}.json(n).pth,'pth',nfname);
        end
    end        
    
    fprintf('done!\n')
end

% Skull-strip image data using mask built from GM, WM and CSF classes from
% spm_preproc8 segmentation
%--------------------------------------------------------------------------
if preproc.do_skull_strip  
    fprintf('Skull-stripping... ')            
    
    dat = skull_strip(dat);                 
    
    fprintf('done!\n')
end  

% Clean-up
%--------------------------------------------------------------------------
if ~preproc.do_segment && (preproc.do_skull_strip || preproc.do_bf_correct)
    fprintf('Cleaning-up... ')
    
    for t=1:numel(dat.segmentation)
        for k=1:numel(dat.segmentation{t}.class)
            delete(dat.segmentation{t}.class{k}.nii.dat.fname);
            
            if isfield(dat.segmentation{t}.class{k},'json')
                delete(dat.segmentation{t}.class{k}.json.pth);
            end
        end
    end
    
    dat = rmfield(dat,'segmentation');
    dat = rmfield(dat,'segmentation_map');
    
    fprintf('done!\n')
end

% Create 2D versions
%--------------------------------------------------------------------------
if write_2d
    fprintf('Writing 2D data... ')
    
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

                    if isfield(dat.modality{m}.channel{c},'json')
                        fname       = dat.modality{m}.channel{c}.json(n).pth;
                        [~,nam,ext] = fileparts(fname);
                        nfname      = fullfile(dir_2d,[nam ext]);                
                        copyfile(fname,nfname);   

                        dat.modality{m}.channel{c}.json(n).pth = nfname;

                        spm_json_manager('modify_json_field',dat.modality{m}.channel{c}.json(n).pth,'pth',nfname);
                    end
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

                if isfield(dat.modality{m},'json')
                    fname       = dat.modality{m}.json(n).pth;
                    [~,nam,ext] = fileparts(fname);
                    nfname      = fullfile(dir_2d,[nam ext]);            
                    copyfile(fname,nfname);

                    dat.modality{m}.json(n).pth = nfname;

                    spm_json_manager('modify_json_field',dat.modality{m}.json(n).pth,'pth',nfname);
                end
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

            if isfield(dat.label{r},'json')
                fname       = dat.label{r}.json.pth;
                [~,nam,ext] = fileparts(fname);
                nfname1     = fullfile(dir_2d,[nam ext]);            
                copyfile(fname,nfname1);

                dat.label{r}.json.pth = nfname1;

                spm_json_manager('modify_json_field',dat.label{r}.json.pth,'pth',nfname);
            end
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
               
                if isfield(dat.segmentation{t}.class{k},'json')
                    fname       = dat.segmentation{t}.class{k}.json.pth;
                    [~,nam,ext] = fileparts(fname);
                    nfname1     = fullfile(dir_2d,[nam ext]);            
                    copyfile(fname,nfname1);

                    dat.segmentation{t}.class{k}.json.pth = nfname1;

                    spm_json_manager('modify_json_field',dat.segmentation{t}.class{k}.json.pth,'pth',nfname);
                end
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

                    nfname = spm_imbasics('create_2d_slice',fname,deg,axis_2d);

                    dat.modality{m}.channel{c}.nii(n) = nifti(nfname);

                    if isfield(dat.modality{m}.channel{c},'json')
                        spm_json_manager('modify_json_field',dat.modality{m}.channel{c}.json(n).pth,'pth',nfname);
                    end
                end
            end
        else
            N = numel(dat.modality{m}.nii);
            for n=1:N
                fname = dat.modality{m}.nii(n).dat.fname;

                nfname = spm_imbasics('create_2d_slice',fname,deg,axis_2d);

                dat.modality{m}.nii(n) = nifti(nfname);

                if isfield(dat.modality{m},'json')
                    spm_json_manager('modify_json_field',dat.modality{m}.json(n).pth,'pth',nfname);
                end
            end
        end   
    end
    
    if isfield(dat,'label')
        R = numel(dat.label);
        for r=1:R
            fname  = dat.label{r}.nii.dat.fname;
            
            nfname = spm_imbasics('create_2d_slice',fname,deg,axis_2d);

            dat.label{r}.nii = nifti(nfname);

            if isfield(dat.label{r},'json')
                spm_json_manager('modify_json_field',dat.label{r}.json.pth,'pth',nfname);
            end
        end
    end    
    
    if isfield(dat,'segmentation')
        T = numel(dat.segmentation);
        for t=1:T
            K = numel(dat.segmentation{t}.class);
            for k=1:K
                fname = dat.segmentation{t}.class{k}.nii.dat.fname;
                
                nfname = spm_imbasics('create_2d_slice',fname,deg,axis_2d);                
                
                dat.segmentation{t}.class{k}.nii = nifti(nfname);

                if isfield(dat.segmentation{t}.class{k},'json')
                    spm_json_manager('modify_json_field',dat.segmentation{t}.class{k}.json.pth,'pth',nfname);
                end
           end
        end
    end    
    
    fprintf('done!\n')
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
if nargin<3, val = 1000; end

msk  = spm_misc('msk_modality',img,modality);     
sint = sum(reshape(img(msk),[],1));
nm   = nnz(msk);        
scl  = (val/(sint/nm));
img  = scl*img;
%==========================================================================   

%==========================================================================   
function [ixc,ixn,ref_ix] = get_ix_hr(modality)
prod_vx = [];
if isfield(modality,'channel')
    C = numel(modality.channel);
    for c=1:C
        N = numel(modality.channel{c}.nii);
        for n=1:N
            vx      = spm_misc('vxsize',modality.channel{c}.nii(n).mat);
            prod_vx = [prod_vx prod(vx)];
        end
    end
    [~,ref_ix] = min(prod_vx);

    cnt = 0;
    for c=1:C
        N = numel(modality.channel{c}.nii);
        for n=1:N
            cnt = cnt + 1;
            if cnt==ref_ix
                ixc = c;
                ixn = n;
            end
        end
    end
else
    ixc = 1;
    N = numel(modality.nii);
    for n=1:N
        vx      = spm_misc('vxsize',modality.nii(n).mat);
        prod_vx = [prod_vx prod(vx)];
    end
    [~,ixn] = min(prod_vx);
    ref_ix  = ixn;
end
%==========================================================================   
