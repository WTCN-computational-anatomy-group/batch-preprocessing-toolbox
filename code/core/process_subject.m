function obj = process_subject(obj) 

% Make copies into obj.dir_preproc (in order to not modify the original
% data)
%--------------------------------------------------------------------------
fname = obj.scans{1}{1}.fname;
dir0  = fileparts(fname);
dir0  = strsplit(dir0,filesep);
dir0  = fullfile(obj.dir_preproc,dir0{end - 2});
mkdir(dir0);     

dir_scans = fullfile(dir0,'scans');
mkdir(dir_scans);     

dir_labels = fullfile(dir0,'labels');
mkdir(dir_labels);        

% Copy scans
N = numel(obj.scans);
for n=1:N
    I = numel(obj.scans{n});
    for i=1:I
        fname = obj.scans{n}{i}.fname;
        dir1  = fileparts(fname);
        dir1  = strsplit(dir1,filesep);
        dir1  = fullfile(dir_scans,dir1{end});
        mkdir(dir1);  
        
        copyfile(fname,dir1);
        [~,nam,ext]           = fileparts(fname);
        nfname                = fullfile(dir1,[nam ext]);
        obj.scans{n}{i}.fname = nfname;
    end
end        

if ~isempty(obj.labels)
    % Copy labels        
    copyfile(obj.labels.fname,dir_labels);
    [~,nam,ext]      = fileparts(obj.labels.fname);
    nfname           = fullfile(dir_labels,[nam ext]);
    obj.labels.fname = nfname;
end

% Remove image data outside of the head (air..)
%--------------------------------------------------------------------------
if obj.preproc.do_crop
    for n=1:N
        I = numel(obj.scans{n});
        for i=1:I
            V = obj.scans{n}{i};
            
            spm_impreproc('atlas_crop',V.fname,'cr_',obj.preproc.do_rem_neck); 

            [pth,nam,ext] = fileparts(V.fname);
            delete(V.fname);
            nfname        = fullfile(pth,['cr_' nam ext]);
            V             = spm_vol(nfname);                                  
            
            obj.scans{n}{i} = V;
        end
    end    
end  

% Denoise images
%--------------------------------------------------------------------------
if obj.preproc.do_denoise
    obj = spm_denoise(obj);
end

% Create equally sized images by super-resolution
%--------------------------------------------------------------------------
if obj.preproc.do_superres    
    obj = spm_superres(obj);
end                  

% Reslice to size of image with largest FOV
%--------------------------------------------------------------------------
if obj.preproc.do_reslice && N>1        
    for n=1:N
        V(n) = obj.scans{n}{1};   
    end  
    
    V = spm_impreproc('reslice',V);    
    
    for n=1:N
        obj.scans{n}{1} = V(n);
    end  
end

% Change voxel size of image(s)
%--------------------------------------------------------------------------
if ~isempty(obj.preproc.vx)
    for n=1:N
        V = obj.scans{n}{1};   
    
        spm_impreproc('nm_reorient',V.fname,obj.preproc.vx,1,'vx_');    
    
        [pth,nam,ext]   = fileparts(V.fname);
        delete(V.fname);
        nfname          = fullfile(pth,['vx_' nam ext]);
        obj.scans{n}{1} = spm_vol(nfname);
    end  
end

% Rigidly realign to MNI space
%--------------------------------------------------------------------------
if obj.preproc.do_realign2mni  
    % Just align the first image
    V               = obj.scans{1}{1};            
    M               = spm_impreproc('rigid_align',V.fname);                     
    obj.scans{1}{1} = spm_vol(V.fname);    
         
    % Then change orientation matrices of the rest accordingly
    for n=1:N
        I = numel(obj.scans{n});
        for i=1:I            
            if n==1 && i==1
                continue;
            end
            
            V = obj.scans{n}{i};
            spm_get_space(V.fname,M\V.mat);  
            obj.scans{n}{i} = spm_vol(V.fname);  
        end
    end    
end               

% Co-register images
%--------------------------------------------------------------------------
if obj.preproc.do_coreg && N>1
    cnt = 1;
    for n=1:N
        I = numel(obj.scans{n});
        for i=1:I
            V(cnt) = obj.scans{n}{i};            
            cnt    = cnt + 1;
        end
    end  
    
    V = spm_impreproc('coreg',V);
    
    cnt = 1;
    for n=1:N
        I = numel(obj.scans{n});
        for i=1:I
            obj.scans{n}{i} = V(cnt);            
            cnt             = cnt + 1;
        end
    end  
end

% Segment and/or bias-field correct and/or skull-strip
%--------------------------------------------------------------------------
if any(any(obj.preproc.write_tc==true) == true) || obj.preproc.do_skull_strip || obj.preproc.do_bf_correct    
    dir_seg = fullfile(dir0,'segmentations');
    mkdir(dir_seg);      
    
    for n=1:N
        V(n) = obj.scans{n}{1};   
    end  
    
    write_tc = obj.preproc.write_tc;
    write_bf = obj.preproc.write_bf;
    write_df = obj.preproc.write_df;
    
    if obj.preproc.do_bf_correct
        write_bf(1,2) = true;
    end
    
    if obj.preproc.do_skull_strip
        write_tc(1:3,1) = true;
    end
    
    segment_subject(V,write_tc,write_bf,write_df,dir_seg);
       
    if obj.preproc.do_bf_correct
        % Overwrite image data with bias-corrected versions
        files = spm_select('FPList',dir_seg,'^m.*\.nii$');
        for n=1:N          
            fname   = obj.scans{n}{1}.fname;            
            Nii     = nifti(fname);
            [~,nam] = fileparts(fname);
            for n1=1:N
                [~,nam_bf] = fileparts(files(n1,:));
                if strcmp(nam,nam_bf(2:end))
                    Nii_bf         = nifti(files(n1,:));
                    Nii.dat(:,:,:) = Nii_bf.dat(:,:,:);
                end
            end  
        end
    end
    
    if obj.preproc.do_skull_strip
        % Overwrite image data with skull-stripped versions
        files = spm_select('FPList',dir_seg,'^c[1-3].*\.nii$');
        V0    = spm_vol(files);
        K     = numel(V0);
        msk   = zeros(V0(1).dim,'single');
        for k=1:K
            Nii  = nifti(V0(k).fname);
            resp = single(Nii.dat(:,:,:)); 
            msk  = msk + resp;
        end
        clear resp

        spm_imbasics('smooth_img_in_mem',msk,10);      
        msk = msk>0.5;    

        % Fill holes
        msk = imfill(msk,6,'holes');    
        
        if 0
            % For testing
            split = 6;
            dm0   = V0(1).dim;
            nfigs = floor(dm0(3)/split);
            
            F1 = floor(sqrt(nfigs));
            F2 = ceil(nfigs/F1);      
            figure(666); 
            for f=1:nfigs
                subplot(F1,F2,f);            
                imagesc(msk(:,:,split*f)'); colormap(gray); axis off xy image;
            end
        end
        
        for n=1:N          
            fname          = obj.scans{n}{1}.fname;            
            Nii            = nifti(fname);
            img            = single(Nii.dat(:,:,:));
            img(~msk)      = NaN;               
            Nii.dat(:,:,:) = img;
        end  
    end  
    
    if obj.preproc.make_ml_labels && isempty(obj.labels)
        % Write maximum-likelihoods labels (only if labels are not available)                
        files = spm_select('FPList',dir_seg,'^c.*\.nii$');
        V0    = spm_vol(files);
        K     = numel(V0);
        img   = zeros([V0(1).dim K],'single');
        for k=1:K
            Nii          = nifti(V0(k).fname);
            img(:,:,:,k) = single(Nii.dat(:,:,:));
        end
                
        if K<6
            % Less than the default SPM template number of classes requested
            % -> correct ML labels
            img1 = ones(V0(1).dim,'single');
            img1 = img1 - sum(img,4);
            img  = cat(4,img,img1);
            clear img1
        end
        
        [~,ml_labels] = max(img,[],4);        
        clear img
         
        fname       = obj.scans{1}{1}.fname;  
        [~,nam,ext] = fileparts(fname);
        nfname      = fullfile(dir_labels,['ml-' nam ext]);
        
        Nii      = nifti;
        Nii.dat  = file_array(nfname,size(ml_labels),'uint8',0,1/K,0);
        Nii.mat  = V0(1).mat;
        Nii.mat0 = V0(1).mat;
        Nii.descrip = 'ML-labels';
        create(Nii);
        
        Nii.dat(:,:,:) = ml_labels;
        clear ml_labels
        
        obj.labels = spm_vol(nfname);
    end
    
    % Clean-up
    if ~any(any(obj.preproc.write_tc==true) == true)
        rmdir(dir_seg,'s');
    end
end         

% Create 2D versions
%--------------------------------------------------------------------------
if obj.preproc.write_2d
    fname = obj.scans{1}{1}.fname;
    dir0  = fileparts(fname);
    dir0  = strsplit(dir0,filesep);
    dir0  = fullfile(obj.dir_preproc_2d,dir0{end - 2});
    mkdir(dir0);     

    dir_scans = fullfile(dir0,'scans');
    mkdir(dir_scans);     

    dir_labels = fullfile(dir0,'labels');
    mkdir(dir_labels);        

    % Of scans    
    N = numel(obj.scans);
    for n=1:N
        I = numel(obj.scans{n});
        for i=1:I
            fname = obj.scans{n}{i}.fname;
            dir1  = fileparts(fname);
            dir1  = strsplit(dir1,filesep);
            dir1  = fullfile(dir_scans,dir1{end});
            mkdir(dir1);  

            copyfile(fname,dir1);
            [~,nam,ext] = fileparts(fname);
            nfname      = fullfile(dir1,[nam ext]);
            
            create_2d_slice(nfname,obj.preproc.axis_2d);
        end
    end        

    if ~isempty(obj.labels)
        % Of labels        
        copyfile(obj.labels.fname,dir_labels);
        [~,nam,ext] = fileparts(obj.labels.fname);
        nfname      = fullfile(dir_labels,[nam ext]);

        create_2d_slice(nfname,obj.preproc.axis_2d);     
    end
    
    if any(any(obj.preproc.write_tc==true) == true)                
        % Of segmentations                 
        dir_seg1 = fullfile(dir0,'segmentations');
        mkdir(dir_seg1);                                   
        
        prefix = {'c','wc','mwc'};
        for i=1:numel(prefix)
            files = spm_select('FPList',dir_seg,['^' prefix{i} '.*\.nii$']);
            
            for i1=1:size(files,1)
                copyfile(files(i1,:),dir_seg1);
                [~,nam,ext] = fileparts(files(i1,:));
                nfname      = fullfile(dir_seg1,[nam ext]);

                create_2d_slice(nfname,obj.preproc.axis_2d);            
            end
        end
    end
end      
%==========================================================================

