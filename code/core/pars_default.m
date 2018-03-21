function pars = pars_default(pars,test_level)
if nargin<2, test_level = 0; end

% Data-set specific parameters (m=1,...,M)
%--------------------------------------------------------------------------
if ~isfield(pars,'dat'), 
    error('pars.dat needs to be defined!'); 
end

M = numel(pars.dat);
for m=1:M
    if ~isfield(pars.dat{m},'dir_data'), 
        error('pars.dat.dir_data needs to be defined!'); 
    end

    % General parameters
    %----------------------------------------------------------------------    
    if ~isfield(pars.dat{m},'S')
        pars.dat{m}.S = Inf;
        if test_level==2 || test_level==3, pars.dat{m}.S = 16;
        elseif test_level==1               pars.dat{m}.S = 1;   
        end 
    end    
    if ~isfield(pars.dat{m},'modality')
        pars.dat{m}.modality = 'MRI';
    end    
    if ~isfield(pars.dat{m},'dir_preproc')
        pars.dat{m}.dir_preproc = '';
    end    
    
    % Pre-processing parameters
    %----------------------------------------------------------------------
    if ~isfield(pars.dat{m},'preproc')
        pars.dat{m}.preproc = struct;
    end
    if ~isfield(pars.dat{m}.preproc,'tol_dist')
        pars.dat{m}.preproc.tol_dist = 4;
    end
    if ~isfield(pars.dat{m}.preproc,'tol_vx')
        pars.dat{m}.preproc.tol_vx = 5;
    end
    if ~isfield(pars.dat{m}.preproc,'do_coreg')
        pars.dat{m}.preproc.do_coreg = false;
    end
    if ~isfield(pars.dat{m}.preproc,'do_reslice')
        pars.dat{m}.preproc.do_reslice = false;
    end    
    if ~isfield(pars.dat{m}.preproc,'do_realign2mni')
        pars.dat{m}.preproc.do_realign2mni = false;
    end
    if ~isfield(pars.dat{m}.preproc,'do_crop')
        pars.dat{m}.preproc.do_crop = false;
    end
    if ~isfield(pars.dat{m}.preproc,'do_rem_neck')
        pars.dat{m}.preproc.do_rem_neck = false;
    end
    if ~isfield(pars.dat{m}.preproc,'do_skull_strip')
        pars.dat{m}.preproc.do_skull_strip = false;
    end
    if ~isfield(pars.dat{m}.preproc,'do_superres')
        pars.dat{m}.preproc.do_superres = false;
    end    
    if ~isfield(pars.dat{m}.preproc,'do_denoise')
        pars.dat{m}.preproc.do_denoise = false;
    end
    if ~isfield(pars.dat{m}.preproc,'write_2d')
        pars.dat{m}.preproc.write_2d = false;
    end
    if ~isfield(pars.dat{m}.preproc,'axis_2d')
        pars.dat{m}.preproc.axis_2d = 3;
    end    
    if ~isfield(pars.dat{m}.preproc,'vx')
        pars.dat{m}.preproc.vx = [];
    end
    if ~isfield(pars.dat{m}.preproc,'write_tc')
        pars.dat{m}.preproc.write_tc = false(6,4);        
    end    
    if ~isfield(pars.dat{m}.preproc,'do_bf_correct')
        pars.dat{m}.preproc.do_bf_correct = false;        
    end    
    if ~isfield(pars.dat{m}.preproc,'write_bf')
        pars.dat{m}.preproc.write_bf = false(1,2);
    end    
    if ~isfield(pars.dat{m}.preproc,'write_df')
        pars.dat{m}.preproc.write_df = false(1,2);
    end    
    if ~isfield(pars.dat{m}.preproc,'make_ml_labels')
        pars.dat{m}.preproc.make_ml_labels = false;
    end 
    
    % Super-resolution parameters
    %----------------------------------------------------------------------
    if ~isfield(pars.dat{m}.preproc,'superres')
        pars.dat{m}.preproc.superres = struct;
    end
    if ~isfield(pars.dat{m}.preproc.superres,'verbose')
        pars.dat{m}.preproc.superres.verbose = false;
    end        
    if ~isfield(pars.dat{m}.preproc.superres,'vx')
        pars.dat{m}.preproc.superres.vx = [1 1 1];
    end    
    if ~isfield(pars.dat{m}.preproc.superres,'proj_mat')
        pars.dat{m}.preproc.superres.proj_mat = 'sinc';
    end    
    if ~isfield(pars.dat{m}.preproc.superres,'trunc_ct')
        pars.dat{m}.preproc.superres.trunc_ct = [];
    end  
    
    if ~isfield(pars.dat{m}.preproc.superres,'admm')
        pars.dat{m}.preproc.superres.admm = struct;
    end
    if ~isfield(pars.dat{m}.preproc.superres.admm,'rho')
        pars.dat{m}.preproc.superres.admm.rho = 1e5;
    end        
    if ~isfield(pars.dat{m}.preproc.superres.admm,'niter')
        pars.dat{m}.preproc.superres.admm.niter = 50;
    end       
    if ~isfield(pars.dat{m}.preproc.superres.admm,'tol')
        pars.dat{m}.preproc.superres.admm.tol = 1e-4;
    end  
    if ~isfield(pars.dat{m}.preproc.superres.admm,'verbose')
        pars.dat{m}.preproc.superres.admm.verbose = false;
    end  
    if ~isfield(pars.dat{m}.preproc.superres.admm,'mu')
        pars.dat{m}.preproc.superres.admm.mu = 10;
    end  
    if ~isfield(pars.dat{m}.preproc.superres.admm,'alpha')
        pars.dat{m}.preproc.superres.admm.alpha = 2;
    end     
    if ~isfield(pars.dat{m}.preproc.superres.admm,'cgs_niter')
        pars.dat{m}.preproc.superres.admm.cgs_niter = 10;
    end  
    if ~isfield(pars.dat{m}.preproc.superres.admm,'cgs_tol')
        pars.dat{m}.preproc.superres.admm.cgs_tol = 1e-3;
    end   
    if ~isfield(pars.dat{m}.preproc.superres.admm,'est_rho')
        pars.dat{m}.preproc.superres.admm.est_rho = true;
    end       
    
    % Denoising parameters
    %----------------------------------------------------------------------
    if ~isfield(pars.dat{m}.preproc,'denoise')
        pars.dat{m}.preproc.denoise = struct;
    end
    if ~isfield(pars.dat{m}.preproc.denoise,'verbose')
        pars.dat{m}.preproc.denoise.verbose = false;
    end    
    
    if ~isfield(pars.dat{m}.preproc.denoise,'admm')
        pars.dat{m}.preproc.denoise.admm = struct;
    end
    if ~isfield(pars.dat{m}.preproc.denoise.admm,'rho')
        pars.dat{m}.preproc.denoise.admm.rho = 1e5;
    end        
    if ~isfield(pars.dat{m}.preproc.denoise.admm,'niter')
        pars.dat{m}.preproc.denoise.admm.niter = 100;
    end       
    if ~isfield(pars.dat{m}.preproc.denoise.admm,'tol')
        pars.dat{m}.preproc.denoise.admm.tol = 1e-4;
    end  
    if ~isfield(pars.dat{m}.preproc.denoise.admm,'verbose')
        pars.dat{m}.preproc.denoise.admm.verbose = false;
    end      
    if ~isfield(pars.dat{m}.preproc.denoise.admm,'mu')
        pars.dat{m}.preproc.denoise.admm.mu = 10;
    end  
    if ~isfield(pars.dat{m}.preproc.denoise.admm,'alpha')
        pars.dat{m}.preproc.denoise.admm.alpha = 2;
    end     
    if ~isfield(pars.dat{m}.preproc.denoise.admm,'est_rho')
        pars.dat{m}.preproc.denoise.admm.est_rho = true;
    end       
end
%==========================================================================