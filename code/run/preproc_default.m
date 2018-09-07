function [job,holly] = preproc_default(job,test_level)
if nargin<2, test_level = 0; end

if isstring(job) || ischar(job)
    % Parameters are in a JSON-file
    job = spm_jsonread(job);
end

% General parameters
%----------------------------------------------------------------------    
if ~isfield(job,'S')
    job.S = Inf;
elseif isempty(job.S)
    job.S = Inf;
elseif strcmpi(job.S,'inf')
    job.S = Inf;
end    
if test_level==1
    job.S = min(1,job.S); 
elseif test_level==2
    job.S = min(8,job.S); 
elseif test_level==3
    job.S = min(16,job.S); 
end     
if ~isfield(job,'dir_preproc')
    error('~isfield(pars,''dir_preproc'')')
end    
if ~isfield(job,'write_3d')
    job.write_3d = true;
end
if ~isfield(job,'write_2d')
    job.write_2d = false;
end
if job.write_2d && ~isfield(job,'dir_2d')
    error('~isfield(pars,''dir_2d'')')
elseif ~isfield(job,'dir_2d')
    job.dir_2d = '';
end   
if ~isfield(job,'axis_2d')
    job.axis_2d = 3;
end    

% Pre-processing parameters
%----------------------------------------------------------------------
if ~isfield(job,'preproc')
    job.preproc = struct;
end
if ~isfield(job.preproc,'tol_dist')
    job.preproc.tol_dist = 4;
end
if ~isfield(job.preproc,'tol_vx')
    job.preproc.tol_vx = 5;
end
if ~isfield(job.preproc,'reset_origin')
    job.preproc.reset_origin = false;
end
if ~isfield(job.preproc,'do_coreg')
    job.preproc.do_coreg = false;
end
if ~isfield(job.preproc,'do_reslice')
    job.preproc.do_reslice = false;
end    
if ~isfield(job.preproc,'do_ds_inplane')
    job.preproc.do_ds_inplane = false;
end    
if ~isfield(job.preproc,'do_ds_throughplane')
    job.preproc.do_ds_throughplane = false;
end    
if ~isfield(job.preproc,'do_realign2mni')
    job.preproc.do_realign2mni = false;
end
if ~isfield(job.preproc,'do_crop')
    job.preproc.do_crop = false;
end
if ~isfield(job.preproc,'do_rem_neck')
    job.preproc.do_rem_neck = 0;
end
if ~isfield(job.preproc,'do_superres')
    job.preproc.do_superres = false;
end    
if ~isfield(job.preproc,'do_denoise')
    job.preproc.do_denoise = false;
end    
if ~isfield(job.preproc,'vx')
    job.preproc.vx = [];
end
if ~isfield(job.preproc,'deg')
    job.preproc.deg = 0;
end
if ~isfield(job.preproc,'do_normalise_intensities')
    job.preproc.do_normalise_intensities = false;
end 
if ~isfield(job.preproc,'do_bf_correct')
    job.preproc.do_bf_correct = false;        
end    
if ~isfield(job.preproc,'do_skull_strip')
    job.preproc.do_skull_strip = false;
end
if ~isfield(job.preproc,'do_segment')
    job.preproc.do_segment = false;
end
if ~isfield(job.preproc,'part_labels')
    job.preproc.part_labels = {};
end

% Super-resolution parameters
%----------------------------------------------------------------------
if ~isfield(job.preproc,'superres')
    job.preproc.superres = struct;
end
if ~isfield(job.preproc.superres,'verbose')
    job.preproc.superres.verbose = false;
end        
if ~isfield(job.preproc.superres,'vx')
    job.preproc.superres.vx = [1 1 1];
end    
if ~isfield(job.preproc.superres,'proj_mat')
    job.preproc.superres.proj_mat = 'sinc';
end    
if ~isfield(job.preproc.superres,'no_lambda')
    job.preproc.superres.no_lambda = {'IR'};
end     

if ~isfield(job.preproc.superres,'admm')
    job.preproc.superres.admm = struct;
end
if ~isfield(job.preproc.superres.admm,'rho')
    job.preproc.superres.admm.rho = 1;
end        
if ~isfield(job.preproc.superres.admm,'niter')
    job.preproc.superres.admm.niter = 50;
end       
if ~isfield(job.preproc.superres.admm,'tol')
    job.preproc.superres.admm.tol = 1e-5;
end  
if ~isfield(job.preproc.superres.admm,'verbose')
    job.preproc.superres.admm.verbose = false;
end  
if ~isfield(job.preproc.superres.admm,'mu')
    job.preproc.superres.admm.mu = 10;
end  
if ~isfield(job.preproc.superres.admm,'alpha')
    job.preproc.superres.admm.alpha = 2;
end     
if ~isfield(job.preproc.superres.admm,'cgs_niter')
    job.preproc.superres.admm.cgs_niter = 20;
end  
if ~isfield(job.preproc.superres.admm,'cgs_tol')
    job.preproc.superres.admm.cgs_tol = 1e-3;
end   
if ~isfield(job.preproc.superres.admm,'est_rho')
    job.preproc.superres.admm.est_rho = false;
end       

% Denoising parameters
%----------------------------------------------------------------------
if ~isfield(job.preproc,'denoise')
    job.preproc.denoise = struct;
end
if ~isfield(job.preproc.denoise,'verbose')
    job.preproc.denoise.verbose = false;
end    
if ~isfield(job.preproc.denoise,'lambda_ct')
    job.preproc.denoise.lambda_ct = 0.3;
end    

if ~isfield(job.preproc.denoise,'admm')
    job.preproc.denoise.admm = struct;
end
if ~isfield(job.preproc.denoise.admm,'rho')
    job.preproc.denoise.admm.rho = 1;
end        
if ~isfield(job.preproc.denoise.admm,'niter')
    job.preproc.denoise.admm.niter = 50;
end       
if ~isfield(job.preproc.denoise.admm,'tol')
    job.preproc.denoise.admm.tol = 1e-5;
end  
if ~isfield(job.preproc.denoise.admm,'verbose')
    job.preproc.denoise.admm.verbose = false;
end      
if ~isfield(job.preproc.denoise.admm,'mu')
    job.preproc.denoise.admm.mu = 10;
end  
if ~isfield(job.preproc.denoise.admm,'alpha')
    job.preproc.denoise.admm.alpha = 2;
end     
if ~isfield(job.preproc.denoise.admm,'est_rho')
    job.preproc.denoise.admm.est_rho = false;
end   

% Holly stuff
%----------------------------------------------------------------------
if ~isfield(job,'holly')
    holly = struct;
else
    holly = job.holly;
end

if     test_level==1, holly.mode = 'for';
elseif test_level==2, holly.mode = 'parfor';
end

holly = distribute_default(holly);

if size(holly.translate,1)>1
    holly.translate = holly.translate';
end
%==========================================================================
