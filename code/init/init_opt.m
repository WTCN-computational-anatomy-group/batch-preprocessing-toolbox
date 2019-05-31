function opt = init_opt(job,dir_preproc,dir_2d)
if nargin < 3, dir_2d = ''; end
    
opt.dir_preproc = dir_preproc;
opt.dir_2d      = dir_2d;
opt.axis_2d     = job.axis_2d;
opt.write_2d    = job.write_2d;
opt.preproc     = job.preproc;
opt.segment     = job.segment;
%==========================================================================

