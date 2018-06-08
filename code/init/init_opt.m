function opt = init_opt(job,dir_preproc,dir_2d)
opt = struct;

opt.dir_preproc = dir_preproc;
opt.dir_2d      = dir_2d;
opt.axis_2d     = job.axis_2d;
opt.write_2d    = job.write_2d;
opt.preproc     = job.preproc;
%==========================================================================

