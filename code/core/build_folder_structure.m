function [dir_preproc,dir_2d] = build_folder_structure(job)
dir_preproc = job.dir_preproc;  
dir_preproc = strsplit(dir_preproc,filesep); 

nam = dir_preproc{end};    
nam = append_nam(nam,job);

dir_preproc = fullfile(dir_preproc{1:end - 1},nam);

if strcmp(job.dir_preproc(1),filesep)
    dir_preproc = [filesep dir_preproc];
end

if exist(dir_preproc,'dir') == 7, rmdir(dir_preproc,'s'); end; mkdir(dir_preproc); 

if job.write_2d
    dir_2d = job.dir_2d;  
    dir_2d = strsplit(dir_2d,filesep); 

    nam = dir_2d{end};    
    nam = append_nam(nam,job);

    dir_2d = fullfile('/',dir_2d{1:end - 1},nam);
    
    if exist(dir_2d,'dir'), rmdir(dir_2d,'s'); end; mkdir(dir_2d);     
else
    dir_2d = '';
end
%==========================================================================

%==========================================================================
function nam = append_nam(nam,job)
if ~isempty(job.channel)
    nam1 = '';
    for i=1:numel(job.channel)
        nam1 = [nam1 job.channel{i}];
    end
    nam = [nam '_' nam1];
end
if job.preproc.reset_origin
    nam = [nam '-ro'];
end        
if job.preproc.do_realign2mni
    nam = [nam '-ra'];
end        
if job.preproc.do_crop
    nam = [nam '-cr'];
    if job.preproc.do_rem_neck
        nam = [nam '-rn'];
    end    
end       
if job.preproc.do_coreg
    nam = [nam '-reg'];
end    
if job.preproc.do_superres
    if job.preproc.mc_superres
        nam = [nam '-mcsr'];
    else
        nam = [nam '-scsr'];
    end                           
end         
if job.preproc.do_denoise
    if job.preproc.mc_denoise
        nam = [nam '-mcden'];
    else
        nam = [nam '-scden'];
    end
end    
if job.preproc.do_ds_inplane
    nam = [nam '-ds'];
end   
if job.preproc.do_ds_throughplane
    nam = [nam '-dsz'];
end        
if job.preproc.resize.do
    nam = [nam '-res'];
end        
if job.preproc.do_bf_correct
    nam = [nam '-bf'];
end          
if job.preproc.do_skull_strip
    nam = [nam '-ss'];
end  
if job.preproc.do_normalise_intensities
    nam = [nam '-ni'];
end   
if job.segment.do
    nam = [nam '-seg'];
end
% deg = ['-deg' num2str(job.preproc.deg_ro) num2str(job.preproc.deg_res) num2str(job.preproc.deg_vx)];
% nam = [nam deg];
%==========================================================================   