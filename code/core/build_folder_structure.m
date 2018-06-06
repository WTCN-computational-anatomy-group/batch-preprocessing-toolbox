function [dir_preproc,dir_2d] = build_folder_structure(job)
dir_preproc = job.dir_preproc;  
dir_preproc = strsplit(dir_preproc,filesep); 

nam = dir_preproc{end};    
nam = append_nam(nam,job);

dir_preproc = fullfile('/',dir_preproc{1:end - 1},nam);

if exist(dir_preproc,'dir'), rmdir(dir_preproc,'s'); end; mkdir(dir_preproc); 

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
if job.preproc.do_ds_inplane
    nam = [nam '-ds'];
end
if job.preproc.do_reslice 
    nam = [nam '-res'];        
end                  
if ~isempty(job.preproc.vx)
    nam = [nam '-vx'];        
end                      
if job.preproc.do_superres
    nam = [nam '-sr'];
end            
if job.preproc.do_denoise
    nam = [nam '-den'];
end    
if job.preproc.do_bf_correct
    nam = [nam '-bf'];
end          
if job.preproc.do_skull_strip
    nam = [nam '-ss'];
end  
if job.preproc.normalise_intensities
    nam = [nam '-ni'];
end   
if job.segment.write_tc
    nam = [nam '-seg'];
end
%==========================================================================   