function obj = init_obj_preproc(pars)
M   = numel(pars.dat);
obj = cell(1,M);
cnt = 1;
for m=1:M           
    scans = pars.dat{m}.scans;
    N     = numel(scans{1});
    
    if scans{1}{1}{1}.dim(3)==1
        pars.dat{m}.preproc.write_2d = false;
    end
    
    %----------------------------------------------------------------------    
    if isempty(pars.dat{m}.dir_preproc)
        dir_preproc          = fileparts(scans{1}{1}{1}.fname);
        dir_preproc          = strsplit(dir_preproc,filesep);
        dir_preproc{end - 3} = [dir_preproc{end - 3} '_preproc'];
        if pars.dat{m}.preproc.do_realign2mni
            dir_preproc{end - 3} = [dir_preproc{end - 3} '-ra'];
        end        
        if pars.dat{m}.preproc.do_crop
            dir_preproc{end - 3} = [dir_preproc{end - 3} '-cr'];
            if pars.dat{m}.preproc.do_rem_neck
                dir_preproc{end - 3} = [dir_preproc{end - 3} '-rn'];
            end    
        end       
        if pars.dat{m}.preproc.do_coreg && N>1 
            dir_preproc{end - 3} = [dir_preproc{end - 3} '-reg'];
        end    
        if pars.dat{m}.preproc.do_reslice && N>1 
            dir_preproc{end - 3} = [dir_preproc{end - 3} '-res'];        
        end                  
        if pars.dat{m}.preproc.vx
            dir_preproc{end - 3} = [dir_preproc{end - 3} '-vx'];        
        end                      
        if pars.dat{m}.preproc.do_superres
            dir_preproc{end - 3} = [dir_preproc{end - 3} '-sr'];
        end            
        if pars.dat{m}.preproc.do_denoise
            dir_preproc{end - 3} = [dir_preproc{end - 3} '-den'];
        end    
        if pars.dat{m}.preproc.do_bf_correct
            dir_preproc{end - 3} = [dir_preproc{end - 3} '-bf'];
        end          
        if pars.dat{m}.preproc.do_skull_strip
            dir_preproc{end - 3} = [dir_preproc{end - 3} '-ss'];
        end  
        if pars.dat{m}.preproc.normalise_intensities
            dir_preproc{end - 3} = [dir_preproc{end - 3} '-ni'];
        end   
        dir_preproc_2d = fullfile('/',dir_preproc{2:end - 4},['2d_' dir_preproc{end - 3}]);
        dir_preproc    = fullfile('/',dir_preproc{2:end - 3});    
    else                
        dir_preproc    = pars.dat{m}.dir_preproc;
        dir_preproc1   = strsplit(dir_preproc,filesep);
        dir_preproc_2d = fullfile('/',dir_preproc1{2:end - 1},['2d_' dir_preproc1{end}]);
    end
    
    if exist(dir_preproc,'dir'), rmdir(dir_preproc,'s'); end; mkdir(dir_preproc);  
    
    if pars.dat{m}.preproc.write_2d
        if exist(dir_preproc_2d,'dir'), rmdir(dir_preproc_2d,'s'); end; mkdir(dir_preproc_2d); 
    end
     
    S      = numel(scans);
    obj{m} = cell(1,S);
    for s=1:S                
        obj1 = struct;
                         
        obj1.scans          = scans{s};
        obj1.labels         = pars.dat{m}.labels{s};      
        obj1.modality       = pars.dat{m}.modality;
        obj1.s              = cnt;                                  
        obj1.dir_preproc    = dir_preproc;
        obj1.dir_preproc_2d = dir_preproc_2d;
        
        obj1.preproc = pars.dat{m}.preproc;

        obj{m}{s} = obj1;        
        cnt       = cnt + 1;
    end
end
%==========================================================================