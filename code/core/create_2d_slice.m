%==========================================================================
function create_2d_slice(fname,axis_2d)
V = spm_vol(fname);

spm_impreproc('nm_reorient',fname,vxsize(V.mat),1,'ro_');          
[pth,nam,ext] = fileparts(fname);
nfname        = fullfile(pth,['ro_' nam ext]);
delete(fname);

V  = spm_vol(nfname);
dm = V.dim;

if axis_2d==1
    d1 = floor(dm(1)/2) + 1;
    bb = [d1 d1;-inf inf;-inf inf];
elseif axis_2d==2
    d1 = floor(dm(2)/2) + 1;
    bb = [-inf inf;d1 d1;-inf inf];
elseif axis_2d==3 
    d1 = floor(dm(3)/2) + 1;
    bb = [-inf inf;-inf inf;d1 d1];
end                

spm_impreproc('subvol',V,bb','2d_');      
[pth,nam,ext] = fileparts(nfname);   
nfname1       = fullfile(pth,['2d_' nam ext]);   
delete(nfname);

spm_impreproc('reset_origin',nfname1);
%==========================================================================