function Nii = downsample_inplane(Nii,vx1)
mat0 = Nii.mat;             

% Get down-sampling factor
vx0       = spm_misc('vxsize',mat0);   
d         = ((vx0 < vx1).*vx0)./vx1;
d(d == 0) = 1;   

if sum(d)==3
    % Do not downsample if in-plane res Xhat equals in-plane res Y
    warning('do_dsinp::false')
    return
end

% Smooth and resample in-plane                
D      = diag([d, 1]);          
mat_ds = mat0/D;
vx_ds  = spm_misc('vxsize',mat_ds);

X   = Nii.dat(:,:,:);     
dm0 = size(X);    

% fwhm = max(vx_ds./vx0 - 1,0.01);        
% spm_imbasics('smooth_img_in_mem',X,fwhm);                                                 

% Resample using 1st order b-splines             
C          = spm_bsplinc(X,[1 1 1 0 0 0]);            
[x1,y1,z1] = get_sampling_grid(D,dm0);                  
X          = spm_bsplins(C,x1,y1,z1,[1 1 1 0 0 0]);
X(~isfinite(X)) = 0;

% Make sure that the downsampled image has the correct scaling
scl = prod(vx_ds./vx0); 
X   = scl*X;                

fname         = Nii.dat.fname;
[pth,nam,ext] = fileparts(fname);
nfname        = fullfile(pth,['dsip_' nam ext]);

spm_misc('create_nii',nfname,X,mat_ds,Nii.dat.dtype,Nii.descrip);
delete(fname);
Nii = nifti(nfname);
%==========================================================================

%==========================================================================
function [x1,y1,z1] = get_sampling_grid(M,dm)
T          = eye(4)/M;   
dm         = floor(M(1:3,1:3)*dm')';
[x0,y0,z0] = ndgrid(1:dm(1),...
                    1:dm(2),...
                    1:dm(3));

x1 = T(1,1)*x0 + T(1,2)*y0 + T(1,3)*z0 + T(1,4);
y1 = T(2,1)*x0 + T(2,2)*y0 + T(2,3)*z0 + T(2,4);
z1 = T(3,1)*x0 + T(3,2)*y0 + T(3,3)*z0 + T(3,4);  
%==========================================================================