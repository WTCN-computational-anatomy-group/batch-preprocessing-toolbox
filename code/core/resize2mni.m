function Nii = resize2mni(Nii,deg,vx,keep_neck,V_all)

fname = Nii.dat.fname;
V     = spm_vol(fname);
M     = V.mat;
vx1   = spm_misc('vxsize',M);
    
if ~isempty(vx)
    
    if numel(vx) == 1, vx = vx*ones(1,3); end
    
    
    vx1 = vx1./vx;
    
    [img,mat] = resample_img(Nii,[],vx1,deg);
    
    spm_misc('create_nii',fname,img,mat,Nii.dat.dtype,Nii.descrip,Nii.dat.offset,Nii.dat.scl_slope,Nii.dat.scl_inter);    
end

% Get SPM TPM bounding box
pth_spm_tpm = fullfile(spm('dir'),'tpm','TPM.nii,');
V_ref       = spm_vol(pth_spm_tpm);
% bb          = vox_bb(V_ref(1).mat,V_ref(1).dim);

bb = world_bb(V_ref(1).mat,V_ref(1).dim);
bb = [1 1 1; 1./vx1.*diff(bb)];


if keep_neck && ~isempty(V_all)
    % Make sure neck is not cropped    
    bb = adjust_bb_neck(bb,V_all,vx(1));
end

spm_impreproc('subvol',V,bb,'r_',0,true);
%==========================================================================

%==========================================================================
function [img,mat,dm] = resample_img(Nii,Mask,samp,deg,bc)
% Resample an image using deg interpolation, with bc boundary conditions.
% If samp < 1, does down-sampling; if samp > 1, does up-sampling.
%_______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 3, samp = 1; end
if nargin < 4, deg  = 0; end
if nargin < 5, bc   = 0; end

if numel(samp) == 1, samp = samp*ones([1 3]); end
if numel(deg)  == 1, deg  = deg*ones([1 3]);  end
if numel(bc)   == 1, bc   = bc*ones([1 3]);   end

if numel(Nii.dat.dim) ~= 3
    error('numel(Nii.dat.dim) ~= 3')
end

% Input image properties
img  = Nii.dat(:,:,:);    
mat0 = Nii.mat;
dm0  = size(img);

% Mask
if nargin < 2 || isempty(Mask), Mask = true(dm0); end

img(~Mask) = 0;

if samp(1) == 1 || samp(1) == 0
    mat = mat0;
    dm  = dm0;
    
    return
end

% Output image properties
D    = diag([samp 1]);
mat  = mat0/D;
dm   = floor(D(1:3,1:3)*dm0')';

% Make interpolation grid
[x0,y0,z0] = ndgrid(1:dm(1),1:dm(2),1:dm(3));

T = mat0\mat;    

x1 = T(1,1)*x0 + T(1,2)*y0 + T(1,3)*z0 + T(1,4);
y1 = T(2,1)*x0 + T(2,2)*y0 + T(2,3)*z0 + T(2,4);
z1 = T(3,1)*x0 + T(3,2)*y0 + T(3,3)*z0 + T(3,4);

% Resample
img                         = spm_bsplins(img,x1,y1,z1,[deg bc]);    
img(~isfinite(img) | img<0) = 0;
%==========================================================================

%==========================================================================
function bb = vox_bb(mat,dm)
%  world-bb -- get bounding box in world (mm) coordinates

vx = spm_misc('vxsize',mat);
d  = dm(1:3);
% corners in voxel-space
c = [ 1    1    1    1
    1    1    d(3) 1
    1    d(2) 1    1
    1    d(2) d(3) 1
    d(1) 1    1    1
    d(1) 1    d(3) 1
    d(1) d(2) 1    1
    d(1) d(2) d(3) 1 ]';
% % corners in world-space
tc = c;

% bounding box (world) min and max
mn = min(tc,[],2)';
mx = max(tc,[],2)';
bb = [mn; mx];
bb = bb(:,1:3);
bb(2,:) = vx.*bb(2,:);
%==========================================================================

%==========================================================================
function bb = adjust_bb_neck(bb,V,vx)
[mat,dm] = max_bb_orient(V,vx);
bb_max   = world_bb(mat,dm);
bb(1,3)  = bb_max(1,3);    
%==========================================================================

%==========================================================================
% function Nii = resize2mni(Nii,deg,vx,keep_neck,V)
% % FORMAT Nii = resize2mni(Nii,deg,vx,keep_neck,V)
% %
% % ...
% %
% %__________________________________________________________________________
% % Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
% 
% if nargin < 2, deg       = 1; end
% if nargin < 3, vx        = 1; end
% if nargin < 4, keep_neck = true; end
% if nargin < 5, V         = []; end
%             
% if numel(vx) == 1
%     vx = vx*ones(1,3);
% end
% 
% % Read input
% pth_img     = Nii.dat.fname;
% pth_spm_tpm = fullfile(spm('dir'),'tpm','TPM.nii,');
% 
% % Get SPM TPM bounding box
% V_ref = spm_vol(pth_spm_tpm);
% bb    = world_bb(V_ref(1).mat,V_ref(1).dim);
% if keep_neck && ~isempty(V)
%     % Make sure neck is not cropped    
%     bb = adjust_bb_neck(bb,V,vx);
% end
% 
% % Do resizing
% resize_img(pth_img,vx,bb,deg);
%  
% % Update NIfTI
% Nii = nifti(pth_img);
% %==========================================================================
% 
% %==========================================================================
% function resize_img(imnames, Voxdim, BB, deg, ismask)
% %  resize_img -- resample images to have specified voxel dims and BBox
% % resize_img(imnames, voxdim, bb, ismask)
% %
% 
% % Check spm version:
% if exist('spm_select','file') % should be true for spm5
%     spm5 = 1;
% elseif exist('spm_get','file') % should be true for spm2
%     spm5 = 0;
% else
%     error('Can''t find spm_get or spm_select; please add SPM to path')
% end
% 
% spm_defaults;
% 
% % prompt for missing arguments
% if ( ~exist('imnames','var') || isempty(char(imnames)) )
%     if spm5
%         imnames = spm_select(inf, 'image', 'Choose images to resize');
%     else
%         imnames = spm_get(inf, 'img', 'Choose images to resize');
%     end
% end
% % check if inter fig already open, don't close later if so...
% Fint = spm_figure('FindWin', 'Interactive'); Fnew = [];
% if ( ~exist('Voxdim', 'var') || isempty(Voxdim) )
%     Fnew = spm_figure('GetWin', 'Interactive');
%     Voxdim = spm_input('Vox Dims (NaN for "as input")? ',...
%         '+1', 'e', '[nan nan nan]', 3);
% end
% if ( ~exist('BB', 'var') || isempty(BB) )
%     Fnew = spm_figure('GetWin', 'Interactive');
%     BB = spm_input('Bound Box (NaN => original)? ',...
%         '+1', 'e', '[nan nan nan; nan nan nan]', [2 3]);
% end
% if ~exist('ismask', 'var')
%     ismask = false;
% end
% if isempty(ismask)
%     ismask = false;
% end
% 
% % reslice images one-by-one
% vols = spm_vol(imnames);
% for V=vols'
%     % (copy to allow defaulting of NaNs differently for each volume)
%     voxdim = Voxdim;
%     bb = BB;
%     % default voxdim to current volume's voxdim, (from mat parameters)
%     if any(isnan(voxdim))
%         vprm = spm_imatrix(V.mat);
%         vvoxdim = vprm(7:9);
%         voxdim(isnan(voxdim)) = vvoxdim(isnan(voxdim));
%     end
%     voxdim = voxdim(:)';
% 
%     mn = bb(1,:);
%     mx = bb(2,:);
%     % default BB to current volume's
%     if any(isnan(bb(:)))
%         vbb = world_bb(V);
%         vmn = vbb(1,:);
%         vmx = vbb(2,:);
%         mn(isnan(mn)) = vmn(isnan(mn));
%         mx(isnan(mx)) = vmx(isnan(mx));
%     end
% 
%     if sum(bb(:,3)) == 0
%         offset = 0;
%         mn(2)  =  mn(2) + offset;
%         mx(2)  =  mx(2) + offset;
%     end
%     
%     % voxel [1 1 1] of output should map to BB mn
%     % (the combination of matrices below first maps [1 1 1] to [0 0 0])
%     mat = spm_matrix([mn 0 0 0 voxdim])*spm_matrix([-1 -1 -1]);
%     % voxel-coords of BB mx gives number of voxels required
%     % (round up if more than a tenth of a voxel over)
%     imgdim = ceil(mat \ [mx 1]' - 0.1)';
% 
%     % output image
%     VO            = V;
%     [pth,nam,ext] = fileparts(V.fname);
%     VO.fname      = fullfile(pth,['r_' nam ext]);
%     VO.dim(1:3)   = imgdim(1:3);
%     VO.mat        = mat;
%     VO = spm_create_vol(VO);
%     spm_progress_bar('Init',imgdim(3),'reslicing...','planes completed');
%     for i = 1:imgdim(3)
%         M = inv(spm_matrix([0 0 -i])*inv(VO.mat)*V.mat);
%         img = spm_slice_vol(V, M, imgdim(1:2), deg);
%         if ismask
%             img = round(img);
%         end
%         spm_write_plane(VO, img, i);
%         spm_progress_bar('Set', i)
%     end
%     spm_progress_bar('Clear');
% end
% % call spm_close_vol if spm2
% if ~spm5
%     spm_close_vol(VO);
% end
% if (isempty(Fint) && ~isempty(Fnew))
%     % interactive figure was opened by this script, so close it again.
%     close(Fnew);
% end
% disp('Done.')
% %==========================================================================
% 
%==========================================================================
function bb = world_bb(mat,dm)
%  world-bb -- get bounding box in world (mm) coordinates

d = dm(1:3);
% corners in voxel-space
c = [ 1    1    1    1
    1    1    d(3) 1
    1    d(2) 1    1
    1    d(2) d(3) 1
    d(1) 1    1    1
    d(1) 1    d(3) 1
    d(1) d(2) 1    1
    d(1) d(2) d(3) 1 ]';
% corners in world-space
tc = mat(1:3,1:4)*c;

% bounding box (world) min and max
mn = min(tc,[],2)';
mx = max(tc,[],2)';
bb = [mn; mx];
%==========================================================================
% 
% %==========================================================================
% function bb = adjust_bb_neck(bb,V,vx)
% [mat,dm] = max_bb_orient(V,vx);
% bb_max   = world_bb(mat,dm);
% bb(1,3)  = bb_max(1,3);    
% %==========================================================================
% 
%==========================================================================
function [mat,dm] = max_bb_orient(V,vx)
% Calculate orientation matrix and dimensions from maximum bounding-box
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

mn = [ Inf  Inf  Inf]';
mx = [-Inf -Inf -Inf]';
for i=1:numel(V)
        d = V(i).dim;
        if numel(d)==2, d(3) = 1; end

        t = uint8(0:7);
        c = diag(d+1)*double([bitshift(bitand(t,bitshift(uint8(1),1-1)),1-1)
                              bitshift(bitand(t,bitshift(uint8(1),2-1)),1-2)
                              bitshift(bitand(t,bitshift(uint8(1),3-1)),1-3)]);
        c = bsxfun(@plus,V(i).mat(1:3,1:3)*c,V(i).mat(1:3,4));
        mx = max(mx,max(c,[],2));
        mn = min(mn,min(c,[],2));
end
mat = spm_matrix(mn-1)*diag([vx 1])*spm_matrix(-[1 1 1]);
dm  = ceil((mat\[mx'+1 1]')');
dm  = dm(1:3);
%==========================================================================