function SmoSuf = ComputeSmoSuf_Gauss(x,y,vx,modality)
if isempty(y),  y  = zeros(size(x),'single'); end
if isempty(vx), vx = [1 1 1]; end

if strcmp(modality,'CT')
    msk0 = isfinite(x(:)) & x(:)>-1010  & x(:)<-990;
end

res = x - y;

% Compute sum of absolute gradients
SmoSuf    = zeros(1,8);
tmp       = res;
msk       = isfinite(tmp(:)) & msk0;
SmoSuf(1) = sum(msk(:));
SmoSuf(2) = sum(tmp(msk(:)).^2);

[Dx,Dy,Dz] = spm_imbasics('grad',res,vx);

tmp       = Dx;
msk       = isfinite(tmp(:)) & msk0;
SmoSuf(3) = sum(msk(:));
SmoSuf(4) = sum(tmp(msk(:)).^2);
clear Dx tmp msk

tmp       = Dy;
msk       = isfinite(tmp(:)) & msk0;
SmoSuf(5) = sum(msk(:));
SmoSuf(6) = sum(tmp(msk(:)).^2);

clear Dy tmp msk
tmp       = Dz;
if numel(tmp)>1
    msk   = isfinite(tmp(:)) & msk0;
else
    msk   = isfinite(tmp(:));
end
SmoSuf(7) = sum(msk(:));
SmoSuf(8) = sum(tmp(msk(:)).^2);
clear Dz tmp msk