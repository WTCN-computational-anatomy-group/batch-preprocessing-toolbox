function tau = estimate_ct_tau(x,y,vx)

msk      = isfinite(x) & x>-1010  & x<-990; % Mask all but air
x       = x + 1000; % Adjust mean
x(~msk) = NaN;     

if isempty(y)
    y = zeros(size(x),'single');
else
    y = y + 1000;
end
y(~msk) = NaN;

res = x - y;

% Compute gradients
d         = size(res);
if numel(d)==2, d(3) = 1; end
[id{1:3}] = ndgrid(single(1:d(1)),single(1:d(2)),single(1:d(3)));
id        = cat(4,id{:});
fc        = spm_diffeo('bsplinc',res,[1 1 1 1 1 1]);
[~,Dx,Dy,Dz] = spm_diffeo('bsplins',fc,id,[1 1 1 1 1 1]);
clear id

res(~isfinite(res)) = [];
Dx(~isfinite(Dx) | Dx==0) = [];
Dy(~isfinite(Dy) | Dy==0) = [];
Dz(~isfinite(Dz) | Dz==0) = [];

% Estimate FWHM
fwhm = sqrt(4*log(2))*sum(abs(res(:)))./[sum(abs(Dx(:))) sum(abs(Dy(:))) sum(abs(Dz(:)))];
fwhm(~isfinite(fwhm)) = 0;

ff  = compute_ff_Gauss(fwhm,vx);
tau = 1/ff;

