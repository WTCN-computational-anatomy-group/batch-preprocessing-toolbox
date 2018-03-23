function [ff,s] = compute_ff_Gauss(fwhm,vx,sk)
% Fudge Factor - to (approximately) account for non-independence of voxels.
% Note that variances add, and that Var[a*x + b*y] = a^2*Var[x] + b^2*Var[y]
% Therefore the variance of i.i.d. noise after Gaussian smoothing is equal
% to the sum of the Gaussian function squared times the original variance.
% A Gaussian is given by g=sqrt(2*pi*s^2)^(-1/2)*exp(-0.5*x.^2/s^2);
% After squaring, this is (2*pi*s^2)^(-1)*exp(-x.^2/s^2), which is a scaled
% Gaussian. Letting s2 = 2/sqrt(2), this is equal to
% (4*pi*s^2)^(-1/2)*(2*pi*s2^2)^(-1/2)*exp(-0.5*x.^2/s2^2), from which
% the (4*pi*s^2)^(-1/2) factor comes from.
if nargin<2, vx = [1 1 1]; end
if nargin<3, sk = [1 1 1]; end

fwhm = fwhm + mean(vx); 
s    = fwhm/sqrt(8*log(2));                 % Standard deviation
ff   = prod(4*pi*(s./vx./sk).^2 + 1)^(1/2); 
%==========================================================================