f1 = '/data/mbrud/populations/preproc/CROMIS-test-ra-cr-rn-vx-den/vx_cr_ro_isosCROMIS2ICH_01024-0004-00003-000001.nii';
f2 = '/data/mbrud/populations/preproc/DELIRIUM-test-ra-cr-rn-vx-den/vx_cr_ro_1000_s02024126-0007-00001-000001.nii';
f3 = '/data/mbrud/populations/preproc/ROB-test-ra-cr-rn-vx-den/vx_cr_ro_P2012.nii';
m  = 'CT';

figure(111)

subplot(131)
[sd,X] = estimate_sd(f1,m);
d3 = floor(size(X,3)/2) + 1;
imagesc(X(:,:,d3)',[0 100]); axis off image xy; colormap(gray);
title(['sd=' num2str(sd)])

subplot(132)
[sd,X] = estimate_sd(f2,m);
d3 = floor(size(X,3)/2) + 1;
imagesc(X(:,:,d3)',[0 100]); axis off image xy; colormap(gray);
title(['sd=' num2str(sd)])

subplot(133)
[sd,X] = estimate_sd(f3,m);
d3 = floor(size(X,3)/2) + 1;
imagesc(X(:,:,d3)',[0 100]); axis off image xy; colormap(gray);
title(['sd=' num2str(sd)])