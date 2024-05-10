

% Usage BUTTERWORTHBPF(I,DO,D1,N)
% Example

lowerBnd=2; % lower frequency bound
upperBnd=10; % upper frequency bound
filterOrder=8; % filter order
image = CI_results1occ; % image in
pltFig=0;

filtered_image = butterworthbpf_v2(image,lowerBnd,upperBnd,filterOrder,pltFig);


figure;
imagesc(image);
colormap('gray');
axis image
colorbar;
title('original image');

figure;
imagesc(filtered_image);
colormap('gray');
axis image
colorbar;
title('filtered image');