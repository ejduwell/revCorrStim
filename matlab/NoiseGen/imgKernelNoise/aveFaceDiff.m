
%% Parameters

imageDir1 = '/home/eduwell/Matlab_Sandbox/imgKernelNoise/baseImgs/londonFaces/neutral_front';
imageDir2 = '/home/eduwell/Matlab_Sandbox/imgKernelNoise/baseImgs/londonFaces/smiling_front';
imTag='*.jpg';
dimsOut=[200,200];

%% Compute Average Image For Each Image Dir

[avgImage1, indAlndImgs1] = averageFace_v3(imageDir1,dimsOut,imTag);
[avgImage2, indAlndImgs2] = averageFace_v3(imageDir2,dimsOut,imTag);

%% Compute Difference Image

dImg=(avgImage1-avgImage2);

%% Plot Ave Images and Diff

figure
subplot(3,1,1);
imagesc(avgImage1);
colormap("gray");
axis equal
axis tight
colorbar
title('Average Img1')

subplot(3,1,2);
imagesc(avgImage2);
colormap("gray");
axis equal
axis tight
colorbar
title('Average Img2')

subplot(3,1,3);
colormap("gray");
imagesc(dImg);
colormap("gray");
axis equal
axis tight
colorbar
title('Difference Img')

%% Compute/Plot an Individual Im1, Im2 and Diff

% grab random index
rng("shuffle");
numImages=size(indAlndImgs1,3);
randIdx=randi(numImages);

% pull images
ind1=indAlndImgs1(:,:,randIdx);
ind2=indAlndImgs2(:,:,randIdx);

% compute differnce image
indDifImg=(ind1-ind2);

% plot
figure
subplot(3,1,1);
imagesc(ind1);
colormap("gray");
axis equal
axis tight
colorbar
title('Individual Img1')

subplot(3,1,2);
imagesc(ind2);
colormap("gray");
axis equal
axis tight
colorbar
title('Individual Img2')

subplot(3,1,3);
colormap("gray");
imagesc(indDifImg);
colormap("gray");
axis equal
axis tight
colorbar
title('Difference Img')
