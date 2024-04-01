
%% Parameters

imgInPath="/home/eduwell/Matlab_Sandbox/imgKernelNoise/baseImgs/londonFaces/smiling_front/103_08.jpg";
%imgInPath="/home/eduwell/Matlab_Sandbox/imgKernelNoise/face/031_03.jpg";
%imgInPath="/home/eduwell/Matlab_Sandbox/imgKernelNoise/face/107_08.jpg";
%imgInPath="/home/eduwell/Matlab_Sandbox/imgKernelNoise/baseImgs/BaseImgSF_10L.png";

%imgInPath='/home/eduwell/Matlab_Sandbox/imgKernelNoise/baseImgs/revCorrStim/genaTexBI_occ_0_ori_m__CR_0_CRpw_1_CRpl_0_CRal_1_CRob_2_CRdm_3214_T_0spheres_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_0_lum2_182.1429_LWT_0.5_Lal_0_5_Ocon_2.png';

%Specify desired size?
selectSize=1; % if 1, this means that we will use/resize the input image the selected size below
desiredSize=[256,256]; % (MUST BE SQUARE IMAGE AND DIMZ MUST BE POWER OF 2)

% specifiy number of samples of each kernel subsample image you want to be
% used per output noise frame..
%Note: length must be = to the image dimension power of two minus one (ie if 512x512-->9, if 256x256-->8)
imgsPerKrnl=[1,1,1,1,1,1,1,1,1]; % adjusts the total number of subtile images from each respective kernel incorporated into the noise
krnlImgWgts=[1,1,1,1,1,1,1,1,1]; % adjusts relative weights applied to subtile images from each respective kernel


totKimz=sum(imgsPerKrnl); % count up total number of intermediate images..

smoothTiles=[0,0,0,0,0,0,0,0,0]; % apply smoothing kernel to intermediate subtile images?
gSmthK1sz=[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]; % size of gaussian smoothing kernel for subtiles..
smoothFinal=0; % apply smoothing kernel to final image?
gSmthK2sz=0.5; % size of gaussian smoothing kernel for final image..

BI_WT=0.35;
N_WT=1-BI_WT;
%% Subsample the image with kernelSubSampler

[adjImg, kSampleArrayz, kernelDim_ImFrac, kernelSqrDims, kernelPosCnts] = kernelSubSampler_v2(imgInPath,selectSize,desiredSize);

%% Generate intermediate "Tile Images" by tiling the image space with kernel subsamples
rng("shuffle"); % Make sure the random number generator is shuffled..

% get the number of kernels.. 
nKernels=length(kernelSqrDims);
kNumberz=linspace(1,nKernels,nKernels); % set kNumberz equal to list 1-->nKernels..

kernelFrames=zeros(size(adjImg,1),size(adjImg,2),totKimz); % preallocate stack for intermediate kernel frames..
countr=1; % initialize countr..

% Loop through the kernels and generate kernel subsample frames for each..
for ii = 1:nKernels

kNumber=kNumberz(ii);

for nn = 1:imgsPerKrnl(ii)
% pull info for this kernel..
kDms=[kernelSqrDims(kNumber),kernelSqrDims(kNumber)];
kFrkn=kernelDim_ImFrac(kNumber);

tileImg=uint8(zeros(size(adjImg,1),size(adjImg,2))); % preallocate..
imsz=size(adjImg,1);

for yy = 1:kFrkn
    yStrt=1+((imsz/kFrkn)*(yy-1));
    for xx = 1:kFrkn

        xStrt=1+((imsz/kFrkn)*(xx-1));

        % select a random kernel window...
        kIndx=randi(kernelPosCnts(kNumber)); % grab a random kernel position index..
        
        % get coordinates for kernel at that index
        kCoords=kSampleArrayz{1,kNumber}(kIndx,:);

        % use coordinates to grab sub-image sample from base image
        kSubImg=uint8(adjImg(kCoords(2):kCoords(4),kCoords(1):kCoords(3)));
        %kSubImg=kSampleArrayz{1,kNumber}(:,:,kIndx); % grab rand
        
        % select a random flip condition..
        flipCon1=randi(2); % horz flip
        flipCon2=randi(2); % vert flip
        flipCon3=randi(2); % invert

        if flipCon1==2
        % Flip the image horizontally
        kSubImg = fliplr(kSubImg);
        end
        
        if flipCon2==2
        % Flip the image vertically
        kSubImg = flipud(kSubImg);
        end

        if flipCon3==2
        % Flip the image both horizontally and vertically
        kSubImg = kSubImg';
        end

        tileImg(yStrt:(yStrt+kDms(1)-1),xStrt:(xStrt+kDms(1)-1))=kSubImg;
    end

end

% Convert to double and scale from 0-1
tileImg = double(tileImg);
maxVal=max(max(tileImg));
tileImg = tileImg./maxVal; % divide all by max value to rescale range from 0-1

% smooth tile image if requested..
if smoothTiles(ii)==1
tileImg = imgaussfilt(tileImg,gSmthK1sz(ii));
end

% apply requested weighting
tileImg = tileImg.*krnlImgWgts(ii);

% Apply random shifts to jitter the position of the tile boundaries 
vShftOptions={'up','down'};
hShftOptions={'left','right'};
shftCon1=randi(2); % vert shft random idx selector
shftCon2=randi(2); % horz shft random idx selector
vShftDir=vShftOptions{shftCon1};
hShftDir=hShftOptions{shftCon2};
vShftVal=randi(imsz);
hShftVal=randi(imsz);
tileImg = shiftImage(tileImg, vShftDir, vShftVal); % apply vertical shift
tileImg = shiftImage(tileImg, hShftDir, hShftVal); % apply horizontal shift

% save tile image frame into kernelFrames stack
kernelFrames(:,:,countr)=tileImg;
countr=countr+1;

end

end
%==========================================================================

%% Combine the stack of kernel subsample images linearly to form noise image

nSubImgs=size(kernelFrames,3); % get total number of intermediate kernel images..
noiseImg=sum(kernelFrames,3)./nSubImgs; % take weighted average of the stack..
% smooth final noise image if requested..
if smoothFinal==1
noiseImg = imgaussfilt(noiseImg,gSmthK2sz);
end

noiseImg=uint8(rescale(noiseImg,0,255)); % scale values of final image to range from 0 to 255 and convert to uint8..

%% Generate other noise type images for comparison
whiteNoiseImg=uint8(randi(255,[size(adjImg,1),size(adjImg,1)]));

%% Attempt to reconstruct the original base image using xCorKrnlNzReArrange for each noise type

% The general hypothesis being tested here:
%--------------------------------------------------------------------------
% Noise types who's "feature spectra" more closely match that of the the
% base image should contain more of the relevant "building blocks" which
% make up the base image. You should, therefore, be able to more effectively
% "reconstruct" the base image using features (kernel windows) from noise
% image frames which are more rich in the base image features.
%--------------------------------------------------------------------------
% How to test?:
%--------------------------------------------------------------------------
% Much like the approach used to generate the tailored kernel noise: 
%
%   - Use a series of kernel window sizes to chop the noise up into arrays of
%   smaller and smaller square chunks/kernel windows. 
%       o Only this time, the kernel will not be allowed to slide 
%         continuously across the x/y dimensions of the images, but will
%         instead sample a regularly spaced grid such that each pixel only 
%         gets sampled/windowed once per kernel. 
%
%   - For each kernel window in each location for each kernel: we will 
%     cross-correlate that window with all possible locations on the base 
%     image and place it in the place it where it correlates best/highest.
%
%   - The contents of the kernel window will then be weighted (multiplied)
%     by the max correlation value it achieves there and added to that
%     location in the output image.
%
%   - A "tally image" will also be created for each individual kernel and
%     also combined across kernels which keep track (or 'tally') pixel-wise 
%     of how many kernels fell in each location for their "max corrlation"
%     location in the cross-correlation operation described above.
%
%   - These tally images are then used to compute the weighted average
%     value across kernel windows for each pixel.
%
%   - This final "tally weighted" image for each kernel should reflect the 
%     best reconstruction of the base image achievable by rearranging using 
%     'that-sized' kernel chunks/windows of the noise.
%
%   - If hypothesis above holds up: noise images more rich in the base 
%     image's features should produce output images which more closesly
%     resemble the base image than those "less rich" in the base image's
%     features.
%   
%   ...
%
%   - Follow up hypothesis 1: The extent to which images produced with 
%     different sized kernels resemble the base image should reflect the 
%     extent to which features of rougly that size (that kernel window size) 
%     present in the base image are present in and reconstructable from the 
%     noise..
%
%   - Follow up hypothesis 2: The extent to which images produced with 
%     different sized kernels resemble the base image should directly
%     translate to/predict the number of trials that will be required using
%     that noise to generate a usable CI. Greater correspondance of
%     reconstructed BI from a noise frame should equate to less trials
%     required for producing a CI of any given/specified quality.
%--------------------------------------------------------------------------

MinMaxCor=1; % specify whether you want to rearrange kernel windows to fall at locations min or max correlation value
nKernReArr=6; % set number of kernels you want to re-arrange..

% Run xCorKrnlNzReArrange_v2 on Tailored Kernel Noise
tic
[xCorAdjNzFrames, KrnlXCorData, combTallyImg, indTallyImgz] = xCorKrnlNzReArrange_v2(adjImg, noiseImg, kernelDim_ImFrac, kernelSqrDims,MinMaxCor,nKernReArr);
toc

% Run xCorKrnlNzReArrange_v2 on White Noise
tic
[xCorAdjNzFrames2, KrnlXCorData2, combTallyImg2, indTallyImgz2] = xCorKrnlNzReArrange_v2(adjImg, whiteNoiseImg, kernelDim_ImFrac, kernelSqrDims,MinMaxCor,nKernReArr);
toc

% Generate "tally weighted" images for all kernels combined
% for both tailored noise and white noise
% Generate ssim images/values between these and the base image too..
xCorOptNzImg = weightByTallyImg(combTallyImg,xCorAdjNzFrames);
xCorOptNzImg2 = weightByTallyImg(combTallyImg2,xCorAdjNzFrames2);

[xCorOptNzImg_ssimVal, xCorOptNzImg_ssimImg] = ssim(xCorOptNzImg,adjImg);
[xCorOptNzImg2_ssimVal, xCorOptNzImg2_ssimImg] = ssim(xCorOptNzImg2,adjImg);

% Now generate "tally weighted" images for each individual kernel
% for both tailored noise and white noise
%
% Also: compute the ssim between each of the resulting images and the 
% base image..
%--------------------------------------------------------------------------
% TAILORED NOISE:
indTallyImgReconz1=uint8(zeros(size(indTallyImgz,1),size(indTallyImgz,2),size(indTallyImgz,3))); % preallocate
TlyImgRcn1_ssimImz=zeros(size(indTallyImgz,1),size(indTallyImgz,2),size(indTallyImgz,3)); % preallocate
TlyImgRcn1_ssimValz=zeros(size(indTallyImgz,3,1)); % preallocate
for qq=1:size(indTallyImgz,3)
    tlyImTmp=indTallyImgz(:,:,qq); % grab this kernel's individual tally img
    xCorOptImgTmp=xCorAdjNzFrames(:,:,qq); % grab this kernel's sum img
    imTmp = weightByTallyImg(tlyImTmp,xCorOptImgTmp); % make weighted ave for this kernel
    indTallyImgReconz1(:,:,qq)=imTmp;   % save in output array..
    % Compute SSIM between reconstructed img and base img..
    [imTmp_ssimVal, imTmp_ssimImg]=ssim(imTmp,adjImg);
    TlyImgRcn1_ssimImz(:,:,qq)=imTmp_ssimImg;
    TlyImgRcn1_ssimValz(qq,1)=imTmp_ssimVal;
end
% WHITE NOISE:
indTallyImgReconz2=uint8(zeros(size(indTallyImgz2,1),size(indTallyImgz2,2),size(indTallyImgz2,3))); % preallocate
TlyImgRcn2_ssimImz=zeros(size(indTallyImgz2,1),size(indTallyImgz2,2),size(indTallyImgz2,3)); % preallocate
TlyImgRcn2_ssimValz=zeros(size(indTallyImgz2,3,1)); % preallocate
for qq=1:size(indTallyImgz2,3)
    tlyImTmp=indTallyImgz2(:,:,qq); % grab this kernel's individual tally img
    xCorOptImgTmp=xCorAdjNzFrames2(:,:,qq); % grab this kernel's sum img
    imTmp = weightByTallyImg(tlyImTmp,xCorOptImgTmp); % make weighted ave for this kernel
    indTallyImgReconz2(:,:,qq)=imTmp;   % save in output array..
    % Compute SSIM between reconstructed img and base img..
    [imTmp_ssimVal, imTmp_ssimImg]=ssim(imTmp,adjImg);
    TlyImgRcn2_ssimImz(:,:,qq)=imTmp_ssimImg;
    TlyImgRcn2_ssimValz(qq,1)=imTmp_ssimVal;
end
%--------------------------------------------------------------------------

%% Make BI/Noise Combo Images

noiseBICombo=((BI_WT.*adjImg+N_WT.*noiseImg)./2);
noiseBICombo=uint8(rescale(noiseBICombo,0,255));
whtnoiseBICombo=((BI_WT.*adjImg+N_WT.*whiteNoiseImg)./2);
whtnoiseBICombo=uint8(rescale(whtnoiseBICombo,0,255));
xCorOptNzImgBICombo=((BI_WT.*adjImg+N_WT.*xCorOptNzImg)./2);
xCorOptNzImgBICombo=uint8(rescale(xCorOptNzImgBICombo,0,255));

%% Plot Images

% PLOT (1)
% Base Image with White Noise vs. Tailored Kernel Noise Plot
%##########################################################################
figure
subplot(2,6,[1,2]);
imshow(adjImg);
title("Base Image")

subplot(2,6,[3,4]);
imshow(whiteNoiseImg);
title("White Noise")

subplot(2,6,[5,6]);
imshow(noiseImg);
title("Tailored 'Kernel' Noise")

subplot(2,6,[8,9]);
imshow(whtnoiseBICombo);
title("Base/White Noise Weighted Average")

subplot(2,6,[10,11]);
imshow(noiseBICombo);
title("Base/Tailored Noise Weighted Average")
%##########################################################################

% PLOT (2)
% Reconstruction of Base Image with Tailored vs. White Noise using all
% kernel windows combined..
%##########################################################################
figure
subplot(2,6,[1,2]);
imshow(adjImg);
title("Base Image")

subplot(2,6,[3,4]);
imshow(whiteNoiseImg);
title("White Noise")

subplot(2,6,[5,6]);
imshow(noiseImg);
title("Tailored 'Kernel' Noise")


subplot(2,6,[8,9]);
imshow(xCorOptNzImg2);
title("BI Reconstructed with White Noise")

subplot(2,6,[10,11]);
imshow(xCorOptNzImg);
title("BI Reconstructed with Tailored Kernel Noise")
%##########################################################################

% PLOT (3)
% Reconstruction of Base Image with Tailored vs. White Noise using each
% individual kernel window..
%##########################################################################
figure
pChkSz=round(nKernReArr/3);
pBegin=1;
pEnd=((pBegin+pChkSz)-1);
subplot(3,nKernReArr,[pBegin,pEnd]);
imshow(adjImg);
title("Base Image")

pBegin=(pEnd+1);
pEnd=((pBegin+pChkSz)-1);
subplot(3,nKernReArr,[pBegin,pEnd]);
imshow(whiteNoiseImg);
title("White Noise")

pBegin=(pEnd+1);
pEnd=((pBegin+pChkSz)-1);
subplot(3,nKernReArr,[pBegin,pEnd]);
imshow(noiseImg);
title("Tailored 'Kernel' Noise")

pstart=pEnd; % initialize..
for uu = 1:nKernReArr
    pstart=(pstart+1);
    subplot(3,nKernReArr,pstart);
    im2Show=indTallyImgReconz2(:,:,uu);
    imshow(im2Show);
    titleTmp=strcat("White Nz BI Recon With K",num2str(uu));
    title(titleTmp)
end

for uu = 1:nKernReArr
    pstart=(pstart+1);
    subplot(3,nKernReArr,pstart);
    im2Show=indTallyImgReconz1(:,:,uu);
    imshow(im2Show);
    titleTmp=strcat("Tailored Nz BI Recon With K",num2str(uu));
    title(titleTmp)
end

% PLOT (4)
% SSIM of Reconstructed Base Images using Tailored vs. White Noise for each
% individual kernel window size..
%##########################################################################
figure
pChkSz=round(nKernReArr/3);
pBegin=1;
pEnd=((pBegin+pChkSz)-1);
subplot(3,nKernReArr,[pBegin,pEnd]);
imshow(adjImg);
title("Base Image")

pBegin=(pEnd+1);
pEnd=((pBegin+pChkSz)-1);
subplot(3,nKernReArr,[pBegin,pEnd]);
imshow(xCorOptNzImg2_ssimImg);
title("White Noise SSIM to BI: All Kernels Combined")

pBegin=(pEnd+1);
pEnd=((pBegin+pChkSz)-1);
subplot(3,nKernReArr,[pBegin,pEnd]);
imshow(xCorOptNzImg_ssimImg);
title("Tailored Noise SSIM to BI: All Kernels Combined")

pstart=pEnd; % initialize..
for uu = 1:nKernReArr
    pstart=(pstart+1);
    subplot(3,nKernReArr,pstart);
    im2Show=TlyImgRcn2_ssimImz(:,:,uu);
    imshow(im2Show);
    titleTmp=strcat("SSIM White w/ K",num2str(uu)," and BI");
    title(titleTmp)
end

for uu = 1:nKernReArr
    pstart=(pstart+1);
    subplot(3,nKernReArr,pstart);
    im2Show=TlyImgRcn1_ssimImz(:,:,uu);
    imshow(im2Show);
    titleTmp=strcat("SSIM Tailored w/ K",num2str(uu)," and BI");
    title(titleTmp)
end

% PLOT (5)
% Bar Plots of Overall SSIM of Reconstructed Base Images using Tailored 
% vs. White Noise for each individual kernel window size..
%##########################################################################

% Preallocate..
xlabelz = num2cell(zeros(1, nKernReArr+1));

ssimValz2Plot=zeros(2,nKernReArr+1);

for ff = 1:nKernReArr+1

    % for the last one label as "all-Ks-combined"
    % and plot the combined kernel overall ssim values..
    if ff == (nKernReArr+1)
        xlabelz(1,ff)={'all-Ks-combined'};
        ssimValz2Plot(1,ff)=xCorOptNzImg2_ssimVal;
        ssimValz2Plot(2,ff)=xCorOptNzImg_ssimVal;
    else
        lablTmp=strcat('kernel',num2str(ff));
        xlabelz(1,ff)={lablTmp};
        ssimValz2Plot(1,ff)=TlyImgRcn2_ssimValz(ff);
        ssimValz2Plot(2,ff)=TlyImgRcn1_ssimValz(ff);
    end
end

figure

xlabelz = categorical(xlabelz);
bar(xlabelz,ssimValz2Plot)
title("SSIM Between BI and Recon-BIs Created with White vs. Tailored Noise by Kernel Size");
ylabel("Overall Structural Similarity with BI");
xlabel("Kernel Size Used During Recon Proceedure");
%##########################################################################
