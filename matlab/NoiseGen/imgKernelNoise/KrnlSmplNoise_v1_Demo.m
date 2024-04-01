
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
imgsPerKrnl=[1,1,1,1,1,1,1,0,0]; % adjusts the total number of subtile images from each respective kernel incorporated into the noise
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

whiteNoiseImg=uint8(randi(255,[size(adjImg,1),size(adjImg,1)]));
noiseImg=uint8(rescale(noiseImg,0,255)); % scale values of final image to range from 0 to 255 and convert to uint8..
noiseBICombo=((BI_WT.*adjImg+N_WT.*noiseImg)./2);
noiseBICombo=uint8(rescale(noiseBICombo,0,255));
whtnoiseBICombo=((BI_WT.*adjImg+N_WT.*whiteNoiseImg)./2);
whtnoiseBICombo=uint8(rescale(whtnoiseBICombo,0,255));


subplot(2,6,[1,2]);
imshow(adjImg);
title("Base Image")
subplot(2,6,[3,4]);
imshow(whiteNoiseImg);
title("White Noise")
subplot(2,6,[5,6]);
imshow(noiseImg);
title("Tailored 'Kernel' Noise")
subplot(2,6,[7,9]);
imshow(whtnoiseBICombo);
title("Base/White Noise Weighted Average")
subplot(2,6,[10,12]);
imshow(noiseBICombo);
title("Base/Tailored Noise Weighted Average")
