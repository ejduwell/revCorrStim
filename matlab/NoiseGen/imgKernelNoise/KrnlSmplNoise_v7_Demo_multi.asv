
%% Parameters

%imgInPath="/home/eduwell/Matlab_Sandbox/imgKernelNoise/baseImgs/londonFaces/smiling_front/103_08.jpg";
%imgInPath="/home/eduwell/Matlab_Sandbox/imgKernelNoise/baseImgs/londonFaces/neutral_front/103_03.jpg";
%imgInPath="/home/eduwell/Matlab_Sandbox/imgKernelNoise/baseImgs/londonFaces/smiling_front/103_08.jpg";
%imgInPath="/home/eduwell/Matlab_Sandbox/imgKernelNoise/face/031_03.jpg";
%imgInPath="/home/eduwell/Matlab_Sandbox/imgKernelNoise/face/107_08.jpg";


% Path to base image used to generate kernel noise..
imgInPath="/home/eduwell/Matlab_Sandbox/imgKernelNoise/baseImgs/BaseImgSF_6L.png"; 
%imgInPath='/home/eduwell/Matlab_Sandbox/imgKernelNoise/baseImgs/SimDifImgs/difBIs/difBI_6Hz_2_L-R.png';
%imgInPath='/home/eduwell/Matlab_Sandbox/imgKernelNoise/baseImgs/revCorrStim/genaTexBI_occ_0_ori_m__CR_0_CRpw_1_CRpl_0_CRal_1_CRob_2_CRdm_3214_T_0spheres_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_0_lum2_182.1429_LWT_0.5_Lal_0_5_Ocon_2.png';

%Specify desired size?
selectSize=1; % if 1, this means that we will use/resize the input image the selected size below
desiredSize=[256,256]; % (MUST BE SQUARE IMAGE AND DIMZ MUST BE POWER OF 2)

% specifiy number of samples of each kernel subsample image you want to be
% used per output noise frame..
%Note: length must be = to the image dimension power of two minus one (ie if 512x512-->9, if 256x256-->8)
%imgsPerKrnl=[1,1,1,1,1,1,1,1,1]; % adjusts the total number of subtile images from each respective kernel incorporated into the noise
imgsPerKrnl=[1,1,1,1,1,1,0,0,0]; % adjusts the total number of subtile images from each respective kernel incorporated into the noise
krnlImgWgts=[1,1,1,1,1,1,1,1,1]; % adjusts relative weights applied to subtile images from each respective kernel
totKimz=sum(imgsPerKrnl); % count up total number of intermediate images..

smoothTiles=[0,0,0,0,0,0,0,0,0]; % apply smoothing kernel to intermediate subtile images?
gSmthK1sz=[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]; % size of gaussian smoothing kernel for subtiles..
smoothFinal=0; % apply smoothing kernel to final image?
gSmthK2sz=0.5; % size of gaussian smoothing kernel for final image..
BI_WT=0.35;
N_WT=1-BI_WT;

nReps=256; % number of noise image copies/reps you want to run..

nWorkers=14; % number of cpu nodes used in the parallelized portions..
%parProfile='Processes';
parProfile='HPC Cluster';
ClstrTimeStr='--time=00-01:00:00';
mlabComDirTmp="/home/eduwell/Matlab_Sandbox/imgKernelNoise/tmp";
%nReps=1; % number of noise image copies/reps you want to run..
%nWorkers=1; % number of cpu nodes used in the parallelized portions..

% BI reconstruction from rearranging noise kernel windows parameters:
reconBIfromNoise=1; % if 1 runs the Base Image reconstruction from noise portion. If 0, doesn't,
nRepsBIgen=256; % number of samples you want to use/reps you want to run in the base image generation from noise (must be <= nReps)
MinMaxCor=1; % specify whether you want to rearrange kernel windows to fall at locations min or max correlation value
nKernReArr=6; % set number of kernels you want to re-arrange..
pltBIrecon=1; % controls whether BI recon figs are created (1=plot, 0=don't)

% CI Generation Parameters
makeCIz=1; % if 1, runs the Classification Image generation portion...
CIVer="corr"; % "corr" = use correlation, "ssim" = use structural similarity
pltCIs=1; % controls whether CI figs are created (1=plot, 0=don't)
pltCI_hist=1; % controls whether plots of histograms w/ distributions of correlation/ssim values used of noise to BIs are produced(1=plot, 0=don't)
pltNzSmpls=1; % controls whether representative noise samples are plotted for each noise type (1=plot, 0=don't)

baseImgL="/home/eduwell/Matlab_Sandbox/imgKernelNoise/baseImgs/BaseImgSF_6L.png"; % path to image for BI 1 (for CI proceedure)
baseImgR="/home/eduwell/Matlab_Sandbox/imgKernelNoise/baseImgs/BaseImgSF_6R.png"; % path to image for BI 2 (for CI proceedure)
%baseImgL="/home/eduwell/Matlab_Sandbox/imgKernelNoise/baseImgs/londonFaces/smiling_front/103_08.jpg";%smiling
%baseImgR="/home/eduwell/Matlab_Sandbox/imgKernelNoise/baseImgs/londonFaces/neutral_front/103_03.jpg";%neutral

%% Subsample the image with kernelSubSampler

nzImz=uint8(zeros(desiredSize(1),desiredSize(2),nReps));

tic
for kk = 1:nReps

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
nzImz(:,:,kk)=noiseImg;
end
disp(" ")
disp("Making Kernel Noise Imz Took:")
toc

%% Generate other noise type images for comparison
tic
WhtNzImz=uint8(zeros(desiredSize(1),desiredSize(2),nReps));
for kk = 1:nReps
whiteNoiseImg=uint8(randi(255,[size(adjImg,1),size(adjImg,1)]));
WhtNzImz(:,:,kk)=whiteNoiseImg;
end
disp(" ")
disp("Making White Noise Imz Took:")
toc

%% Attempt to reconstruct the original base image using xCorKrnlNzReArrange for each noise type
if reconBIfromNoise == 1
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
%   - The extent to which the output images 'resemble' or correspond 
%     to the base image can be quantified by computing the structural
%     similarity (ssim) between them. This will give both an overall
%     summary value and the ssim pixel-wise across the image array to 
%     indicate which regions/image features are driving the overall ssim
%     for each kernel's "reconstructed base image" built from 
%     kernel-windowed noise tiles.
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

xCorOptNzImg_all=cell(1,nReps);
xCorOptNzImg2_all=cell(1,nReps);

xCorOptNzImg_ssimVal_all=cell(1,nReps);
xCorOptNzImg_ssimImg_all=cell(1,nReps);

xCorOptNzImg2_ssimVal_all=cell(1,nReps);
xCorOptNzImg2_ssimImg_all=cell(1,nReps);

indTallyImgReconz1_all=cell(1,nReps);
TlyImgRcn1_ssimImz_all=cell(1,nReps);
TlyImgRcn1_ssimValz_all=cell(1,nReps);

indTallyImgReconz2_all=cell(1,nReps);
TlyImgRcn2_ssimImz_all=cell(1,nReps);
TlyImgRcn2_ssimValz_all=cell(1,nReps);

% Open Local Parpool..
%parpool(parProfile,nWorkers);
%parpool('Processes',nWorkers);
startDir=pwd; % save starting location..
cd(mlabComDirTmp); % go to the tmp communication dir
addpath(pwd);
% Start the parpool in here so annoying files don't clog up main dir:
% Select cluster
c = parcluster(parProfile);
  
% Set time
%c.AdditionalProperties.AdditionalSubmitArgs = '--time=00-01:00:00';
c.AdditionalProperties.AdditionalSubmitArgs = ClstrTimeStr;

%for kk = 1:nReps
tic

%parfor (kk = 1:nReps,nWorkers)
parfor (kk = 1:nRepsBIgen,c)
% Run xCorKrnlNzReArrange_v2 on Tailored Kernel Noise
%tic
[xCorAdjNzFrames, KrnlXCorData, combTallyImg, indTallyImgz] = xCorKrnlNzReArrange_v3(adjImg, nzImz(:,:,kk), kernelDim_ImFrac, kernelSqrDims,MinMaxCor,nKernReArr);
%toc

% Run xCorKrnlNzReArrange_v2 on White Noise
%tic
[xCorAdjNzFrames2, KrnlXCorData2, combTallyImg2, indTallyImgz2] = xCorKrnlNzReArrange_v3(adjImg, WhtNzImz(:,:,kk), kernelDim_ImFrac, kernelSqrDims,MinMaxCor,nKernReArr);
%toc

% Generate "tally weighted" images for all kernels combined
% for both tailored noise and white noise
% Generate ssim images/values between these and the base image too..
xCorOptNzImg = weightByTallyImg(combTallyImg,xCorAdjNzFrames);
xCorOptNzImg2 = weightByTallyImg(combTallyImg2,xCorAdjNzFrames2);

[xCorOptNzImg_ssimVal, xCorOptNzImg_ssimImg] = ssim(xCorOptNzImg,adjImg);
[xCorOptNzImg2_ssimVal, xCorOptNzImg2_ssimImg] = ssim(xCorOptNzImg2,adjImg);

% Save into combined output..
xCorOptNzImg_all{kk}=xCorOptNzImg;
xCorOptNzImg2_all{kk}=xCorOptNzImg2;
xCorOptNzImg_ssimVal_all{kk}=xCorOptNzImg_ssimVal;
xCorOptNzImg_ssimImg_all{kk}=xCorOptNzImg_ssimImg;
xCorOptNzImg2_ssimVal_all{kk}=xCorOptNzImg2_ssimVal;
xCorOptNzImg2_ssimImg_all{kk}=xCorOptNzImg2_ssimImg;

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

% Save into combined output
indTallyImgReconz1_all{kk}=indTallyImgReconz1;
TlyImgRcn1_ssimImz_all{kk}=TlyImgRcn1_ssimImz;
TlyImgRcn1_ssimValz_all{kk}=TlyImgRcn1_ssimValz;
indTallyImgReconz2_all{kk}=indTallyImgReconz2;
TlyImgRcn2_ssimImz_all{kk}=TlyImgRcn2_ssimImz;
TlyImgRcn2_ssimValz_all{kk}=TlyImgRcn2_ssimValz;
end
disp(" ")
disp("Making Recon-BIs From Noise Took:")
toc
% Close down local parpool
delete(gcp('nocreate'));


%% Make Weighted Average Versions of Images Across nReps..

xCorOptNzImg2_all_mean=uint8(zeros(size(adjImg,1),size(adjImg,2),nRepsBIgen));
for kk=1:nRepsBIgen
xCorOptNzImg2_all_mean(:,:,kk)=xCorOptNzImg2_all{kk};
end
xCorOptNzImg2_all_mean=uint8(sum(xCorOptNzImg2_all_mean,3)./size(xCorOptNzImg2_all_mean,3));

xCorOptNzImg_all_mean=uint8(zeros(size(adjImg,1),size(adjImg,2),nRepsBIgen));
for kk=1:nRepsBIgen
xCorOptNzImg_all_mean(:,:,kk)=xCorOptNzImg_all{kk};
end
xCorOptNzImg_all_mean=uint8(sum(xCorOptNzImg_all_mean,3)./size(xCorOptNzImg_all_mean,3));

xCorOptNzImg_ssimImg_all_mean=zeros(size(adjImg,1),size(adjImg,2),size(xCorOptNzImg_ssimImg_all{1,1},3),nRepsBIgen);
for kk=1:nRepsBIgen
xCorOptNzImg_ssimImg_all_mean(:,:,:,kk)=xCorOptNzImg_ssimImg_all{kk};
end
xCorOptNzImg_ssimImg_all_mean=(sum(xCorOptNzImg_ssimImg_all_mean,4)./size(xCorOptNzImg_ssimImg_all_mean,4));

xCorOptNzImg2_ssimImg_all_mean=zeros(size(adjImg,1),size(adjImg,2),size(xCorOptNzImg2_ssimImg_all{1,1},3),nRepsBIgen);
for kk=1:nRepsBIgen
xCorOptNzImg2_ssimImg_all_mean(:,:,:,kk)=xCorOptNzImg2_ssimImg_all{kk};
end
xCorOptNzImg2_ssimImg_all_mean=(sum(xCorOptNzImg2_ssimImg_all_mean,4)./size(xCorOptNzImg2_ssimImg_all_mean,4));


indTallyImgReconz1_all_mean=zeros(size(adjImg,1),size(adjImg,2),size(indTallyImgReconz1_all{1,1},3),nRepsBIgen);
for kk=1:nRepsBIgen
indTallyImgReconz1_all_mean(:,:,:,kk)=indTallyImgReconz1_all{kk};
end
indTallyImgReconz1_all_mean=uint8((sum(indTallyImgReconz1_all_mean,4)./size(indTallyImgReconz1_all_mean,4)));

indTallyImgReconz2_all_mean=zeros(size(adjImg,1),size(adjImg,2),size(indTallyImgReconz2_all{1,1},3),nRepsBIgen);
for kk=1:nRepsBIgen
indTallyImgReconz2_all_mean(:,:,:,kk)=indTallyImgReconz2_all{kk};
end
indTallyImgReconz2_all_mean=uint8((sum(indTallyImgReconz2_all_mean,4)./size(indTallyImgReconz2_all_mean,4)));

TlyImgRcn1_ssimImz_all_mean=zeros(size(adjImg,1),size(adjImg,2),size(TlyImgRcn1_ssimImz_all{1,1},3),nRepsBIgen);
for kk=1:nRepsBIgen
TlyImgRcn1_ssimImz_all_mean(:,:,:,kk)=TlyImgRcn1_ssimImz_all{kk};
end
TlyImgRcn1_ssimImz_all_mean=(sum(TlyImgRcn1_ssimImz_all_mean,4)./size(TlyImgRcn1_ssimImz_all_mean,4));

TlyImgRcn2_ssimImz_all_mean=zeros(size(adjImg,1),size(adjImg,2),size(TlyImgRcn2_ssimImz_all{1,1},3),nRepsBIgen);
for kk=1:nRepsBIgen
TlyImgRcn2_ssimImz_all_mean(:,:,:,kk)=TlyImgRcn2_ssimImz_all{kk};
end
TlyImgRcn2_ssimImz_all_mean=(sum(TlyImgRcn2_ssimImz_all_mean,4)./size(TlyImgRcn2_ssimImz_all_mean,4));

TlyImgRcn1_ssimValz_all_mean=zeros(size(TlyImgRcn1_ssimValz_all{1,1},1),nRepsBIgen);
for kk=1:nRepsBIgen
TlyImgRcn1_ssimValz_all_mean(:,kk)=TlyImgRcn1_ssimValz_all{kk}(:,1);
end
TlyImgRcn1_ssimValz_all_mean_std=std(TlyImgRcn1_ssimValz_all_mean,0,2);
TlyImgRcn1_ssimValz_all_mean=mean(TlyImgRcn1_ssimValz_all_mean,2);

TlyImgRcn2_ssimValz_all_mean=zeros(size(TlyImgRcn2_ssimValz_all{1,1},1),nRepsBIgen);
for kk=1:nRepsBIgen
TlyImgRcn2_ssimValz_all_mean(:,kk)=TlyImgRcn2_ssimValz_all{kk}(:,1);
end
TlyImgRcn2_ssimValz_all_mean_std=std(TlyImgRcn2_ssimValz_all_mean,0,2);
TlyImgRcn2_ssimValz_all_mean=mean(TlyImgRcn2_ssimValz_all_mean,2);

xCorOptNzImg_ssimVal_all_mean=zeros(1,nRepsBIgen);
for kk=1:nRepsBIgen
xCorOptNzImg_ssimVal_all_mean(1,kk)=xCorOptNzImg_ssimVal_all{kk};
end
xCorOptNzImg_ssimVal_all_mean_std=std(xCorOptNzImg_ssimVal_all_mean,0,2);
xCorOptNzImg_ssimVal_all_mean=mean(xCorOptNzImg_ssimVal_all_mean,2);

xCorOptNzImg2_ssimVal_all_mean=zeros(1,nRepsBIgen);
for kk=1:nRepsBIgen
xCorOptNzImg2_ssimVal_all_mean(1,kk)=xCorOptNzImg2_ssimVal_all{kk};
end
xCorOptNzImg2_ssimVal_all_mean_std=std(xCorOptNzImg2_ssimVal_all_mean,0,2);
xCorOptNzImg2_ssimVal_all_mean=mean(xCorOptNzImg2_ssimVal_all_mean,2);

%% Make BI/Noise Combo Images

noiseBICombo=((BI_WT.*adjImg+N_WT.*noiseImg)./2);
noiseBICombo=uint8(rescale(noiseBICombo,0,255));
whtnoiseBICombo=((BI_WT.*adjImg+N_WT.*whiteNoiseImg)./2);
whtnoiseBICombo=uint8(rescale(whtnoiseBICombo,0,255));
xCorOptNzImgBICombo=((BI_WT.*adjImg+N_WT.*xCorOptNzImg_all{1})./2);
xCorOptNzImgBICombo=uint8(rescale(xCorOptNzImgBICombo,0,255));

%% Plot Images
if pltBIrecon==1
% PLOT (2)
% Reconstruction of Base Image with Tailored vs. White Noise using all
% kernel windows combined..
%##########################################################################
figure
hold on
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
imshow(xCorOptNzImg2_all_mean);
title("BI Reconstructed with White Noise")

subplot(2,6,[10,11]);
imshow(xCorOptNzImg_all_mean);
title("BI Reconstructed with Tailored Kernel Noise")
hold off
shg
%##########################################################################

% PLOT (3)
% Reconstruction of Base Image with Tailored vs. White Noise using each
% individual kernel window..
%##########################################################################
figure
shg
hold on
%pChkSz=round(nKernReArr/3);

nKernReArr2=nKernReArr;
if nKernReArr<6
nKernReArr2=nKernReArr2*2;
fillerNum=1;
else
fillerNum=0;
end
IS_EVEN = ~mod(nKernReArr2,2);
pMid = round(nKernReArr2/2);


if IS_EVEN == 0
subplot(3,nKernReArr2,1);
imshow(adjImg);
title("Base Image")
else
subplot(3,nKernReArr2,[1,2]);
imshow(adjImg);
title("Base Image")
end

if IS_EVEN == 0
subplot(3,nKernReArr2,pMid);
imshow(whiteNoiseImg);
title("White Noise")
else
subplot(3,nKernReArr2,[pMid,pMid+1]);
imshow(whiteNoiseImg);
title("White Noise") 
end

if IS_EVEN == 0
subplot(3,nKernReArr2,nKernReArr2);
imshow(noiseImg);
title("Tailored 'Kernel' Noise")
else
subplot(3,nKernReArr2,[nKernReArr2-1,nKernReArr2]);
imshow(noiseImg);
title("Tailored 'Kernel' Noise")
end

pstart=nKernReArr2+1; % initialize..
for uu = 1:nKernReArr
    
    subplot(3,nKernReArr2,[pstart,(pstart+fillerNum)]);
    im2Show=indTallyImgReconz2_all_mean(:,:,uu);
    imshow(im2Show);
    titleTmp=strcat("White Nz BI Recon With K",num2str(uu));
    title(titleTmp)
    pstart=(pstart+1+fillerNum);
end

for uu = 1:nKernReArr
    
    subplot(3,nKernReArr2,[pstart,(pstart+fillerNum)]);
    im2Show=indTallyImgReconz1_all_mean(:,:,uu);
    imshow(im2Show);
    titleTmp=strcat("Tailored Nz BI Recon With K",num2str(uu));
    title(titleTmp)
    pstart=(pstart+1+fillerNum);
end
hold off
shg

% PLOT (4)
% SSIM of Reconstructed Base Images using Tailored vs. White Noise for each
% individual kernel window size..
%##########################################################################
figure
shg
hold on

nKernReArr2=nKernReArr;
if nKernReArr<6
nKernReArr2=nKernReArr2*2;
fillerNum=1;
else
fillerNum=0;
end
IS_EVEN = ~mod(nKernReArr2,2);
pMid = round(nKernReArr2/2);

if IS_EVEN == 0
subplot(3,nKernReArr2,1);
imshow(adjImg);
title("Base Image")
else
subplot(3,nKernReArr2,[1,2]);
imshow(adjImg);
title("Base Image")
end

if IS_EVEN == 0
subplot(3,nKernReArr2,pMid);
imshow(uint8(rescale(xCorOptNzImg2_ssimImg_all_mean,0,255)));
title("White Noise SSIM to BI: All Kernels Combined")
else
subplot(3,nKernReArr2,[pMid,pMid+1]);
imshow(uint8(rescale(xCorOptNzImg2_ssimImg_all_mean,0,255)));
title("White Noise SSIM to BI: All Kernels Combined")  
end

if IS_EVEN == 0
subplot(3,nKernReArr2,nKernReArr2);
imshow(uint8(rescale(xCorOptNzImg_ssimImg_all_mean,0,255)));
title("Tailored Noise SSIM to BI: All Kernels Combined")
else
subplot(3,nKernReArr2,[nKernReArr2-1,nKernReArr2]);
imshow(uint8(rescale(xCorOptNzImg_ssimImg_all_mean,0,255)));
title("Tailored Noise SSIM to BI: All Kernels Combined")
end

pstart=nKernReArr2+1; % initialize..
for uu = 1:nKernReArr

    subplot(3,nKernReArr2,[pstart,(pstart+fillerNum)]);
    im2Show=TlyImgRcn2_ssimImz_all_mean(:,:,uu);
    imshow(im2Show);
    titleTmp=strcat("SSIM White w/ K",num2str(uu)," and BI");
    title(titleTmp)
    pstart=(pstart+1+fillerNum);
end

for uu = 1:nKernReArr

    subplot(3,nKernReArr2,[pstart,(pstart+fillerNum)]);
    im2Show=TlyImgRcn1_ssimImz_all_mean(:,:,uu);
    imshow(im2Show);
    titleTmp=strcat("SSIM Tailored w/ K",num2str(uu)," and BI");
    title(titleTmp)
    pstart=(pstart+1+fillerNum);
end
hold off
shg

% PLOT (5)
% Bar Plots of Overall SSIM of Reconstructed Base Images using Tailored 
% vs. White Noise for each individual kernel window size..
%##########################################################################

% Preallocate..
xlabelz = num2cell(zeros(1, (nKernReArr+1)));
%nTix=1:(nKernReArr+1):1;
nTix=linspace(1,nKernReArr+1,nKernReArr+1);
ssimValz2Plot=zeros(2,nKernReArr+1);
ssimValz2Plot_std=zeros(2,nKernReArr+1);

countr=1;
for ff = 1:nKernReArr+1

    % for the last one label as "all-Ks-combined"
    % and plot the combined kernel overall ssim values..
    if ff == (nKernReArr+1)
        xlabelz(1,countr)={'all-Ks-combined'};
        ssimValz2Plot(1,ff)=xCorOptNzImg2_ssimVal_all_mean;
        ssimValz2Plot(2,ff)=xCorOptNzImg_ssimVal_all_mean;

        ssimValz2Plot_std(1,ff)=xCorOptNzImg2_ssimVal_all_mean_std;
        ssimValz2Plot_std(2,ff)=xCorOptNzImg_ssimVal_all_mean_std;

    else
        lablTmp=strcat('kernel',num2str(ff));
        xlabelz(1,countr)={lablTmp};
        ssimValz2Plot(1,ff)=TlyImgRcn2_ssimValz_all_mean(ff);
        ssimValz2Plot(2,ff)=TlyImgRcn1_ssimValz_all_mean(ff);

        ssimValz2Plot_std(1,ff)=TlyImgRcn2_ssimValz_all_mean_std(ff);
        ssimValz2Plot_std(2,ff)=TlyImgRcn1_ssimValz_all_mean_std(ff);
    end
    countr=countr+1;
end


figure
shg
hold on
%   Symmetric Example:
h = barwitherr(ssimValz2Plot_std', ssimValz2Plot');% Plot with errorbars
%ax=gca;
%ax.XTickLabel=xlabelz;
x0=10;
y0=10;
width=1600;
height=1000;
set(gcf,'position',[x0,y0,width,height]);
ax=gca;
xticks( nTix );
set(gca,'XTickLabel',xlabelz)
xtickangle(35);
% labels = string(ax.XAxis.TickLabels); % extract
% labels(2:2:end) = nan; % remove every other one
% ax.XAxis.TickLabels = labels; % set

%xtickangle(45)
legend('White Noise','Tailored Kernel Noise')
ylabel("Overall Structural Similarity with BI");
xlabel("Kernel Size Used During Recon Proceedure");
set(h(1),'FaceColor',[0.9 0.9 0.9]);
set(h(2),'FaceColor',[0.5 0.5 0.5]);
%bar(xlabelz,ssimValz2Plot);
title("SSIM Between BI and Recon-BIs Created with White vs. Tailored Noise by Kernel Size");
hold off

shg
%hold off
%##########################################################################

cd(startDir); % return to startpoint..

end
end

%% Now generate CIs with the white and tailored noise samples..
if makeCIz == 1
tic
%% CI Generation Parameters

rng("shuffle");


%% Read in base images
BI_L=imread(baseImgL);
BI_R=imread(baseImgR);

% Convert to desired dimensions..
BI_L = imgReformater(BI_L,desiredSize);
BI_R = imgReformater(BI_R,desiredSize);

%% Combine White and Tailored Noise Images Into Same Array

noiseImsComb=cell(1,2);
noiseImsComb(1,1)={nzImz};
noiseImsComb(1,2)={WhtNzImz};

%% For each noise image, compute structural similarity between noise image and each base image

for ww=1:size(noiseImsComb,2)

ssmat=cell(6,size(noiseImsComb{1,ww},3));

for uu = 1:size(noiseImsComb{1,ww},3)
    %img = noiseImsComb{1,uu};
    img = noiseImsComb{1,ww}(:,:,uu);
    %[ssimvalL,ssimmapL]  = ssim(img,BI_L);
    %[ssimvalR,ssimmapR]  = ssim(img,BI_R);
    
    if CIVer == "ssim"
    ssimvalL = ssim(img,BI_L);
    ssimvalR = ssim(img,BI_R);
    end

    if CIVer == "corr"
    ssimvalL = corr2(img,BI_L);
    ssimvalR = corr2(img,BI_R);
    end

    ssimvalDiff=ssimvalL-ssimvalR;
    
    ssmat{1,uu}=ssimvalL;
    ssmat{2,uu}=ssimvalR;
    %ssmat{3,uu}=ssimmapL;
    %ssmat{4,uu}=ssimmapR;
    ssmat{5,uu}=ssimvalDiff;
    %ssimvalDiffMap=ssimmapL-ssimmapR;
    %ssmat{6,uu}=ssimvalDiffMap;

end

%% Get the max and min ssim vals for both the the right and left base images
ssimvalL_min=min(cell2mat(ssmat(1,:)));
ssimvalL_max=max(cell2mat(ssmat(1,:)));
ssimvalR_min=min(cell2mat(ssmat(2,:)));
ssimvalR_max=max(cell2mat(ssmat(2,:)));

ssimvalDiff_min=min(cell2mat(ssmat(5,:)));
ssimvalDiff_max=max(cell2mat(ssmat(5,:)));

%% Use Linspace to set up a 100 value vector ranging from ssimvalDiff_min to ssimvalDiff_max

pMin=-max(max(abs(cell2mat(ssmat(5,:)))));
pMax=max(max(abs(cell2mat(ssmat(5,:)))));

x = linspace(-10, 10, 100);
%x = linspace(-5, 5, 100);
x2 = linspace(pMin, pMax, 100);

% Psychometric Curve..
y1= 1./(1+exp(-x));

difValProbMat=zeros(2,100); % intitialize
difValProbMat(1,:)=linspace(ssimvalDiff_min,ssimvalDiff_max,100); % ssim difference
difValProbMat(2,:)=y1; % corresponding response probability (Right vs. Left);

figure
plot(difValProbMat(1,:),difValProbMat(2,:));

%% Combine noiseImsComb and ssmat arrays images into same array
noiseImsComb_cell=cell(1,size(ssmat,2));
for ll = 1:size(ssmat,2)
noiseImsComb_cell{1,ll}=noiseImsComb{1,ww}(:,:,ll);
end
%noiseImsComb_ssmat=vertcat(noiseImsComb,ssmat);
noiseImsComb_ssmat=vertcat(noiseImsComb_cell,ssmat);

%% Loop through the ssimvalDiff values for each noise image and assign a response based on probability vector values
respMat=zeros(2,size(noiseImsComb_cell,2));
respOps=[0,1]; % 0 = Left, 1 = R;
for zz=1:size(noiseImsComb_cell,2)

    % get ssim difference val for this noise image
    ssimDif=ssmat{5,zz};

    % find the closest matching value in the difValProbMat array (row 1) and its
    % index position.
    [closestVal,idx] = findClosest(ssimDif, difValProbMat(1,:));

    % use the index to grab the corresponding probablility assigned in the
    % second row of difValProbMat. Use this to compute pRight and pLeft
    pLeft=difValProbMat(2,idx);
    pRight=1-pLeft;
%     pRight=difValProbMat(2,idx);
%     pLeft=1-pRight;

    % feed these probabilites into randsample to select a random response
    resp = randsample(respOps,1,true,[pLeft, pRight]);

    % save response in the output respmat
    respMat(1,zz)=resp;

    % feed these probabilites into randsample to select a random condition
    tcon = randsample(respOps,1,true,[pLeft, pRight]);

    respMat(2,zz)=tcon;

end

%% Add response data to noiseImsComb_ssmat
respMat=num2cell(respMat);
noiseImsComb_ssmat=vertcat(noiseImsComb_ssmat,respMat);

%% Assign images to CI_11 --> CI_22 based in response/tcon combo..

% initialize
CI_11raw={};
CI_12raw={};
CI_21raw={};
CI_22raw={};


for hh = 1:size(noiseImsComb_ssmat,2)

    if respMat{1,hh}==0 && respMat{2,hh}==0
        CI_11raw{end+1}=noiseImsComb_cell{1,hh}; 

    end

    if respMat{1,hh}==0 && respMat{2,hh}==1
        CI_21raw{end+1}=noiseImsComb_cell{1,hh};
     
    end

    if respMat{1,hh}==1 && respMat{2,hh}==0
        CI_12raw{end+1}=noiseImsComb_cell{1,hh};

    end

    if respMat{1,hh}==1 && respMat{2,hh}==1
        CI_22raw{end+1}=noiseImsComb_cell{1,hh};        
    end

end

% Normal Images
CI_11rawmat=zeros(size(CI_11raw{1,1},1),size(CI_11raw{1,1},2),length(CI_11raw));
for jj = 1:length(CI_11raw)
CI_11rawmat(:,:,jj)=CI_11raw{1,jj};
end
clear jj;

CI_21rawmat=zeros(size(CI_21raw{1,1},1),size(CI_21raw{1,1},2),length(CI_21raw));
for jj = 1:length(CI_21raw)
CI_21rawmat(:,:,jj)=CI_21raw{1,jj};
end
clear jj;

CI_12rawmat=zeros(size(CI_12raw{1,1},1),size(CI_12raw{1,1},2),length(CI_12raw));
for jj = 1:length(CI_12raw)
CI_12rawmat(:,:,jj)=CI_12raw{1,jj};
end
clear jj;

CI_22rawmat=zeros(size(CI_22raw{1,1},1),size(CI_22raw{1,1},2),length(CI_22raw));
for jj = 1:length(CI_22raw)
CI_22rawmat(:,:,jj)=CI_22raw{1,jj};
end
clear jj;


%% --- Compute CIs --- %

%"normal" version
CI_results1raw = (mean(CI_11rawmat,3) + mean(CI_21rawmat,3)) - (mean(CI_12rawmat,3) + mean(CI_22rawmat,3));
CI_results2raw = (mean(CI_12rawmat,3) + mean(CI_22rawmat,3)) - (mean(CI_11rawmat,3) + mean(CI_21rawmat,3));

%make smoothed version..
CI_results1rawSM=imgaussfilt(CI_results1raw,2);
CI_results2rawSM=imgaussfilt(CI_results2raw,2);


if pltCIs==1
%figure;
f = figure;
f.Position = [100 100 900 750];
subplot(2,2,1);
imagesc(BI_L);
colormap("gray");
axis equal
axis tight
colorbar
title("Orig Base Im1")

subplot(2,2,2);
imagesc(BI_R);
colormap("gray");
axis equal
axis tight
colorbar
title("Orig Base Im2")

subplot(2,2,3);
imagesc(CI_results1rawSM);
colormap("gray");
axis equal
axis tight
colorbar
title("CI 1")

subplot(2,2,4);
imagesc(CI_results2rawSM);
colormap("gray");
axis equal
axis tight
colorbar
title("CI 2")
end

%% Histogram of SSIM Difference Values

if pltCI_hist==1
figure
histogram(cell2mat(ssmat(5,:)));
end

end
disp(" ")
disp("Making CIs took:")
toc

%% Show samples of noises..

if pltNzSmpls==1
knSample=nzImz(:,:,1);
wnSample=WhtNzImz(:,:,1);

figure
subplot(2,1,1);
imagesc(knSample);
colormap("gray");
axis equal
axis tight
colorbar
title("Kernel Noise Sample")

subplot(2,1,2);
imagesc(wnSample);
colormap("gray");
axis equal
axis tight
colorbar
title("White Noise Noise Sample")
end

end
