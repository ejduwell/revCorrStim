%% Parameters
difVersion=0;
kernel_size=32;
%nWorkers=8;
%kernel_sizes=linspace(6,120,3);


% Open Local Parpool..
% parpool('Processes',nWorkers);



% NOISE SELECTION

%matFile="1726705_20231010T003550_AML_127.mat"; % white noise 5000
matFile="3616781_20231009T223655_AML_127.mat"; % glob noise 5000
%matFile="9330047_20231010T125132_AML_127.mat"; % glob noise 10000
%matFile="5305486_20231010T003408_AML_127.mat"; % grating noise 5000 sf18
%matFile="9738144_20231010T002719_AML_127.mat"; % grating noise 5000 sf10
%matFile="7487775_20231010T002026_AML_127.mat"; % grating noise 5000 sf2


rng("shuffle");
fSpaceRad=20;
filter="hf";

% BASE IMAGES
% 18 cycles per image
%baseImgL="/media/eduwell/Ubuntu/home/eduwell/Matlab_Sandbox/ssim_CIgen/baseImg/BaseImgSF_18L.png";
%baseImgR="/media/eduwell/Ubuntu/home/eduwell/Matlab_Sandbox/ssim_CIgen/baseImg/BaseImgSF_18R.png";

% 10 cycles per image
%baseImgL="/media/eduwell/Ubuntu/home/eduwell/Matlab_Sandbox/ssim_CIgen/baseImg/BaseImgSF_10L.png";
%baseImgR="/media/eduwell/Ubuntu/home/eduwell/Matlab_Sandbox/ssim_CIgen/baseImg/BaseImgSF_10R.png";

% 4 cycles per image
%baseImgL="/media/eduwell/Ubuntu/home/eduwell/Matlab_Sandbox/ssim_CIgen/baseImg/BaseImgSF_4L.png";
%baseImgR="/media/eduwell/Ubuntu/home/eduwell/Matlab_Sandbox/ssim_CIgen/baseImg/BaseImgSF_4R.png";

% 
% % 2 cycles per image
baseImgL="/media/eduwell/Ubuntu/home/eduwell/Matlab_Sandbox/ssim_CIgen/baseImg/BaseImgSF_2L.png";
baseImgR="/media/eduwell/Ubuntu/home/eduwell/Matlab_Sandbox/ssim_CIgen/baseImg/BaseImgSF_2R.png";

% baseImgL="/media/eduwell/Ubuntu/home/eduwell/Matlab_Sandbox/ssim_CIgen/baseImg/gabLeft.png";
% baseImgR="/media/eduwell/Ubuntu/home/eduwell/Matlab_Sandbox/ssim_CIgen/baseImg/gabRight.png";

% FACES
%baseImgL="/media/eduwell/Ubuntu/home/eduwell/Matlab_Sandbox/ssim_CIgen/face/107_03_rfmt.jpg";
%baseImgR="/media/eduwell/Ubuntu/home/eduwell/Matlab_Sandbox/ssim_CIgen/face/107_08_rfmt.jpg";
% baseImgL="/media/eduwell/Ubuntu/home/eduwell/Matlab_Sandbox/ssim_CIgen/face/031_03_rfmt.jpg";
% baseImgR="/media/eduwell/Ubuntu/home/eduwell/Matlab_Sandbox/ssim_CIgen/face/031_08_rfmt.jpg";

%REVCORR
baseImgL="/media/eduwell/Ubuntu/home/eduwell/Matlab_Sandbox/ssim_CIgen/baseImg/revCorrStim/left_rfmt.png";
baseImgR="/media/eduwell/Ubuntu/home/eduwell/Matlab_Sandbox/ssim_CIgen/baseImg/revCorrStim/right_rfmt.png";
%% Read In Noise Images
data = load(matFile);
noiseIms=data.noiseImgMat;

%% Read in base images
BI_L=imread(baseImgL);
BI_R=imread(baseImgR);

%% Take 2dFFT of of each

fftFilter=0; % make the fft mask filtered noise images?
noiseImsFFT=cell(2,size(noiseIms,2));

if fftFilter==1
for ii =1:size(noiseIms,2)

    % take 2dfft ...
    %noiseImsFFT{1,ii}=fft2(noiseIms{1,ii});
    
    % run fftMaskFilter
    noiseImsFFT{1,ii}= fftMaskFilter(noiseIms{1,ii},fSpaceRad,filter);

    % For randomization of order later..
    noiseImsFFT{2,ii}=rand(1); % assign random number

end
end

%% Combine normal and fftMask filtered images into same array
noiseImsComb=vertcat(noiseIms,noiseImsFFT);


%% For each noise image, compute structural similarity between noise image and each base image
ssmat=cell(6,size(noiseImsComb,2));

%parfor (uu = 1:size(noiseImsComb,2),nWorkers)
for uu = 1:size(noiseImsComb,2)
    img = noiseImsComb{1,uu};
    cL = corr2(img,BI_L);
    cR = corr2(img,BI_R);
%     corrImL = localCorrelation_v3(img, BI_L, kernel_sizes);
%     corrImR = localCorrelation_v3(img, BI_R, kernel_sizes);
    corrImL = localCorrelation_v2(img, BI_L, kernel_size);
    corrImR = localCorrelation_v2(img, BI_R, kernel_size);

    %[ssimvalL,ssimmapL]  = ssim(img,BI_L);
    %[ssimvalR,ssimmapR]  = ssim(img,BI_R);
    corValDiff=cL-cR;
    
    ssmat{1,uu}=cL;
    ssmat{2,uu}=cR;
    ssmat{3,uu}=corrImL;
    ssmat{4,uu}=corrImR;
    ssmat{5,uu}=corValDiff;
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

% Psychometric Curve Condition..
y1= 1./(1+exp(-x));
%UNCOMMENT FOR NULL CONDITION (50% chance response across board..)
%y1=0.5.*(ones(1,100));

difValProbMat=zeros(2,100); % intitialize

%difValProbMat(1,:)=linspace(ssimvalDiff_min,ssimvalDiff_max,100); % ssim difference
difValProbMat(1,:)=linspace(ssimvalDiff_min,ssimvalDiff_max,100); % ssim difference

%difValProbMat(2,:)=linspace(0,1,100); % corresponding response probability (Right vs. Left);
difValProbMat(2,:)=y1; % corresponding response probability (Right vs. Left);

figure
plot(difValProbMat(1,:),difValProbMat(2,:));

%% Combine noiseImsComb and ssmat arrays images into same array
noiseImsComb_ssmat=vertcat(noiseImsComb,ssmat);

%% Loop through the ssimvalDiff values for each noise image and assign a response based on probability vector values
respMat=zeros(2,size(noiseImsComb,2));
respOps=[0,1]; % 0 = Left, 1 = R;
for zz=1:size(noiseImsComb,2)

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

    if resp == 0
        ssimvalDiffMap=ssmat{3,zz}-ssmat{4,zz};
        ssmat{6,zz}=ssimvalDiffMap;
    elseif resp == 1
        ssimvalDiffMap=ssmat{4,zz}-ssmat{3,zz};
        ssmat{6,zz}=ssimvalDiffMap;
    end
end

%% Add response data to noiseImsComb_ssmat
respMat=num2cell(respMat);
noiseImsComb_ssmat=vertcat(noiseImsComb_ssmat,respMat);

%% Assign images to CI_11 --> CI_22 based in response/tcon combo..

% initialize
CI_11={};
CI_12={};
CI_21={};
CI_22={};

CI_11raw={};
CI_12raw={};
CI_21raw={};
CI_22raw={};

if difVersion==1
for hh = 1:size(noiseImsComb_ssmat,2)

    if respMat{1,hh}==0 && respMat{2,hh}==0
        CI_11{end+1}=ssmat{6,hh};
        %CI_11{end+1}=ssmat{3,hh}; 
        CI_11raw{end+1}=noiseImsComb{1,hh}; 

    end

    if respMat{1,hh}==0 && respMat{2,hh}==1
        CI_21{end+1}=ssmat{6,hh}; 
        %CI_21{end+1}=ssmat{3,hh};
        CI_21raw{end+1}=noiseImsComb{1,hh};
     
    end

    if respMat{1,hh}==1 && respMat{2,hh}==0
        CI_12{end+1}=ssmat{6,hh};
        %CI_12{end+1}=ssmat{4,hh};
        CI_12raw{end+1}=noiseImsComb{1,hh};

    end

    if respMat{1,hh}==1 && respMat{2,hh}==1
        CI_22{end+1}=ssmat{6,hh};
        %CI_22{end+1}=ssmat{4,hh};
        CI_22raw{end+1}=noiseImsComb{1,hh};
        
    end

end

else

for hh = 1:size(noiseImsComb_ssmat,2)

    if respMat{1,hh}==0 && respMat{2,hh}==0
        %CI_11{end+1}=ssmat{6,hh};
        CI_11{end+1}=ssmat{3,hh}; 
        CI_11raw{end+1}=noiseImsComb{1,hh}; 

    end

    if respMat{1,hh}==0 && respMat{2,hh}==1
        %CI_21{end+1}=ssmat{6,hh}; 
        CI_21{end+1}=ssmat{3,hh};
        CI_21raw{end+1}=noiseImsComb{1,hh};
     
    end

    if respMat{1,hh}==1 && respMat{2,hh}==0
        %CI_12{end+1}=ssmat{6,hh};
        CI_12{end+1}=ssmat{4,hh};
        CI_12raw{end+1}=noiseImsComb{1,hh};

    end

    if respMat{1,hh}==1 && respMat{2,hh}==1
        %CI_22{end+1}=ssmat{6,hh};
        CI_22{end+1}=ssmat{4,hh};
        CI_22raw{end+1}=noiseImsComb{1,hh};
        
    end

end
end

% Normal Images
CI_11mat=zeros(size(CI_11{1,1},1),size(CI_11{1,1},2),length(CI_11));
CI_11rawmat=zeros(size(CI_11{1,1},1),size(CI_11{1,1},2),length(CI_11));
for jj = 1:length(CI_11)
CI_11mat(:,:,jj)=CI_11{1,jj};
CI_11rawmat(:,:,jj)=CI_11raw{1,jj};
end
clear jj;


CI_21mat=zeros(size(CI_21{1,1},1),size(CI_21{1,1},2),length(CI_21));
CI_21rawmat=zeros(size(CI_21{1,1},1),size(CI_21{1,1},2),length(CI_21));
for jj = 1:length(CI_21)
CI_21mat(:,:,jj)=CI_21{1,jj};
CI_21rawmat(:,:,jj)=CI_21raw{1,jj};
end
clear jj;

CI_12mat=zeros(size(CI_12{1,1},1),size(CI_12{1,1},2),length(CI_12));
CI_12rawmat=zeros(size(CI_12{1,1},1),size(CI_12{1,1},2),length(CI_12));
for jj = 1:length(CI_12)
CI_12mat(:,:,jj)=CI_12{1,jj};
CI_12rawmat(:,:,jj)=CI_12raw{1,jj};
end
clear jj;

CI_22mat=zeros(size(CI_22{1,1},1),size(CI_22{1,1},2),length(CI_22));
CI_22rawmat=zeros(size(CI_22{1,1},1),size(CI_22{1,1},2),length(CI_22));
for jj = 1:length(CI_22)
CI_22mat(:,:,jj)=CI_22{1,jj};
CI_22rawmat(:,:,jj)=CI_22raw{1,jj};
end
clear jj;




%% --- Compute CIs --- %

% ssim version
CI_results1 = (mean(CI_11mat,3) + mean(CI_21mat,3)) - (mean(CI_12mat,3) + mean(CI_22mat,3));
CI_results2 = (mean(CI_12mat,3) + mean(CI_22mat,3)) - (mean(CI_11mat,3) + mean(CI_21mat,3));
% CI_results1=CI_results1*mask;
% CI_results2=CI_results2*mask;

%make smoothed version..
CI_results1SM=imgaussfilt(CI_results1,2);
CI_results2SM=imgaussfilt(CI_results2,2);

%"normal" version
CI_results1raw = (mean(CI_11rawmat,3) + mean(CI_21rawmat,3)) - (mean(CI_12rawmat,3) + mean(CI_22rawmat,3));
CI_results2raw = (mean(CI_12rawmat,3) + mean(CI_22rawmat,3)) - (mean(CI_11rawmat,3) + mean(CI_21rawmat,3));
% CI_results1raw=CI_results1raw*mask;
% CI_results2raw=CI_results2raw*mask;

%make smoothed version..
CI_results1rawSM=imgaussfilt(CI_results1raw,2);
CI_results2rawSM=imgaussfilt(CI_results2raw,2);

%figure;
f = figure;
f.Position = [100 100 2000 750];
subplot(2,5,1);
imagesc(BI_L);
colormap("gray");
axis equal
axis tight
colorbar
title("Orig Base Im1")

subplot(2,5,6);
imagesc(BI_R);
colormap("gray");
axis equal
axis tight
colorbar
title("Orig Base Im2")

subplot(2,5,2);
%imagesc(CI_results1raw);
imagesc(CI_results1rawSM);
colormap("gray");
axis equal
axis tight
colorbar
title("Raw CI 1")

%figure;
subplot(2,5,7);
%imagesc(CI_results2raw);
imagesc(CI_results2rawSM);
colormap("gray");
axis equal
axis tight
colorbar
title("Raw CI 2")


%figure;
subplot(2,5,3);
%imagesc(CI_results1);
imagesc(CI_results1SM);
colormap("gray");
axis equal
axis tight
colorbar
title("LocalCorr CI 1")

%figure;
subplot(2,5,8);
%imagesc(CI_results2);
imagesc(CI_results2SM);
colormap("gray");
axis equal
axis tight
colorbar
title("LocalCorr CI 2")


%figure;
subplot(2,5,4);
imagesc(mean(CI_11mat,3));
colormap("gray");
axis equal
axis tight
colorbar
title("LocalCorr CI-11 (mean)")

%figure;
subplot(2,5,5);
imagesc(mean(CI_21mat,3));
colormap("gray");
axis equal
axis tight
colorbar
title("LocalCorr CI-21 (mean)")

%figure;
subplot(2,5,9);
imagesc(mean(CI_12mat,3));
colormap("gray");
axis equal
axis tight
colorbar
title("LocalCorr CI-12 (mean)")

%figure;
subplot(2,5,10);
imagesc(mean(CI_22mat,3));
colormap("gray");
axis equal
axis tight
colorbar
title("LocalCorr CI-22 (mean)")

%% Nth percentile thresholding..

% Define the percentile
n = 90; % Change this to your desired value

%CI_results1_IC = imcomplement(CI_results1);
%CI_results2_IC = imcomplement(CI_results2);

% Determine the threshold
threshold_value1 = prctile(CI_results1(:), n);
threshold_value2 = prctile(CI_results2(:), n);

% Threshold the image
CI_results1_Nth = CI_results1 > threshold_value1;
CI_results2_Nth = CI_results2 > threshold_value2;

% Display the original and thresholded image
figure;
subplot(2,2,1);
imagesc(BI_L);
colormap("gray");
axis equal
axis tight
title("Orig Base Im1")

subplot(2,2,3);
imagesc(BI_R);
colormap("gray");
axis equal
axis tight
title("Orig Base Im2")

subplot(2, 2, 2);
imshow(CI_results1_Nth);
axis equal
axis tight
title(strcat("CIresults1: ",num2str(n),"th"," percentile"));

subplot(2, 2, 4);
imshow(CI_results2_Nth);
axis equal
axis tight
title(strcat("CIresults2: ",num2str(n),"th"," percentile"));



%% Histogram of SSIM Difference Values
figure
histogram(cell2mat(ssmat(5,:)));

%% Close down local parpool
delete(gcp('nocreate'));% Close down local parpool