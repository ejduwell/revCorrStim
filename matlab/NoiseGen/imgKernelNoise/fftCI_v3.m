%% Parameters
matFile="189542_20231009T185227.mat";
rng("shuffle");
fSpaceRad=20;
filter="hf";
baseImgL="/MATLAB Drive/sandbox/fftFiltered_CIgen/baseImg/gabLeft.png";
baseImgR="/MATLAB Drive/sandbox/fftFiltered_CIgen/baseImg/gabRight.png";

%% Read In Noise Images
data = load(matFile);
noiseIms=data.noiseImgMat;

%% Read in base images
BI_L=imread(baseImgL);
BI_R=imread(baseImgR);

%% Take 2dFFT of of each
noiseImsFFT=cell(2,size(noiseIms,2));

for ii =1:size(noiseIms,2)

    % take 2dfft ...
    %noiseImsFFT{1,ii}=fft2(noiseIms{1,ii});
    
    % run fftMaskFilter
    noiseImsFFT{1,ii}= fftMaskFilter(noiseIms{1,ii},fSpaceRad,filter);

    % For randomization of order later..
    noiseImsFFT{2,ii}=rand(1); % assign random number

end


%% Combine normal and fftMask filtered images into same array
noiseImsComb=vertcat(noiseIms,noiseImsFFT);


%% For each noise image, compute structural similarity between noise image and each base image
ssmat=cell(6,size(noiseImsComb,2));

for uu = 1:size(noiseImsComb,2)
    img = noiseImsComb{1,uu};
    [ssimvalL,ssimmapL]  = ssim(img,BI_L);
    [ssimvalR,ssimmapR]  = ssim(img,BI_R);
    ssimvalDiff=ssimvalL-ssimvalR;
    ssimvalDiffMap=ssimmapL-ssimmapR;
    ssmat{1,uu}=ssimvalL;
    ssmat{2,uu}=ssimvalR;
    ssmat{3,uu}=ssimmapL;
    ssmat{4,uu}=ssimmapR;
    ssmat{5,uu}=ssimvalDiff;
    ssmat{6,uu}=ssimvalDiffMap;

end

%% Get the max and min ssim vals for bothe the right and left base images
ssimvalL_min=min(cell2mat(ssmat(1,:)));
ssimvalL_max=max(cell2mat(ssmat(1,:)));
ssimvalR_min=min(cell2mat(ssmat(2,:)));
ssimvalR_max=max(cell2mat(ssmat(2,:)));

ssimvalDiff_min=min(cell2mat(ssmat(5,:)));
ssimvalDiff_max=max(cell2mat(ssmat(5,:)));

%% Use Linspace to set up a 100 value vector ranging from ssimvalDiff_min to ssimvalDiff_max
difValProbMat=zeros(2,100); % intitialize
difValProbMat(1,:)=linspace(ssimvalDiff_min,ssimvalDiff_max,100); % ssim difference
difValProbMat(2,:)=linspace(0,1,100); % corresponding response probability (Right vs. Left);

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
    pRight=difValProbMat(2,idx);
    pLeft=1-pRight;

    % feed these probabilites into randsample to select a random response
    resp = randsample(respOps,1,true,[pLeft, pRight]);

    % save response in the output respmat
    respMat(1,zz)=resp;

    % randomly decide whether base image is right v left.. (no skewing..
    % 50/50 prob..
    tcon = randsample(respOps,1,true,[0.5, 0.5]);
    respMat(2,zz)=tcon;
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

for hh = 1:size(noiseImsComb_ssmat,2)

    if respMat{hh,1}==0 && respMat{hh,2}==0
        CI_11{end+1}=ssmat{6,hh};       
    end

    if respMat{hh,1}==0 && respMat{hh,2}==1
        CI_12{end+1}=ssmat{6,hh};       
    end

    if respMat{hh,1}==1 && respMat{hh,2}==0
        CI_21{end+1}=ssmat{6,hh};
    end

    if respMat{hh,1}==1 && respMat{hh,2}==1
        CI_22{end+1}=ssmat{6,hh};
    end

end


% Normal Images
CI_11mat=zeros(size(CI_11{1,1},1),size(CI_11{1,1},2),length(CI_11));
for jj = 1:length(CI_11)
CI_11mat(:,:,jj)=CI_11{1,jj};
end
clear jj;


CI_21mat=zeros(size(CI_21{1,1},1),size(CI_21{1,1},2),length(CI_21));
for jj = 1:length(CI_21)
CI_21mat(:,:,jj)=CI_21{1,jj};
end
clear jj;

CI_12mat=zeros(size(CI_12{1,1},1),size(CI_12{1,1},2),length(CI_12));
for jj = 1:length(CI_12)
CI_12mat(:,:,jj)=CI_12{1,jj};
end
clear jj;

CI_22mat=zeros(size(CI_22{1,1},1),size(CI_22{1,1},2),length(CI_22));
for jj = 1:length(CI_22)
CI_22mat(:,:,jj)=CI_22{1,jj};
end
clear jj;




%% --- Compute CIs --- %
% Normal images
CI_results1 = (mean(CI_11mat,3) + mean(CI_21mat,3)) - (mean(CI_12mat,3) + mean(CI_22mat,3));
CI_results2 = (mean(CI_12mat,3) + mean(CI_22mat,3)) - (mean(CI_11mat,3) + mean(CI_21mat,3));

% % FFT filtered Images
% CI_resultsFFT_1 = (mean(CI_11FFTmat,3) + mean(CI_21FFTmat,3)) - (mean(CI_12FFTmat,3) + mean(CI_22FFTmat,3));
% CI_resultsFFT_2 = (mean(CI_12FFTmat,3) + mean(CI_22FFTmat,3)) - (mean(CI_11FFTmat,3) + mean(CI_21FFTmat,3));

figure;
imshow(CI_results1);
figure;
imshow(CI_results2);

