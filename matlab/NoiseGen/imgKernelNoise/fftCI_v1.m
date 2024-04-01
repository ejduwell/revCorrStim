%% Parameters
matFile="189542_20231009T185227.mat";
rng("shuffle");


%% Read In Noise Images
data = load(matFile);

noiseIms=data.noiseImgMat;


%% Build mask for frequency space..
tmpImg="";

%% Take 2dFFT of of each
noiseImsFFT=cell(2,size(noiseIms,2));

for ii =1:size(noiseIms,2)

    % take 2dfft ...
    noiseImsFFT{1,ii}=fft2(noiseIms{1,ii});

    % For randomization of order later..
    noiseImsFFT{2,ii}=rand(1); % assign random number

end



%% Assign images to CI_11 --> CI_22 randomly..

noiseImsComb=vertcat(noiseIms,noiseImsFFT);

% sort to set random order..
noiseImsComb=noiseImsComb';
noiseImsComb=sortrows(noiseImsComb,3);
noiseImsComb=noiseImsComb';

nImgs=length(noiseImsComb);
nImgPerGrp=nImgs/4;

% Normal Images
CI_11 = noiseImsComb(1,1:nImgPerGrp);
CI_11mat=zeros(size(CI_11{1,1},1),size(CI_11{1,1},2),length(CI_11));
for jj = 1:length(CI_11)
CI_11mat(:,:,jj)=CI_11{1,jj};
end
clear jj;

CI_21 = noiseImsComb(1,((1*nImgPerGrp)+1):(((1*nImgPerGrp))+nImgPerGrp));
CI_21mat=zeros(size(CI_11{1,1},1),size(CI_11{1,1},2),length(CI_11));
for jj = 1:length(CI_21)
CI_21mat(:,:,jj)=CI_21{1,jj};
end
clear jj;

CI_12 = noiseImsComb(1,((2*nImgPerGrp)+1):(((2*nImgPerGrp))+nImgPerGrp));
CI_12mat=zeros(size(CI_12{1,1},1),size(CI_12{1,1},2),length(CI_12));
for jj = 1:length(CI_12)
CI_12mat(:,:,jj)=CI_12{1,jj};
end
clear jj;

CI_22 = noiseImsComb(1,((3*nImgPerGrp)+1):(((3*nImgPerGrp))+nImgPerGrp));
CI_22mat=zeros(size(CI_22{1,1},1),size(CI_22{1,1},2),length(CI_22));
for jj = 1:length(CI_22)
CI_22mat(:,:,jj)=CI_22{1,jj};
end
clear jj;


% fft of Images
CI_11FFT = noiseImsComb(2,1:nImgPerGrp);
CI_11FFTmat=zeros(size(CI_11FFT{1,1},1),size(CI_11FFT{1,1},2),length(CI_11FFT));
for jj = 1:length(CI_11FFT)
CI_11FFTmat(:,:,jj)=CI_11FFT{1,jj};
end
clear jj;

CI_21FFT = noiseImsComb(2,((1*nImgPerGrp)+1):(((1*nImgPerGrp))+nImgPerGrp));
CI_21FFTmat=zeros(size(CI_21FFT{1,1},1),size(CI_21FFT{1,1},2),length(CI_21FFT));
for jj = 1:length(CI_11FFT)
CI_21FFTmat(:,:,jj)=CI_21FFT{1,jj};
end
clear jj;

CI_12FFT = noiseImsComb(2,((2*nImgPerGrp)+1):(((2*nImgPerGrp))+nImgPerGrp));
CI_12FFTmat=zeros(size(CI_12FFT{1,1},1),size(CI_12FFT{1,1},2),length(CI_12FFT));
for jj = 1:length(CI_12FFT)
CI_12FFTmat(:,:,jj)=CI_12FFT{1,jj};
end
clear jj;

CI_22FFT = noiseImsComb(2,((3*nImgPerGrp)+1):(((3*nImgPerGrp))+nImgPerGrp));
CI_22FFTmat=zeros(size(CI_22FFT{1,1},1),size(CI_22FFT{1,1},2),length(CI_22FFT));
for jj = 1:length(CI_22FFT)
CI_22FFTmat(:,:,jj)=CI_22FFT{1,jj};
end
clear jj;


%% --- Compute CIs --- %
% Normal images
CI_results1 = (mean(CI_11mat,3) + mean(CI_21mat,3)) - (mean(CI_12mat,3) + mean(CI_22mat,3));
CI_results2 = (mean(CI_12mat,3) + mean(CI_22mat,3)) - (mean(CI_11mat,3) + mean(CI_21mat,3));
% FFT of Images
CI_resultsFFT_1 = (mean(CI_11FFTmat,3) + mean(CI_21FFTmat,3)) - (mean(CI_12FFTmat,3) + mean(CI_22FFTmat,3));
CI_resultsFFT_2 = (mean(CI_12FFTmat,3) + mean(CI_22FFTmat,3)) - (mean(CI_11FFTmat,3) + mean(CI_21FFTmat,3));

CI_resultsFFT_1INV=ifft2(CI_resultsFFT_1);
CI_resultsFFT_2INV=ifft2(CI_resultsFFT_2);
