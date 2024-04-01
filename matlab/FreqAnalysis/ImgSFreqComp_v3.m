function RMSE = ImgSFreqComp_v3(noiseImPars)
%noise = ["noiseData.NoiseType = 'white';","noiseData.mu = 0;","noiseData.std = 1;","noiseData.imgSize   = 512;"];   % Size of noise image (enter 2 values if x & y dims are different, ex: [256,512] )


%% Get base and noise images
[baseIm, noiseData] = noiseBaseImDescFile();

% --- Required Input1 (base image) --- %
% im1_in = "/Users/eduwell/OneDrive - mcw.edu/Documents/SNAP/projects/Rev_Corr/code/RevCorr/working_dir/out_test/revcorrBI_occ_0_ori_z_CR1CRpw1CRpl75CRal1CRob2CRdm3214_T1spheresTr11_5Tr20_5T_al1_L1l130l230LWT0.7L_al0_7.png";
im1 = imread(baseIm);
if length(size(im1))>2 % check if it has a color channel..
im1 = rgb2gray(im1); % convert to greyscale..
end
im1_size = size(im1);

% --- Required Input2 (noise image) --- %
% Unpack the noise vector describing the desired noise/execute strings
% inside using evalStrVect
%noiseData = evalStrVect(noiseImPars);
%noiseData.NoiseType = 'gabor';  % Noise type (white)
%noiseData.mu = 0;
%noiseData.std = 1;
%noiseData.imgSize   = 512;      % Size of noise image (enter 2 values if x & y dims are different, ex: [256,512] )

if noiseImPars.imgSize == "match"
%Make the size match the base image..
noiseImPars.imgSize   = [im1_size(1), im1_size(2)];      % Size of noise image (enter 2 values if x & y dims are different, ex: [256,512] )
end


% --- Call function --- %
[im2,~] = revCorrNoise_v2(noiseImPars); % myNoise returns noise field 


%% Generate Spatial Frequency Power Spectra for Images 1 and 2
% Image 1 (base image)
% #########################################################################
% make powerspectrum plot
f1 = figure;
[sf_im1, nps_im1] = powerspectrum(im1,[1,1],400);
sfData_im1 = [sf_im1, nps_im1];
semilogy(sf_im1, nps_im1);
%plot(sf_im1, log10(nps_im1));
ylim([0 1])
% plot original image too ..
figure 
imshow(im1)
%imagesc(im1)
%colormap("gray");
% #########################################################################

% Image 2 (noise image)
% #########################################################################
% make powerspectrum plot
f2 = figure;
[sf_im2, nps_im2] = powerspectrum(im2,[1,1],400);
sfData_im2 = [sf_im2 nps_im2];
semilogy(sf_im2, nps_im2);
ylim([0 1])
% plot original image too ..
figure
%imshow(im2)
imagesc(im2)
colormap("gray");
% #########################################################################

%% Compare the Power Spectra
% Calculate RMSE between the two powerspectra..
V1 = nps_im1;
V2 = nps_im2;
RMSE = sqrt(mean((V1-V2).^2));
disp("Mean Square Error Between im1 and im2 Spatial Frequency Power Spectra:")
disp(num2str(RMSE))

close all % close all the open figure windows.. so they don't stack up during optimization..

