
%% Create the noise images
% --- Required Input1 --- %
% im1_in = "/Users/eduwell/OneDrive - mcw.edu/Documents/SNAP/projects/Rev_Corr/code/RevCorr/working_dir/out_test/revcorrBI_occ_0_ori_z_CR1CRpw1CRpl75CRal1CRob2CRdm3214_T1spheresTr11_5Tr20_5T_al1_L1l130l230LWT0.7L_al0_7.png";
% im1 = imread(im1_in);
% im1 = rgb2gray(im1); % convert to greyscale..
noiseData.NoiseType = 'white';  % Noise type (white)
noiseData.mu = 0;
noiseData.std = 1;
noiseData.imgSize   = 512;      % Size of noise image (enter 2 values if x & y dims are different, ex: [256,512] )

% --- Call function --- %
[im1,~] = revCorrNoise(noiseData); % myNoise returns noise field 

% --- Required Input2 --- %
% im2_in = "/Users/eduwell/OneDrive - mcw.edu/Documents/SNAP/projects/Rev_Corr/code/RevCorr/working_dir/out_test/revcorrBI_occ_0_ori_m_CR0CRpw1CRpl75CRal1CRob2CRdm3214_T1spheresTr11_5Tr20_5T_al1_L1l1120l2120LWT0.7L_al0_7.png";
% im2 = imread(im2_in);
% im2 = rgb2gray(im2); % convert to greyscale..
noiseData.NoiseType = 'white';  % Noise type (white)
noiseData.mu = 0;
noiseData.std = 1;
noiseData.imgSize   = 512;      % Size of noise image (enter 2 values if x & y dims are different, ex: [256,512] )
% --- Call function --- %
[im2,~] = revCorrNoise(noiseData); % myNoise returns noise field 


%% Generate Spatial Frequency Power Spectra for Images 1 and 2
% Image 1
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
%imshow(im1)
imagesc(im1)
colormap("gray");
% #########################################################################

% Image 2
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

clear