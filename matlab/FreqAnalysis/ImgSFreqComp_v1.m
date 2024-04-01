tic

% Read in a standard MATLAB gray scale demo image.
% folder = fileparts(which('cameraman.tif')); % Determine where demo folder is (works with all versions).
% baseFileName = 'cameraman.tif';
% % Get the full filename, with path prepended.
% fullFileName = fullfile(folder, baseFileName);
% im1 = imread(fullFileName);


% --- Required Input1 --- %
noiseData.NoiseType = 'grating';  % Noise type (white)
noiseData.mu = 0;
noiseData.std = 1;
noiseData.imgSize   = 512;      % Size of noise image (enter 2 values if x & y dims are different, ex: [256,512] )

% --- Call function --- %
[im1,~] = revCorrNoise(noiseData); % myNoise returns noise field 

% --- Required Input2 --- %
noiseData.NoiseType = 'grating';  % Noise type (white)
noiseData.mu = 0;
noiseData.std = 1;
noiseData.imgSize   = 512;      % Size of noise image (enter 2 values if x & y dims are different, ex: [256,512] )
% --- Call function --- %
[im2,~] = revCorrNoise(noiseData); % myNoise returns noise field 

% Image 1
% im1_in = "/Users/eduwell/OneDrive - mcw.edu/Documents/SNAP/projects/Rev_Corr/code/RevCorr/working_dir/out_test/revcorrBI_occ_0_ori_z_CR1CRpw1CRpl75CRal1CRob2CRdm3214_T1spheresTr11_5Tr20_5T_al1_L1l130l230LWT0.7L_al0_7.png";
% im1 = imread(im1_in);
% im1 = rgb2gray(im1); % convert to greyscale..
[sf_im1 nps_im1] = powerspectrum(im1,[1,1],400);
sfData_im1 = [sf_im1 nps_im1];
f1 = figure;
semilogy(sf_im1, nps_im1);
%plot(sf_im1, log10(nps_im1));
ylim([0 1])

figure 
%imshow(im1)
imagesc(im1)
colormap("gray");

% Image 2
% im2_in = "/Users/eduwell/OneDrive - mcw.edu/Documents/SNAP/projects/Rev_Corr/code/RevCorr/working_dir/out_test/revcorrBI_occ_0_ori_m_CR0CRpw1CRpl75CRal1CRob2CRdm3214_T1spheresTr11_5Tr20_5T_al1_L1l1120l2120LWT0.7L_al0_7.png";
% im2 = imread(im2_in);
% im2 = rgb2gray(im2); % convert to greyscale..
f2 = figure;
[sf_im2, nps_im2] = powerspectrum(im2,[1,1],400);
sfData_im2 = [sf_im2 nps_im2];
semilogy(sf_im2, nps_im2);

ylim([0 1])

figure
%imshow(im2)
imagesc(im2)
colormap("gray");

% Calculate RMSE between the two powerspectra..
V1 =nps_im1;
V2 = nps_im2;
RMSE = sqrt(mean((V1-V2).^2));
disp("Mean Square Error Between im1 and im2 Spatial Frequency Power Spectra:")
disp(num2str(RMSE))

toc








% figure
% title("input image")
% imshow(im1)
% 
% figure
% % cmap = colormap("jet");
% % title("fft2 of input image")
% title(" 'powerspectrum' function outpot of input image..")
% % imshow(im1_fft);
% imshow(im1_ps)

%% mean 1d radial profile of 2d fft approach
% from: https://www.mathworks.com/matlabcentral/answers/340538-how-to-average-the-2d-spectrum-of-an-image-from-fft2-to-get-1d-spectrum
% Ethan's general understanding.. : 
% ########################################################################
% General goal: we want to plot a 1d "powerspectrum" of the mean magnitude
% of spatial frequencies across a specified range within a 2d image
%
% Starting image is in the position domain/is 2d..
% 1) take the 2d fft of the input image..
%     this should result in a 2d plot of frequency amplitude values..
%     the center of the plot is "low frequency"
%     higher eccentricity positions/the edges represent "higher spatial
%     frequency" amplitudes. Values along a given "radius" are therefore
%     amplitudes corresponding to to the same spatial frequency..
% 2) Because there are multiple amplitude values along a given radius, and
%    because we want out powerspectrum output to be "1D", we need to
%    collapse/combine this in some way.. In this method, I believe they take
%    the mean/average.
% 3) repeat the process above on a series of radii to make a 1d power
% spectrum..


% Things we'll need to figure out..
% How to map a given "radius" to a particular spatial frequency value.
% look in help page for fft2 can we specify output units?
% ########################################################################

% DEMO:
% ########################################################################
% run_demo2 = 0;
% if run_demo2 == 1
% %   2D FFT Demo to get average radial profile of the spectrum of an image.
% clc;    % Clear the command window.
% close all;  % Close all figures (except those of imtool.)
% imtool close all;  % Close all imtool figures.
% clear;  % Erase all existing variables.
% workspace;  % Make sure the workspace panel is showing.
% format longg;
% format compact;
% fontSize = 20;
% 
% % Change the current folder to the folder of this m-file.
% if(~isdeployed)
%   cd(fileparts(which(mfilename)));
% end
%   
% % Check that user has the Image Processing Toolbox installed.
% hasIPT = license('test', 'image_toolbox');
% if ~hasIPT
%   % User does not have the toolbox installed.
%   message = sprintf('Sorry, but you do not seem to have the Image Processing Toolbox.\nDo you want to try to continue anyway?');
%   reply = questdlg(message, 'Toolbox missing', 'Yes', 'No', 'Yes');
%   if strcmpi(reply, 'No')
%     % User said No, so exit.
%     return;
%   end
% end
% 
% % Read in a standard MATLAB gray scale demo image.
% folder = fileparts(which('cameraman.tif')); % Determine where demo folder is (works with all versions).
% baseFileName = 'cameraman.tif';
% % Get the full filename, with path prepended.
% fullFileName = fullfile(folder, baseFileName);
% % Check if file exists.
% if ~exist(fullFileName, 'file')
%   % File doesn't exist -- didn't find it there.  Check the search path for it.
%   fullFileName = baseFileName; % No path this time.
%   if ~exist(fullFileName, 'file')
%     % Still didn't find it.  Alert user.
%     errorMessage = sprintf('Error: %s does not exist in the search path folders.', fullFileName);
%     uiwait(warndlg(errorMessage));
%     return;
%   end
% end
% % Read in image.
% grayImage = imread('cameraman.tif');
% [rows, columns, numberOfColorChannels] = size(grayImage)
% if numberOfColorChannels > 1
%   grayImage = rgb2gray(grayImage);
% end
% 
% % Display original grayscale image.
% subplot(2, 3, 1);
% imshow(grayImage)
% axis on;
% title('Original Gray Scale Image', 'FontSize', fontSize)
% % Enlarge figure to full screen.
% set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% 
% % Perform 2D FFTs
% fftOriginal = fft2(double(grayImage));
% % Move center from (1,1) to (129, 129) (the middle of the matrix).
% shiftedFFT = fftshift(fftOriginal);
% subplot(2, 3, 2);
% scaledFFTr = 255 * mat2gray(real(shiftedFFT));
% imshow(log(scaledFFTr), []);
% title('Log of Real Part of Spectrum', 'FontSize', fontSize)
% subplot(2, 3, 3);
% scaledFFTi = mat2gray(imag(shiftedFFT));
% imshow(log(scaledFFTi), []);
% axis on;
% title('Log of Imaginary Part of Spectrum', 'FontSize', fontSize)
% 
% % Display magnitude and phase of 2D FFTs
% subplot(2, 3, 4);
% shiftedFFTMagnitude = abs(shiftedFFT);
% imshow(log(abs(shiftedFFTMagnitude)),[]);
% axis on;
% colormap gray
% title('Log Magnitude of Spectrum', 'FontSize', fontSize)
% 
% % Get the average radial profile
% midRow = rows/2+1
% midCol = columns/2+1
% maxRadius = ceil(sqrt(129^2 + 129^2))
% radialProfile = zeros(maxRadius, 1);
% count = zeros(maxRadius, 1);
% for col = 1 : columns
%   for row = 1 : rows
%     radius = sqrt((row - midRow) ^ 2 + (col - midCol) ^ 2);
%     thisIndex = ceil(radius) + 1;
%     radialProfile(thisIndex) = radialProfile(thisIndex) + shiftedFFTMagnitude(row, col);
%     count(thisIndex) = count(thisIndex) + 1;
%   end
% end
% % Get average
% radialProfile = radialProfile ./ count;
% subplot(2, 3, 5:6);
% plot(radialProfile, 'b-', 'LineWidth', 2);
% grid on;
% title('Average Radial Profile of Spectrum', 'FontSize', fontSize)
% end
% % ########################################################################