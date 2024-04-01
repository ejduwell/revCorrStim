%% Set parameters
%imPath="/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/noise/texOnly4BI_krnlNz_imzPerKrnl_111111100/noiseSample00001.png";
%imPath="/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/noise/lumOnlyBI_20000frms_krnlNz_imzPerKrnl_111111100/noiseSample00051.png";
%imPath="/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/images/test5_TexOnly-12-Mar-2024-18-48-17/occ/L/BaseIm_occ_0_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3315_T_1checkrz_Tr1_1.3_Tr2_0.76923_Tal_1_L_0_lum1_127_lum2_127_LWT_1_Lal_1_Ocon_1.png";
imPath="/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/noise/whiteNoise512-512/noiseSample_1.jpg";
figSizeX =1500*2;
figSizeY=750*2;
flipRate=0.2;

% lowFreqz = [0,0,0,0]; % Low frequency cutoff
% highFreqz = [10,100,500,1000]; % High frequency cutoff

% traveling band attempt
% bandWidSep=25;
% lowerLim=0;
% upperLim=256;
% stepz=100;
% lowFreqz = linspace(lowerLim,upperLim-bandWidSep,stepz); % Low frequency cutoffs
% highFreqz = linspace(lowerLim+bandWidSep,upperLim,stepz); % High frequency cutoffs

% widen from bottom attempt
lowerLim=0;
upperLim=256;
stepz=1000;
lowFreqz = zeros(1,stepz);
highFreqz = linspace(lowerLim,upperLim,stepz); % High frequency cutoffs

% widen from top attempt
% lowerLim=0;
% upperLim=256;
% stepz=1;
% lowFreqz = linspace(upperLim-1,lowerLim,stepz); 
% highFreqz = (ones(1,stepz)).*upperLim; % High frequency cutoffs

%% Code
inputImage = imread(imPath); % Replace 'yourImage.jpg' with your image file
nFiltImgs=min([length(lowFreqz),length(highFreqz)]);
figure('Position', [10 10 figSizeX figSizeY]); % open figure window
hold on;
for ii = 1:nFiltImgs
lowFreq = lowFreqz(ii); % Low frequency cutoff
highFreq = highFreqz(ii); % High frequency cutoff

filteredImage = applyBandPassFilter(inputImage, @myBandPassFilter,lowFreq,highFreq);


subplot(1,2,1);
imshow(inputImage);
title("original")
subplot(1,2,2);
imshow(filteredImage);
titleStr=strcat("bandpass filtered between highFreq=",num2str(highFreq)," and lowFreq=",num2str(lowFreq));
title(titleStr);
disp(" ");
disp(titleStr);
pause(flipRate);
end
hold off;
disp(" ");
input("Press Enter to Close..");
close all;
clear;