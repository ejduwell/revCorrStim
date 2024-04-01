%% Parameters for running CI sim..

%outDir='/home/eduwell/Matlab_Sandbox/imgKernelNoise/exp2Test/individualFaceBIs/5000';
outDir='/home/eduwell/Matlab_Sandbox/imgKernelNoise/CIquality_expmt_v2/aveNeutral-Smiling_BIs/20000';
outName="structOut_face_20000.mat";
% Note: Other specific pars are set in KrnlSmplNoise_v10_Demo_multi.m
% Be sure to set these to the desired values..

%% Run CI sim

% save start time
strtTime=datetime;

structOut = KrnlSmplNoise_v10_Demo_multi();

% save end time / compute duration..
endTime=datetime;
runTimeDur=endTime-strtTime;
runTimeSecs=seconds(runTimeDur);
timeString=strcat("Simulating CIs took: ",num2str(runTimeSecs)," seconds...");

disp(" ");
disp(" ");
disp("===================================================================");
disp(timeString);
disp("===================================================================");
disp(" ");
disp(" ");

% Save Data
strtDir=pwd; % save start position
cd(outDir); % go to output dir..
save(outName,"structOut","-v7.3"); % save struct..
cd(strtDir); % return to start location

% Clear the structOut variable as its probably really big / we just saved..
clear structOut

%% Load in struct if requested
%path2structFile='/home/eduwell/Matlab_Sandbox/imgKernelNoise/CIquality_expmt_v1/face/cropped/1000/combAveFaceKnoise_smile_v_neutral_aveBIs_1000.mat';
%path2structFile='/home/eduwell/Matlab_Sandbox/imgKernelNoise/CIquality_expmt_v1/gabor/take2/structOut_1000_BI_6hz_v2.mat';
%path2structFile='/home/eduwell/Matlab_Sandbox/imgKernelNoise/exp2Test/aveNeutral-Smiling_BIs/1000/structOut.mat';
%path2structFile='/home/eduwell/Matlab_Sandbox/imgKernelNoise/exp2Test/aveNeutral-Smiling_BIs/5000/structOut_face_5000.mat';
%path2structFile='/home/eduwell/Matlab_Sandbox/imgKernelNoise/exp2Test/aveNeutral-Smiling_BIs/10000/structOut_face_10000.mat';
%path2structFile='/home/eduwell/Matlab_Sandbox/imgKernelNoise/exp2Test/aveNeutral-Smiling_BIs/20000/structOut_face_20000.mat';

path2structFile=strcat(outDir,"/",outName);
s = load(path2structFile); % load in the struct..

%% Set Parameters for Overlay..

path2BI="/home/eduwell/Matlab_Sandbox/imgKernelNoise/pyAlignFaceOutTest6/combined_average_face.jpg"; %path to base image file
%path2BI="/home/eduwell/Matlab_Sandbox/imgKernelNoise/baseImgs/BaseImgSF_6L.png"; %path to base image file

desiredDimz=[256,256]; % bi/ci dims
biIn = imgReformater(imread(path2BI),desiredDimz); % read in and format base image to desired dimensions
biStr="All Cond Ave Face"; %string label indicating what the base image is in titles..
%biStr="BI 1 (left)"; %string label indicating what the base image is in titles..


% noise type/configuration used in sim
whichConfig="config1"; %specify noise configuration in struct you want to visualize
whichNtype="white"; % either kernel noise: "krnlNz" or white noise: "white"
% Specify which CI from this configuration you want to use
whichCI="rawCI2";% either rawCI1, rawCI2, smthCI1, or smthCI2
% pull CI matching specifications above out of the struct/set as ciIn..
ciIn=s.structOut.(whichConfig).CIs.(whichNtype).(whichCI); 

% specify ideal CI to plot..
whichIdealCI="idealCI2"; % specify which ideal CI version you want (either idealCI1 or idealCI2).. should probably match the ci version you choose for whichCI

idealCI=s.structOut.idealCIs.(whichIdealCI); % pull selection from the struct/set as idealCI
%idealCI = s.structOut.(whichConfig).CIs.(whichNtype).idealCIs.(whichIdealCI); % pull selection from the struct/set as idealCI

% specify CI overlay color parameters
clrMap="parula"; % colormap used on colorized/thresholded CI overlays for real and ideal CIs over the base image..
alphVal=0.2; % alpha/transparency value used on thresholded/colorized CI overlay..

% specify weights and thresholds
ciWt=0.5; % weight given to ci in ci/bi overlay
ciThr=90; % percent threshold used for thesholded/colorized CI and ideal CI overlays..

% specify figure dimensions..
figDims= [2000,1000]; % figure dimensions ([xdim,ydim])

% auto-set noiseConfigStr describing noise used in this configuration based
% on whichNtype specified above..
if whichNtype == "krnlNz"
noiseConfigStr=strcat("Kernel Noise Built With: ",mat2str(s.structOut.(whichConfig).krnlNz_krnlsUsed));
elseif whichNtype == "white"
noiseConfigStr="White Noise Control";
else
noiseConfigStr="";
end

%% Call plotCIasOverlay
plotCIasOverlay(biIn,ciIn,ciWt,ciThr, figDims, clrMap,alphVal,idealCI,biStr,whichCI,whichIdealCI,whichConfig,noiseConfigStr)