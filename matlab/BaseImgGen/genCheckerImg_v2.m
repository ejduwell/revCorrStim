

%% Set Parameters
chkrBrdSz=25; % checkerboard size (checks.. this times 2 is the number of checks in x and y direction..) 
sizeArrayOut=5000; % size of array (x and y dimensions in pixels..)

meanLum=0.3; % desired mean luminance of output checkerboard as fraction of total range 0-255..
chckrAmplitude=0.1; % checker amplitude (number from 0 to 1.. 1 means amplitde of checker contrast spans 100% the possible range of black to white, 0.5 is 50%, 0 is 0% )

% Output directory location
outDir="/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/matlab/BaseImgGen/texImgs"; 

%% Autogen filename based on parameters set above..
outName=strcat("checkerBrd_",num2str(sizeArrayOut),"by",num2str(sizeArrayOut),"Array","_",num2str(chkrBrdSz*2),"by",num2str(chkrBrdSz*2),"_Checks_meanLum_",num2str(meanLum),"_chkrAmp_",num2str(chckrAmplitude),".png");

%% 1. Run makeCheckerboard_v2 Command to make black and white checkerboard
checkerboard = makeCheckerboard_v2(chkrBrdSz,sizeArrayOut);

%% 2. Adjust the checker amplitude may taking weighted mean of checkerboard and  and a mono-lum array of equal mean lum.

% compute mean-lum of checkerboard
chkrMeanLum=mean(mean(checkerboard));

% make monoluminance/gray array with luminance equal to the mean computed above
monoLumArray=ones(size(checkerboard,1),size(checkerboard,2)); % initialize as an array of ones
monoLumArray=chkrMeanLum.*(monoLumArray); % multiply by chkrMeanLum

% Take weighted average between monoLumArray and the checkerboard using
% chckrAmplitude to set the relative weighting
mlarrayWt=1-chckrAmplitude; % compute weight for mono-lum array based on chkrMeanLum
%take weighted average of lum-only and tex-only images...
chkrMLarray_wAve = (mlarrayWt*double(monoLumArray) + chckrAmplitude*double(checkerboard));


%% 3. Adjust mean luminance to specified value

chkrMax=max(max(double(checkerboard)));
chkrMin=min(min(double(checkerboard)));
meanLumAbs=meanLum*(chkrMax-chkrMin);


% compute diff between mean chkrMLarray_wAve lum and desired output mean
% this will be the required adjustment to add to all pixels to push mean
% checker lum to desired value
chkrMLarray_wAve_mean=mean(mean(chkrMLarray_wAve));
meanAdjust=meanLumAbs-chkrMLarray_wAve_mean;

% create a matrix of values all equal to meanAdjust
adjustMat=ones(size(checkerboard,1),size(checkerboard,2)); % initialize as an array of ones
adjustMat=adjustMat*meanAdjust;

% sum adjustMat and  chkrMLarray_wAve to adjust to final desired output
% mean luminance
chkrOut=uint8(double(adjustMat)+double(chkrMLarray_wAve));

%% Save Image
startDir=pwd; % save start location
cd(outDir); % enter output directory
imwrite(chkrOut,outName); % save the checkerboard image..
cd(startDir); % return to starting position
