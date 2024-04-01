

%% Set Parameters
chkrBrdSz=25; % checkerboard size (checks.. this times 2 is the number of checks in x and y direction..) 
sizeArrayOut=5000; % size of array (x and y dimensions in pixels..)

% Output directory location
outDir="/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/matlab/BaseImgGen/texImgs"; 

% Autogen filename based on parameters set above..
outName=strcat("checkerBrd_",num2str(sizeArrayOut),"by",num2str(sizeArrayOut),"Array","_",num2str(chkrBrdSz*2),"by",num2str(chkrBrdSz*2),"_Checks",".png");

%% Run makeCheckerboard_v2 Command
checkerboard = makeCheckerboard_v2(chkrBrdSz,sizeArrayOut);

%% Save Image
startDir=pwd; % save start location
cd(outDir); % enter output directory
imwrite(checkerboard,outName); % save the checkerboard image..
cd(startDir); % return to starting position
