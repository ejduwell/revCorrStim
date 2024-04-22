function whiteNoise2File_v1()

%% Parameters

Repz=20000; % sets number of image reps/trials used for each pass..
jobChunks=10; % sets the number of chunks you want to break the job 
              % (total reps) into..

nReps=Repz/jobChunks; % number of noise image copies/reps you want to run per pass..
nWorkers=8; % number of cpu nodes used in the parallelized portions..

%Specify desired size?
selectSize=1; % if 1, this means that we will use/resize the input image the selected size below
desiredSize=[512,512]; % (MUST BE SQUARE IMAGE AND DIMZ MUST BE POWER OF 2)

% Out Directory Base (parent directory where you want your output dirs
% created)
outDirMain="/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/noise";
baseOutDirName="512by512_whiteNoise_20000frms_smpl3";


%% Call Function

rng('shuffle');% make sure random number generator is shuffled

% Go to output parent directory...
strtDir=pwd; % save start location..
cd(outDirMain);

% Set up local parallel computing stuff..
%--------------------------------------------------------------------------
parpool("Processes",nWorkers);
%--------------------------------------------------------------------------

tic
frmCountr=1; % initialize frame countr
for qq=1:jobChunks

% Generate the noise images
WhtNzImz=uint8(zeros(desiredSize(1),desiredSize(2),nReps));
parfor kk = 1:nReps
whiteNoiseImg=uint8(randi(255, desiredSize));
WhtNzImz(:,:,kk)=whiteNoiseImg;
end

OutDirTmp=baseOutDirName;
if ~exist(OutDirTmp, 'dir')
    mkdir(OutDirTmp) % make the output directory..
end
cd(OutDirTmp); % enter it

% write out the noise files
nImzTemp=size(WhtNzImz,3);
for ii=1:nImzTemp
    imNameStr=strcat("noiseSample",num2str(frmCountr,'%05.f'),".png");
    frmCountr=frmCountr+1; % update frame countr for next pass..
    imwrite(WhtNzImz(:,:,ii),imNameStr);
end
cd .. %go back up a level..
end

% Close down local parpool
delete(gcp('nocreate'));

% return to start location
cd(strtDir);

disp(" ")
disp("Making White Noise Imz Took:")
toc

end