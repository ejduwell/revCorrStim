function rxivMatlabCode(baseDir,outDirBase)
 % matlab wrapper for E.J. Duwell's rxivMatlabCode.sh bash function
 %
 % rxivMatlabCode.sh archives all matlab code (files ending in '.m') within
 % a specified directory for future reference. This is intended to make it
 % easy to save the precise state of the code as run on a particular day,
 % experiment, etc. for future reference.
 % 
 % Input parameters:
 % baseDir    : should be the full path to the program base directory in which your matlab code resides ...
 % outDirBase : should be the full path to the output parent/base directory where you want to save the archived code ...

% Get path info on this machine..
pathToThisFile=which("rxivMatlabCode.m");
currentDir = fileparts(pathToThisFile);

% ensure rxivMatlabCode.sh is executable..
chmodCmd=strcat("chmod +wrx ",currentDir,"/rxivMatlabCode_v2.sh ");
system(chmodCmd);

% build bash command to call rxivMatlabCode.sh
rxivCmd=strcat("bash ",currentDir,"/rxivMatlabCode_v2.sh ", "'",baseDir,"'"," ","'",outDirBase,"'");

% run the command in the bash shell..
system(rxivCmd);

end