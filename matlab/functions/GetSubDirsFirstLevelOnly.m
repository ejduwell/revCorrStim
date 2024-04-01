function subDs = GetSubDirsFirstLevelOnly(folder)

% save start location..
startDir=pwd;

% enter main dir of interest
cd(folder); 

% get subdir names in this folder..
d=dir(folder);

isub = [d(:).isdir]; % returns logical vector
subDs = {d(isub).name}';
subDs(ismember(subDs,{'.','..'})) = [];

% return to start location..
cd(startDir);

end