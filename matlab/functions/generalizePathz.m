
function outputArray = generalizePathz(main_dir,inputArray)
% This function adjusts a cell array of path strings (inputArray) to chop 
% the input "main_dir" off the beginning of all path strings and outputs an
% array of adjusted strings.. This was written for the intended use of
% adjusting paths which may have been saved from different computers but
% within a common program directory structure such that they can be
% compared without the initial machine-specific path info prior to the
% program's "main directory" (ie main_dir).

outputArray=cell(size(inputArray,1),1); % preallocate

for ii=1:size(inputArray,1)
strTemp=inputArray{ii,1};
%[~,matchTmp] = strsplit(strTemp,strcat('\s*',main_dir,'\s*'),'DelimiterType','RegularExpression');
[BefAfter,matchTmp] = strsplit(strTemp,strcat(main_dir,'\s*'),'DelimiterType','RegularExpression');
outputArray{ii,1}=BefAfter{1,2};
end