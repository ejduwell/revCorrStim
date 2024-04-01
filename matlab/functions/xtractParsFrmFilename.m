function [parValz, parTagz] = xtractParsFrmFilename(fname,parTagz,pType)
% This function extracts the value of parameters denoted by tags in
% filenames and returns the parameter value detected...
%
% INPUTS:
% fname: file name .. string
% parTag: string tag in filename denoting a parameter value used to produce
% the file contents.. this function assumes all tags follow the convention:
% 'TAG_PARVALUE' where the TAG is the parTag and PARVALUE is the numeric
% value of the parameter denoted by TAG..
%
% OUTPUTS:
% parVal: parameter value detected .. number (can be either integer or floating point)

% separate the filename away from its path and extension...
[~,name_only,~] = fileparts(fname); 

% Substring extraction:
parValz = zeros(1,length(parTagz)); %initialize output array
if pType=="str"
    parValz=string(parValz);
end
my_string = name_only;

for i=1:length(parTagz)
    parTag = parTagz(i);
    if pType=="num"
    pattern = strcat(parTag,'_([0-9]+(?:\.[0-9]+)?).*');
    matches = regexp(my_string, pattern, 'tokens');
    if ~isempty(matches)
        parVal = str2double(matches{1}{1});
    else
        parVal = str2double("");
    end
    end
    
    if pType=="str"

    pattern = strcat(parTag,'_([^_]+)_');
    matches = regexp(my_string, pattern, 'tokens');
    if ~isempty(matches)
        parVal = matches{1}{1};
    else
        parVal = "NA";
    end
    end
    parValz(i)=parVal; %add parVal to output parValz vector
end

end