function imzOut = getImzWithTagStr(imDir,imTag,adjSize,desiredSize)

%--------------------------------------------------------------------------
% Note: all output images will be grayscale..
%--------------------------------------------------------------------------

% Get set of all images file paths in directory containing the tag
%--------------------------------------------------------------------------

% Get a list of all files in the directory
files = dir(imDir);

% Initialize an empty cell array to store the filenames
imFiles = {};

% Iterate over each file in the directory
for i = 1:numel(files)
    % Check if the file ends with ".png"
    if endsWith(files(i).name, imTag, 'IgnoreCase', true)
        % Add the filename to the cell array
        imFiles{end+1} = strcat(files(i).folder,"/",files(i).name);
    end
end
%--------------------------------------------------------------------------

if adjSize == 0
    % If no size adjustment requested, read in sample and use input
    % dimensions
    % read in sample from top of stack to get the dimensions..
    sampleIm=imread(imFiles{1});
    smplDims=size(sampleIm);
    desiredSize=[smplDims(1),smplDims(2)];
elseif adjSize == 1
    % if adjustment requested,leave desiredSize as desiredSize..
end

% Read in the images that contained the tag
imzOut = uint8(zeros(desiredSize(1),desiredSize(2),length(imFiles))); % preallocate
for ii = 1:length(imFiles)
    im=imread(imFiles{ii});
    im=imgReformater(im,desiredSize);
    imzOut(:,:,ii)=im;
end
