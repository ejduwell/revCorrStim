function [matchedMaxCoords, matchedMinCoords, maxCor, minCor] = getBestLocMatch_v2(largeImage,sampleImage)

% First test to be sure all the values in the sampleImage ("template")
% are not the same.. this is a requirement of normxcorr2
% if this is not the case, add a very small value (0.0001) to a random location
% in the template to satisfy normxcorr2's quirky requirement.. 
%--------------------------------------------------------------------------
% Reshape the matrix into a column vector and find unique values
uniqueValues = unique(sampleImage(:));
% Count the number of unique values
nUniqueValz = length(uniqueValues);
uValTest=nUniqueValz>1;
switch uValTest
case 0
    yDtmp=size(sampleImage,1);
    xDtmp=size(sampleImage,2);
    ranXtmp=randi(xDtmp);
    ranYtmp=randi(yDtmp);
    sampleImage(ranYtmp,ranXtmp)=sampleImage(ranYtmp,ranXtmp)+0.0001;
end
%--------------------------------------------------------------------------

% Perform normalized cross-correlation
correlationOutput = normxcorr2(sampleImage, largeImage);

% Create a mask for valid correlation regions
mask = true(size(correlationOutput));
mask(1:size(sampleImage,1)-1, :) = false;
mask(end-size(sampleImage,1)+2:end, :) = false;
mask(:, 1:size(sampleImage,2)-1) = false;
mask(:, end-size(sampleImage,2)+2:end) = false;

% Apply the mask to the correlation output
%correlationOutput(~mask) = -inf;
correlationOutput(~mask) = 0;

% Find the peak in the masked correlation output
maxCor = max(correlationOutput(:));
[yPeak, xPeak] = find(correlationOutput == maxCor);

% Find the mimimum in the masked correlation output
minCor = min(correlationOutput(:));
[yMin, xMin] = find(correlationOutput == minCor);

% Adjust for the offset
yOffset = yPeak - size(sampleImage, 1);
xOffset = xPeak - size(sampleImage, 2);

yOffset2 = yMin - size(sampleImage, 1);
xOffset2 = xMin - size(sampleImage, 2);

% Validate and adjust the indices to ensure they are within the bounds of largeImage
yStart = max(1, yOffset + 1);
yEnd = min(size(largeImage, 1), yOffset + size(sampleImage, 1));
xStart = max(1, xOffset + 1);
xEnd = min(size(largeImage, 2), xOffset + size(sampleImage, 2));

yStart2 = max(1, yOffset2 + 1);
yEnd2 = min(size(largeImage, 1), yOffset2 + size(sampleImage, 1));
xStart2 = max(1, xOffset2 + 1);
xEnd2 = min(size(largeImage, 2), xOffset2 + size(sampleImage, 2));

% Extract the matched regions of max/min correlation..
%matchedRegionMax = largeImage(yStart:yEnd, xStart:xEnd);
%matchedRegionMin = largeImage(yStart2:yEnd2, xStart2:xEnd2);

% Save the Coordinates for the matched region
matchedMaxCoords=[yStart,yEnd,xStart,xEnd];
matchedMinCoords=[yStart2,yEnd2,xStart2,xEnd2];

end