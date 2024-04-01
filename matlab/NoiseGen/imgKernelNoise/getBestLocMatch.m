function getBestLocMatch(largeImage,sampleImage)

% Load the images
% largeImage = imread('path_to_large_image.jpg');
% sampleImage = imread('path_to_sample_image.jpg');

% Ensure the images are in grayscale
% if size(largeImage, 3) == 3
%     largeImage = rgb2gray(largeImage);
% end
% if size(sampleImage, 3) == 3
%     sampleImage = rgb2gray(sampleImage);
% end

% Perform normalized cross-correlation
correlationOutput = normxcorr2(sampleImage, largeImage);

% Create a mask for valid correlation regions
mask = true(size(correlationOutput));
mask(1:size(sampleImage,1)-1, :) = false;
mask(end-size(sampleImage,1)+2:end, :) = false;
mask(:, 1:size(sampleImage,2)-1) = false;
mask(:, end-size(sampleImage,2)+2:end) = false;

% Apply the mask to the correlation output
correlationOutput(~mask) = -inf;

% Find the peak in the masked correlation output
[yPeak, xPeak] = find(correlationOutput == max(correlationOutput(:)));

% Adjust for the offset
yOffset = yPeak - size(sampleImage, 1);
xOffset = xPeak - size(sampleImage, 2);

% Validate and adjust the indices to ensure they are within the bounds of largeImage
yStart = max(1, yOffset + 1);
yEnd = min(size(largeImage, 1), yOffset + size(sampleImage, 1));
xStart = max(1, xOffset + 1);
xEnd = min(size(largeImage, 2), xOffset + size(sampleImage, 2));

% Extract the matched region
matchedRegion = largeImage(yStart:yEnd, xStart:xEnd);

close all
figure;
% display original tile..
subplot(1,3,1);
imshow(uint8(sampleImage));
title('Noise Tile');

% Display the matched region
subplot(1,3,2);
imshow(uint8(matchedRegion));
title('Best Matched Region in BI');

% Draw a rectangle around the matched region on the larger image
subplot(1,3,3);
imshow(uint8(largeImage));
hold on;
rectangle('Position', [xOffset, yOffset, size(sampleImage, 2), size(sampleImage, 1)], ...
          'EdgeColor', 'r');
title('Best Matched Region Within Full BI');

% Release hold
hold off;

end