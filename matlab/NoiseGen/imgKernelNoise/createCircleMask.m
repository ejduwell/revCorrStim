function img = createCircleMask(imgWidth, imgHeight, circleRadius)
    % Create a grid of coordinates
    [x, y] = meshgrid(1:imgWidth, 1:imgHeight);
    
    % Compute the center of the image
    centerX = imgWidth / 2;
    centerY = imgHeight / 2;
    
    % Compute the squared distance to the center for each coordinate
    distSq = (x - centerX).^2 + (y - centerY).^2;
    
    % Generate the binary image: 1 for pixels within the circle's radius, 0 otherwise
    img = distSq <= circleRadius^2;
end
