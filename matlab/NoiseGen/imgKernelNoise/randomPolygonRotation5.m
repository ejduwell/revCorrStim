function imgOut = randomPolygonRotation5(imgIn, nPolygons, nSides, sizeIn)
    % Get the dimensions of the input image
    [rows, cols, channels] = size(imgIn);
    
    % Initialize the output image
    imgOut = imgIn;
    
    % Create a mask for the polygon
    theta = linspace(0, 2*pi, nSides+1);
    xVertices = sizeIn * cos(theta);
    yVertices = sizeIn * sin(theta);
    
    % Calculate the bounding box of the polygon
    minX = floor(min(xVertices));
    maxX = ceil(max(xVertices));
    minY = floor(min(yVertices));
    maxY = ceil(max(yVertices));
    
    % Create a binary mask for the polygon
    [X, Y] = meshgrid(minX:maxX, minY:maxY);
    polygonMask = inpolygon(X, Y, xVertices, yVertices);
    
    % Initialize a mask to keep track of occupied regions
    occupiedMask = false(rows, cols);
    
    % Compute the full set of possible valid locations
    [centerXGrid, centerYGrid] = meshgrid((sizeIn + 1):(cols - sizeIn - 1), (sizeIn + 1):(rows - sizeIn - 1));
    validLocations = [centerXGrid(:), centerYGrid(:)];
    
    % Filter valid locations initially
    isValid = true(size(validLocations, 1), 1);
    for i = 1:size(validLocations, 1)
        centerX = validLocations(i, 1);
        centerY = validLocations(i, 2);
        for x = minX:maxX
            for y = minY:maxY
                if polygonMask(y - minY + 1, x - minX + 1)
                    xi = centerX + x;
                    yi = centerY + y;
                    if xi < 1 || xi > cols || yi < 1 || yi > rows || occupiedMask(yi, xi)
                        isValid(i) = false;
                        break;
                    end
                end
            end
            if ~isValid(i)
                break;
            end
        end
    end
    validLocations = validLocations(isValid, :);
    
    for i = 1:nPolygons
        % If no valid locations are found, exit the loop
        if isempty(validLocations)
            warning('Could not place all polygons without overlap.');
            break;
        end
        
        % Randomly select a valid location
        selectedIndex = randi(size(validLocations, 1));
        centerX = validLocations(selectedIndex, 1);
        centerY = validLocations(selectedIndex, 2);
        
        % Extract the polygon region from the image
        polygonRegion = zeros(size(polygonMask, 1), size(polygonMask, 2), channels, 'like', imgIn);
        for c = 1:channels
            for x = minX:maxX
                for y = minY:maxY
                    if polygonMask(y - minY + 1, x - minX + 1)
                        xi = centerX + x;
                        yi = centerY + y;
                        if xi >= 1 && xi <= cols && yi >= 1 && yi <= rows
                            polygonRegion(y - minY + 1, x - minX + 1, c) = imgIn(yi, xi, c);
                        end
                    end
                end
            end
        end
        
        % Apply random rotation to the polygon region
        angle = randi([1, nSides - 1]) * 360 / nSides;
        rotatedRegion = imrotate(polygonRegion, angle, 'crop');
        rotatedMask = imrotate(polygonMask, angle, 'crop');
        
        % Normalize the rotated mask
        rotatedMask(rotatedMask > 0) = 1;
        
        % Place the rotated polygon region back in the image
        for c = 1:channels
            for x = minX:maxX
                for y = minY:maxY
                    if rotatedMask(y - minY + 1, x - minX + 1)
                        xi = centerX + x;
                        yi = centerY + y;
                        if xi >= 1 && xi <= cols && yi >= 1 && yi <= rows
                            imgOut(yi, xi, c) = rotatedRegion(y - minY + 1, x - minX + 1, c);
                            occupiedMask(yi, xi) = true;
                        end
                    end
                end
            end
        end
        
        % Update the valid locations by removing the ones that are too close to the placed polygon
        buffer = sizeIn; % Define a buffer zone around the polygon
        isValid = (sqrt((validLocations(:, 1) - centerX).^2 + (validLocations(:, 2) - centerY).^2)/2) > buffer;
        validLocations = validLocations(isValid, :);
    end
end
