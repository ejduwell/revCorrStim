function imgOut = randomPolygonRotation2(imgIn, nPolygons, nSides, sizeIn)
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
    
    for i = 1:nPolygons
        isValidPosition = false;
        attempt = 0;
        maxAttempts = 10000;
        
        while ~isValidPosition && attempt < maxAttempts
            % Generate random center for the polygon
            centerX = randi([sizeIn + 1, cols - sizeIn - 1]);
            centerY = randi([sizeIn + 1, rows - sizeIn - 1]);
            
            % Check if the new polygon overlaps with existing ones
            isValidPosition = true;
            for x = minX:maxX
                for y = minY:maxY
                    if polygonMask(y - minY + 1, x - minX + 1)
                        xi = centerX + x;
                        yi = centerY + y;
                        if xi >= 1 && xi <= cols && yi >= 1 && yi <= rows
                            if occupiedMask(yi, xi)
                                isValidPosition = false;
                                break;
                            end
                        else
                            isValidPosition = false;
                            break;
                        end
                    end
                end
                if ~isValidPosition
                    break;
                end
            end
            attempt = attempt + 1;
        end
        
        if ~isValidPosition
            warning('Could not place all polygons without overlap.');
            break;
        end
        
        % Extract the polygon region from the image
        polygonRegion = zeros(size(polygonMask, 1), size(polygonMask, 2), channels, 'like', imgIn);
        for c = 1:channels
            for x = minX:maxX
                for y = minY:maxY
                    if polygonMask(y - minY + 1, x - minX + 1)
                        xi = centerX + x;
                        yi = centerY + y;
                        if xi >= 1 && xi <= cols && yi >= 1 && yi <= rows
                            polygonRegion(y - minY + 1, x - minX + 1, c) = ...
                                imgIn(yi, xi, c);
                        end
                    end
                end
            end
        end
        
        % Apply random rotation to the polygon region
        angle = randi([0, nSides - 1]) * 360 / nSides;
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
    end
end
