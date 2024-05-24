function imgOut = randomCircleRotation(imgIn, nCircles, radIn)
    % Get the dimensions of the input image
    [rows, cols, channels] = size(imgIn);
    
    % Initialize the output image
    imgOut = imgIn;
    
    % Create a mask for the circle
    [X, Y] = meshgrid(-radIn:radIn, -radIn:radIn);
    circleMask = (X.^2 + Y.^2) <= radIn^2;
    
    % Prepare mask for blending
    mask = double(circleMask);
    
    for i = 1:nCircles
        % Generate random center for the circle
        centerX = randi([radIn + 1, cols - radIn - 1]);
        centerY = randi([radIn + 1, rows - radIn - 1]);
        
        % Extract the circular region from the image
        circleRegion = zeros(2*radIn+1, 2*radIn+1, channels, 'like', imgIn);
        for c = 1:channels
            for x = -radIn:radIn
                for y = -radIn:radIn
                    if circleMask(y + radIn + 1, x + radIn + 1)
                        circleRegion(y + radIn + 1, x + radIn + 1, c) = ...
                            imgIn(centerY + y, centerX + x, c);
                    end
                end
            end
        end
        
        % Apply random rotation to the circular region
        angle = randi([0, 360]);
        rotatedRegion = imrotate(circleRegion, angle, 'bilinear', 'crop');
        rotatedMask = imrotate(mask, angle, 'bilinear', 'crop');
        
        % Normalize the rotated mask
        rotatedMask(rotatedMask > 0) = 1;
        
        % Place the rotated circular region back in the image
        for c = 1:channels
            for x = -radIn:radIn
                for y = -radIn:radIn
                    if rotatedMask(y + radIn + 1, x + radIn + 1)
                        imgOut(centerY + y, centerX + x, c) = ...
                            rotatedRegion(y + radIn + 1, x + radIn + 1, c);
                    end
                end
            end
        end
    end
end
