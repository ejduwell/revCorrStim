function shiftedImage = shiftImage(image, direction, n)
    % Validate the direction input
    if ~ismember(direction, {'up', 'down', 'left', 'right'})
        error('Direction must be ''up'', ''down'', ''left'', or ''right''.');
    end

    % Validate the shift size
    if n < 0 || n ~= round(n)
        error('Shift size must be a non-negative integer.');
    end

    % Get the image dimensions
    [rows, cols, ~] = size(image);

    % Perform the shifting based on the direction
    switch direction
        case 'right'
            n = mod(n, cols); % Ensure n is within the image bounds
            shiftedImage = [image(:,end-n+1:end,:), image(:,1:end-n,:)];
        case 'left'
            n = mod(n, cols);
            shiftedImage = [image(:,n+1:end,:), image(:,1:n,:)];
        case 'down'
            n = mod(n, rows);
            shiftedImage = [image(end-n+1:end,:,:); image(1:end-n,:,:)];
        case 'up'
            n = mod(n, rows);
            shiftedImage = [image(n+1:end,:,:); image(1:n,:,:)];
        otherwise
            % Should never reach here due to initial validation
            error('Invalid direction.');
    end
end
