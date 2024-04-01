function mask = create_border_mask(img, distance)
    % Function to create a mask where pixels within a specified distance 
    % from the border of the image are set to 0, and the rest are set to 1.

    % Input:
    %   img: The input image.
    %   distance: The specified distance from the border.

    % Output:
    %   mask: The resultant mask with desired pixels set to 0 and 1.

    [rows, cols, ~] = size(img);  % Assuming the image might be color or grayscale
    mask = ones(rows, cols);

    % Setting the border pixels to 0
    mask(1:distance, :) = 0;
    mask(end-distance+1:end, :) = 0;
    mask(:, 1:distance) = 0;
    mask(:, end-distance+1:end) = 0;
end
