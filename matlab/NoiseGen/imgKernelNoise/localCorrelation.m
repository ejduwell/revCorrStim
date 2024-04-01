function correlation_img = localCorrelation(img1, img2, kernel_size)
    % Ensure the two images are of the same size
    assert(all(size(img1) == size(img2)), 'The two images must have the same dimensions.');
    
    % Pad the images to handle edge cases
    pad_size = floor(kernel_size/2);
    img1_padded = padarray(img1, [pad_size pad_size], 'symmetric');
    img2_padded = padarray(img2, [pad_size pad_size], 'symmetric');
    
    % Initialize the output correlation image
    [rows, cols] = size(img1);
    correlation_img = zeros(rows, cols);
    
    for i = 1:rows
        for j = 1:cols
            % Extract local patches from both images
            patch1 = double(img1_padded(i:i+2*pad_size, j:j+2*pad_size));
            patch2 = double(img2_padded(i:i+2*pad_size, j:j+2*pad_size));
            
            % Compute correlation coefficient for the patches
            c = corrcoef(patch1(:), patch2(:));
            correlation_img(i, j) = c(1, 2); % c is a 2x2 matrix; we want the off-diagonal element
        end
    end
end
