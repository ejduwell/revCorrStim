function imgOut = localCorrelation_v3(img1, img2, kernel_sizes)

nKernels=length(kernel_sizes);
imMat=zeros(size(img1,1),size(img1,2),nKernels);

for qq=1:nKernels

    kernel_size=kernel_sizes(qq);% assign kernel size for this pass.. 

    % Ensure the two images are of the same size
    assert(all(size(img1) == size(img2)), 'The two images must have the same dimensions.');
    
    % Convert images to double
    img1 = double(img1);
    img2 = double(img2);
    
    % Pad the images to handle edge cases
    pad_size = floor(kernel_size/2);
    img1_padded = padarray(img1, [pad_size pad_size], 'symmetric');
    img2_padded = padarray(img2, [pad_size pad_size], 'symmetric');
%     img1_padded = padarray(img1, [pad_size pad_size], 'replicate');
%     img2_padded = padarray(img2, [pad_size pad_size], 'replicate');
    
    % Compute mean images
    mean_filter = ones(kernel_size) / (kernel_size^2);
    mean1 = conv2(img1_padded, mean_filter, 'same');
    mean2 = conv2(img2_padded, mean_filter, 'same');
    
    % Compute local products
    product_img = conv2(img1_padded .* img2_padded, mean_filter, 'same');
    
    % Compute local variances
    var1 = conv2(img1_padded.^2, mean_filter, 'same') - mean1.^2;
    var2 = conv2(img2_padded.^2, mean_filter, 'same') - mean2.^2;
    
    % Compute the local correlation
    numerator = product_img - mean1 .* mean2;
    denominator = sqrt(var1 .* var2);
    
    % Avoid division by zero
    denominator(denominator == 0) = 1;
    
    correlation_img = numerator ./ denominator;

    % Crop the output image to retain only the valid region of the same size as input
    % also: run "real" to only output the real values of the complex image
    correlation_img = real(correlation_img(pad_size+1:end-pad_size, pad_size+1:end-pad_size));
    
    % mask out pixels with the kernel's width distance from the border..
    % doing this to eliminate potential artifacts from kernel near edge of
    % image..
%     mask = create_borkernel_sizesder_mask(correlation_img, kernel_size);
%     correlation_img=correlation_img.*mask;

    % save image in imMat
    imMat(:,:,qq)=correlation_img;

end

imgOut=mean(imMat,3);

end


