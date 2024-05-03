function adjusted_img = adjustMeanLum(img,desiredMean)

% Load your image
%impath="/Users/eduwell/matlabSandBox/texture/textureSynth/testNoise13/wtdAvgSmplz1-25_r36.png";
%img = imread(impath); % Replace with your image path

% Convert the image to grayscale if it's a color image
if size(img, 3) == 3
    img = rgb2gray(img);
end

%--------------------
%figure
%imshow(img)
%--------------------

% Convert the image to double for processing
img_double = double(img);

% Subtract the current mean from all pixels
img_zero_mean = img_double - mean(img_double(:));

% Calculate the current mean
%fprintf('The original mean pixel value is: %f\n', mean(img_double(:)));

% Scale the pixel values to have a range of -127 to 127
img_scaled = img_zero_mean * (desiredMean / max(abs(img_zero_mean(:))));

% Add 127 to all pixel values
adjusted_img_double = img_scaled + desiredMean;

% Convert the pixel values back to uint8
adjusted_img = uint8(adjusted_img_double);

% Verify new mean value
new_mean = mean(adjusted_img(:));
%fprintf('The new mean pixel value is: %f\n', new_mean);

%--------------------
%figure
%imshow(adjusted_img)
%--------------------
end