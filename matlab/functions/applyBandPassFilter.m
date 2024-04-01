function filteredImage = applyBandPassFilter(inputImage, filterFunction,lowFreq,highFreq)
    % Check if the input image is RGB and convert to grayscale
    if size(inputImage, 3) == 3
        inputImage = rgb2gray(inputImage);
    end

    % Convert the image to double precision for FFT
    inputImage = double(inputImage);

    % Get the dimensions of the input image
    [rows, cols] = size(inputImage);

    % Apply the Fourier Transform to the image
    fftImage = fft2(inputImage);
    fftShifted = fftshift(fftImage);

    % Create a matrix to hold the filter
    filter = zeros(rows, cols);

    % Determine the center of the image
    centerX = cols / 2;
    centerY = rows / 2;

    % Generate the band-pass filter
    for i = 1:rows
        for j = 1:cols
            % Calculate the distance from the center
            distance = sqrt((i - centerY)^2 + (j - centerX)^2);
            % Apply the input function to determine the filter value
            filter(i, j) = filterFunction(distance,lowFreq,highFreq);
        end
    end

    % Apply the filter to the frequency domain image
    filteredFFTShifted = fftShifted .* filter;

    % Shift back and apply the inverse Fourier Transform
    filteredFFT = ifftshift(filteredFFTShifted);
    filteredImageRaw = real(ifft2(filteredFFT));

    % Normalize the filtered image to have values between 0 and 255
    %filteredImage = uint8(255 * mat2gray(filteredImageRaw));
    filteredImage=filteredImageRaw; % dont normalized..
end
