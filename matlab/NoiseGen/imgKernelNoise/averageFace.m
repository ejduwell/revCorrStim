% Directory containing the face images
%imageDir = '/home/eduwell/Matlab_Sandbox/imgKernelNoise/baseImgs/londonFaces/neutral_front';
imageDir = '/home/eduwell/Matlab_Sandbox/imgKernelNoise/baseImgs/londonFaces/smiling_front';
imageFiles = dir(fullfile(imageDir, '*.jpg')); % Assuming images are in JPG format

% specify desired output image size
dimsOut=[200, 200];

% Initialize variables for averaging
sumImage = 0;
numImages = length(imageFiles);

% preallocate array to store individual aligned faces
alndFaceArray=zeros(dimsOut(1),dimsOut(2),numImages);

% Load a pre-trained face detector
faceDetector = vision.CascadeObjectDetector;

% Loop through each image
for i = 1:numImages
    % Read image
    img = imread(fullfile(imageDir, imageFiles(i).name));
    
    % Convert image to gray-scale if it has color channels..
    if size(img,3)>1
    img=rgb2gray(img);
    end

    % Detect facial features
    bbox = step(faceDetector, img);
    
    % Select the first detected face (if multiple faces are in the image)
    if ~isempty(bbox)
        bbox = bbox(1, :);
        
        % Crop and resize image to a standard size
        faceImg = imcrop(img, bbox);
        faceImg = imresize(faceImg, dimsOut); % Example size, adjust as needed

        % Convert to double for averaging
        faceImg = im2double(faceImg);

        % store aligned face..
        alndFaceArray(:,:,i)=faceImg; 

        % Sum the images
        sumImage = sumImage + faceImg;
    end
end

%% Compute the average image
avgImage = sumImage / numImages;

%% Display the average face
figure
imshow(avgImage);
title('Average Face');

%% Show a random individual face individual face
rng("shuffle");
randIdx=randi(numImages);

figure
imshow(alndFaceArray(:,:,randIdx));
title('Individual Aligned Face Sample');
