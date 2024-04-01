% Directory containing the face images
%imageDir = '/home/eduwell/Matlab_Sandbox/imgKernelNoise/baseImgs/londonFaces/neutral_front';
imageDir = '/home/eduwell/Matlab_Sandbox/imgKernelNoise/baseImgs/londonFaces/smiling_front';
imageFiles = dir(fullfile(imageDir, '*.jpg')); % Assuming images are in JPG format

% specify desired output image size
dimsOut=[200, 200];

% Initialize variables for averaging
sumImage1 = 0;
sumImage2 = 0;
numImages = length(imageFiles);

% preallocate array to store individual aligned faces
alndImgArray=zeros(dimsOut(1),dimsOut(2),numImages);
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
        sumImage1 = sumImage1 + faceImg;
    end
end

%% Compute the average image
avgImage1 = sumImage1 / numImages;

%% Display the average face
figure
imshow(avgImage1);
title('Average Face1');

%% Now use facial landmarks to align all to average to get better feature alignment

% Initialize an array to store reference landmarks
%referenceLandmarks = []; % This will be filled with the landmarks of the first image
% Detect facial landmarks
referenceLandmarks = detectFacialLandmarks(avgImage1); % Implement this function based on your chosen method

% Loop through each image
for i = 1:numImages

    % Read image
    %img = imread(fullfile(imageDir, imageFiles(i).name));
    img = alndFaceArray(:,:,i);

    % Detect facial landmarks
    landmarks = detectFacialLandmarks(img); % Implement this function based on your chosen method
    
    % Align image based on landmarks
    alignedImg = alignImagesUsingLandmarks(img, referenceLandmarks, landmarks);

    % Convert to double for averaging
    alignedImg = im2double(alignedImg);
    

    % store aligned image..
    alndImgArray(:,:,i)=alignedImg; 

    % Sum the images
    sumImage2 = sumImage2 + alignedImg;
    
end

%% Compute the average image
avgImage2 = sumImage / numImages;

%% Display the average face
figure
imshow(avgImage2);
title('Average Face2');

%% Show a random individual face individual face
rng("shuffle");
randIdx=randi(numImages);

figure
imshow(alndImgArray(:,:,randIdx));
title('Individual Aligned Face Sample');
