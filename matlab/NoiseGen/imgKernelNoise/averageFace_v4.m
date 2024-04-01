function [avgImage,alndFaceArray] = averageFace_v4(imageDir,dimsOut,imTag)

% Directory containing the face images
%imageDir = '/home/eduwell/Matlab_Sandbox/imgKernelNoise/baseImgs/londonFaces/neutral_front';
%imageDir = '/home/eduwell/Matlab_Sandbox/imgKernelNoise/baseImgs/londonFaces/smiling_front';
imageFiles = dir(fullfile(imageDir,imTag)); % Assuming images are in JPG format
%imageFiles = dir(fullfile(imageDir, '*.jpg')); % Assuming images are in JPG format

% specify desired output image size
%dimsOut=[200, 200];

% Initialize variables for averaging
sumImage = 0;
numImages = length(imageFiles);

% preallocate array to store individual aligned faces
alndFaceArray=zeros(dimsOut(1),dimsOut(2),numImages);

% Load a pre-trained face detector
%faceDetector = vision.CascadeObjectDetector;
detector = buildDetector();


%specify pretrain net for classifier that validates face cropping..
%netStr="alexnet";

% Loop through each image
for i = 1:numImages

    % Read image
    img = imread(fullfile(imageDir, imageFiles(i).name));
    
    % Detect facial features
    %bbox = step(faceDetector, img);
    [bbox, bbimg, faces, bbfaces] = detectFaceParts(detector,img,2);
    

    % Select the first detected face (if multiple faces are in the image)
    if ~isempty(faces)

        %bbox = bbox(1, :); 

        % If there is more than 1 box, cycle through and pick the one which
        % most closely correlates with
                      
        % Crop and resize image to a standard size
        %faceImg = imcrop(img, bbox);
        
        % check number of face images..
        nfaces=max(size(faces));
        
        if nfaces > 1
        pauz="";
        else
        faceImg = imresize(faces{1,1}, dimsOut); % Example size, adjust as needed
        end

        % googlenet classifier..
        %[labelOut,labelConf,scores] = NetClassifier(faceImg,netStr);

        % Convert to double for averaging
        faceImg = im2double(faceImg);

        % Convert image to gray-scale if it has color channels..
        if size(faceImg,3)>1
        faceImg=rgb2gray(faceImg);
        end

        % store aligned face..
        alndFaceArray(:,:,i)=faceImg; 

        % Sum the images
        sumImage = sumImage + faceImg;
    end
end

%% Compute the average image
avgImage = sumImage / numImages;

%% Display the average face
% figure
% imshow(avgImage);
% title('Average Face');

%% Show a random individual face individual face
% rng("shuffle");
% randIdx=randi(numImages);
% 
% figure
% imshow(alndFaceArray(:,:,randIdx));
% title('Individual Aligned Face Sample');

end
