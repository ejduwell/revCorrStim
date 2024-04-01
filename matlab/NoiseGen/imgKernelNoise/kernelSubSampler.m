function [adjImg, kSampleArrayz, kernelDim_ImFrac, kernelSqrDims, kernelPosCnts] = kernelSubSampler(imgInPath,selectSize,desiredSize)


%% Parameters
%imgInPath="/home/eduwell/Matlab_Sandbox/imgKernelNoise/face/031_03.jpg";

% Specify desired size?
%selectSize=1; % if 1, this means that we will use/resize the input image the selected size below
%desiredSize=[256,256]; % (MUST BE SQUARE IMAGE AND DIMZ MUST BE POWER OF 2)

%% Read in input image
imgIn=imread(imgInPath);

%% Get input image dimensions
imgInDimz=size(imgIn); % input get image dimension vector
imgInNdims = ndims(imgIn); % get number of dimensions

%% Ensure image is grayscale, square, and that dimensions are a power of two

disp("Checking to make sure that the input image is grayscale, square, and that dimensions are a power of two...")

% If image is not grayscale (imgInNdims>2), convert
if imgInNdims > 2
    disp("Input image is not grayscale.. converting to grayscale")
    imgIn=rgb2gray(imgIn);

    imgInNdims = ndims(imgIn); % re-get number of dimensions
    imgInDimz=size(imgIn); % re-get image dimension vector
end

% check if image is square
sqrTest=imgInDimz(1)==imgInDimz(2);

% check if minimum dimension is an exponent of 2
isPow2 = isPowerOfTwo(min(imgInDimz));

if selectSize==1
    isDesiredSize=min(imgInDimz)==desiredSize(1);
else
    isDesiredSize=1;
end

% If image is not square and/or minimum dimension is not a power of 2,
% adjust the image to make it the nearest square with dimension being
% nearest power of two to minumum dimension..
if sqrTest==0 || isPow2==0 || isDesiredSize==0
    if sqrTest==0
        disp("Input image was not square..")
    end

    if isPow2==0
        disp("Minimum dimension of input image was not a power of two..")    
    end
    
    if isDesiredSize==0
        disp("User specified a desired size and image dimensions don't match..")
    end
    disp("Adjusting image to ensure it is square with x and dimensions equal to either the desired size (if specified) or the nearest power of two")
    nearestPow2 = findNearestPowerOfTwo(min(imgInDimz));
    disp(strcat("Nearest power of two is: ", num2str(nearestPow2)))
    if isDesiredSize == 0
        imgIn = imgReformater(imgIn,desiredSize);
    else
        imgIn = imgReformater(imgIn,[nearestPow2,nearestPow2]);
    end
    imgInNdims = ndims(imgIn); % re-get number of dimensions
    imgInDimz=size(imgIn); % re-get image dimension vector
else
    disp("Image matched required dimensions... continuing")
end



%% Compute series of kernel sizes and position/subsample counts

squareDim=imgInDimz(1); % square dimension of image
nKernels=log2(squareDim); % compute number of kernels..
nKernels=nKernels-1; %subtract 1 as we don't want to use the last kernel (this will just be a pixel...)

% compute kernel img fraction vector
kernelDim_ImFrac=ones(1,nKernels);
kernelDim_ImFrac=kernelDim_ImFrac.*2;
expnts=linspace(1,nKernels,nKernels);
kernelDim_ImFrac=kernelDim_ImFrac.^expnts;

% compute kernel square dimension vector
tmpVec=ones(1,nKernels).*squareDim;
kernelSqrDims=tmpVec./kernelDim_ImFrac;

% compute kernel position/subsample count vector
kernelPosCnts = (kernelSqrDims-tmpVec).^2;


%% Preallocate matrices to store image subsamples for each kernel

kSampleArrayz=cell(1,nKernels);
for ii=1:nKernels
    kSamplesMat=uint8(zeros(kernelSqrDims(ii),kernelSqrDims(ii),kernelPosCnts(ii)));
    kSampleArrayz{1,ii}=kSamplesMat;
end
clear kSamplesMat; %clear variable to minimize RAM overhead..

%% Populate kernel subsample matrices for each kernel with their respective image subsamples..

% Note: kernel positions/image samples are stored as a function of the
% upper left kernel corner vertex of the square kernel w/in the image
for vv = 1:nKernels

    % Build image coordinate mat for this pass..
    %----------------------------------------------------------------------
    %
    % First construct full coordinate set
    % #####################################################################
    n = squareDim;
    yMat = repmat((1:n)', 1, n);
    xMat = yMat';
    xVect=xMat(:);
    yVect=yMat(:);
    imCoords=horzcat(yVect,xVect);
    % #####################################################################
    %
    % Then reduce to only acceptable coordinates for this kernel (ones in
    % which the whole kernel is contained within the image..)
    % #####################################################################
    % 1. Define your limits L1 and L2
    xLim=(squareDim-kernelSqrDims(vv));
    yLim=(squareDim-kernelSqrDims(vv));
    L1 = yLim; % set your value for L1
    L2 = xLim; % set your value for L2
    %
    % 2. Create logical arrays for each condition
    condition1 = imCoords(:, 1) <= L1;
    condition2 = imCoords(:, 2) <= L2;
    %
    % 3. Combine conditions
    combinedCondition = condition1 & condition2;
    %
    % 4. Extract rows from M that satisfy both conditions
    imCoords = imCoords(combinedCondition, :);
    % #####################################################################
    %----------------------------------------------------------------------
    
    % Now Add two more columns for the y2 and x2 coordinates..
    y2Col=(imCoords(:,1)+(kernelSqrDims(vv))-1);
    x2Col=(imCoords(:,2)+(kernelSqrDims(vv))-1);
    imCoords=horzcat(imCoords,y2Col,x2Col);
    

    % Finally, use coordinates in imCoords to grab image kernel subset for
    % each position (row in imCoords)..
    for qq = 1:size(imCoords,1)
        kSampleArrayz{1,vv}(:,:,qq)=uint8(imgIn(imCoords(qq,2):imCoords(qq,4),imCoords(qq,1):imCoords(qq,3)));
    end
    
    % (FOR DEBUGGING)
    % Loop through set of kernel images and display each in a figure
%     for ww = 1:size(imCoords,1)
%         cleae
%         pause(0.1);
%     end

end

adjImg=imgIn;

end