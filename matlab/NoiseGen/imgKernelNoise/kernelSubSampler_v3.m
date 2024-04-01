function [adjImg, kSampleArrayz, kernelDim_ImFrac, kernelSqrDims, kernelPosCnts] = kernelSubSampler_v3(imgIn)

%% Get input image dimensions
imgInDimz=size(imgIn); % input get image dimension vector
%imgInNdims = ndims(imgIn); % get number of dimensions

%% Ensure image is grayscale, square, and that dimensions are a power of two
% THIS VERSION NOW ASSUMES THE DIMENSIONS ARE VIABLE COMING IN..
% must run "imgReformater.m" prior to input..

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
%for ii=1:nKernels
%    kSamplesMat=uint8(zeros(kernelSqrDims(ii),kernelSqrDims(ii),kernelPosCnts(ii)));
%    kSampleArrayz{1,ii}=kSamplesMat;
%end
%clear kSamplesMat; %clear variable to minimize RAM overhead..

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
    

    % Finally, save coordinates in imCoords from this image kernel into
    % kSampleArrayz output variable..

    % OLD METHOD OF SAVING ACTUAL IMAGE SUBSAMPLES (commented):
    %for qq = 1:size(imCoords,1)
    %    kSampleArrayz{1,vv}(:,:,qq)=uint8(imgIn(imCoords(qq,2):imCoords(qq,4),imCoords(qq,1):imCoords(qq,3)));
    %end
    
    kSampleArrayz{1,vv}=imCoords; % save coordinates from this kernel/pass

    % (FOR DEBUGGING)
    % Loop through set of kernel images and display each in a figure
%     for ww = 1:size(imCoords,1)
%         clear
%         pause(0.1);
%     end

end

adjImg=imgIn;

end