function [kernelFrames, krnlWindwStax, combTallyImg, indTallyImgz]= xCorKrnlNzReArrange_v2(baseIm, noiseIm, kernelDim_ImFrac, kernelSqrDims,MinMaxCor,nKernReArr)
% put the features in the noise image in the location with the highest
% local correlation value relative to the base image.. ie put it them the most
% similar location in the base image..

% convert to double and make range from 0-1..
baseIm=rescale(double(baseIm),0,1); 
noiseIm=rescale(double(noiseIm),0,1);

rng("shuffle"); % Make sure the random number generator is shuffled..

% get the number of kernels.. 
%nKernels=length(kernelSqrDims);
nKernels=nKernReArr; % trying a pre-specified subset..

kNumberz=linspace(1,nKernels,nKernels); % set kNumberz equal to list 1-->nKernels..

kernelFrames=zeros(size(baseIm,1),size(baseIm,2),nKernels); % preallocate stack for intermediate kernel frames..
krnlWindwStax=cell(2,nKernels);

combTallyImg=zeros(size(baseIm,1),size(baseIm,2));
indTallyImgz=zeros(size(baseIm,1),size(baseIm,2),nKernels);
countr=1; % initialize countr..
% Loop through the kernels and generate kernel subsample frames for each..
for ii = 1:nKernels

kNumber=kNumberz(ii);

% pull info for this kernel..
kDms=[kernelSqrDims(kNumber),kernelSqrDims(kNumber)];
kFrkn=kernelDim_ImFrac(kNumber);

imsz=size(baseIm,1);
nKwindowz=kFrkn*kFrkn;

kStack=zeros(kDms(1),kDms(2),nKwindowz); % preallocate for kernel window stack
kStackCorData=cell(nKwindowz,4); % preallocate for associated correlation data..

countr2=1; % initialize counter number 2..
% build stack of kernel windows for this kernel size..
for yy = 1:kFrkn
    yStrt=1+((imsz/kFrkn)*(yy-1));
    for xx = 1:kFrkn

        xStrt=1+((imsz/kFrkn)*(xx-1));

        % select the next non-overlapping kSubImg tile for this kernel size...
        kSubImg=noiseIm(yStrt:(yStrt+kDms(1)-1),xStrt:(xStrt+kDms(1)-1));

        % Find the best matching position in the base image for this
        % sample..
        [MaxCorCoords, MinCorCoords, maxCor, minCor] = getBestLocMatch_v2(baseIm,kSubImg);

        % add kernel-windowed img to the stack..
        kStack(:,:,countr2)=kSubImg;

        % save corresponding correlation data too
        kStackCorData(countr2,1)={MaxCorCoords};
        kStackCorData(countr2,2)={MinCorCoords}; 
        kStackCorData(countr2,3)={maxCor};
        kStackCorData(countr2,4)={minCor};

        countr2=countr2+1; % update countr2 for next pass..

    end
end

% Save kernel windows and correlation data for this kernel pass..
krnlWindwStax(1,countr)={kStack};
krnlWindwStax(2,countr)={kStackCorData};


% Now, loop through the kernel sub windows in kStack and rearrange to their 
% respective positions of max correlation and with their
% contributions to the final values weighted by correlation value.
kernelFrame=zeros(size(baseIm,1),size(baseIm,2)); % preallocate
for aa=1:size(kStack,3)

    % pull min/max cor value and position data
    %MaxCoordzTmp=kStackCorData{aa,1};
    %MaxCorValTmp=kStackCorData{aa,3};
    %MinCoordzTmp=kStackCorData{aa,2};
    %MinCorValTmp=kStackCorData{aa,4};
    
    if MinMaxCor == 0 % codes for "use min"
        CoordsTmp=kStackCorData{aa,2};
        CorValTmp=kStackCorData{aa,3};
    elseif MinMaxCor == 1 % codes for "use max"
        CoordsTmp=kStackCorData{aa,1};
        CorValTmp=kStackCorData{aa,3};
    end

    % Check if there is more then one set of coordinates in CoordsTmp
    % (this occurs if there were multiple locations that tied for
    % "higheset/lowest" correlation value..)
    % If this is the case, grab one at random
    if size(CoordsTmp,1)>1
        ranIdxTmp=randi(size(CoordsTmp,1)); %pick random row index
        CoordsTmp=CoordsTmp(ranIdxTmp,:);%reassign CoordsTmp to only include that row..
    end

    % pull the corresonding region from the kernelFrame
    matchedRegion = kernelFrame(CoordsTmp(1):CoordsTmp(2), CoordsTmp(3):CoordsTmp(4));
    
    % weight the kernel window values by the max correlation value achieved
    % at this position
    kWinWeighted=((kStack(:,:,aa)).*abs(CorValTmp));

    % set pixel values equal to sum of kWinWeighted and what ever is
    % already there..
    kernelFrame(CoordsTmp(1):CoordsTmp(2), CoordsTmp(3):CoordsTmp(4)) = (kWinWeighted+matchedRegion);
    
    % Finally, record the fact this region matched in the individual kernel
    % and combined tally images
    % (ie add one to the pixels in this region in the tally img)
    tlyTmp=ones(size(matchedRegion,1),size(matchedRegion,2));
    combTallyImg(CoordsTmp(1):CoordsTmp(2), CoordsTmp(3):CoordsTmp(4))=((combTallyImg(CoordsTmp(1):CoordsTmp(2), CoordsTmp(3):CoordsTmp(4)))+tlyTmp);
    indTallyImgz(CoordsTmp(1):CoordsTmp(2), CoordsTmp(3):CoordsTmp(4),ii)=((indTallyImgz(CoordsTmp(1):CoordsTmp(2), CoordsTmp(3):CoordsTmp(4),ii))+tlyTmp);

end

kernelFrames(:,:,ii)=kernelFrame; % save kernel frame..

countr=countr+1; % update countr..

end
% clear arbitrary counter variables..
clear countr
clear countr2
%==========================================================================

end
