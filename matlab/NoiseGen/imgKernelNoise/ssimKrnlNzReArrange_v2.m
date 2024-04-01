function [kernelFrames, krnlWindwStax]= ssimKrnlNzReArrange_v2(baseIm, noiseIm,kSampleArrayz, kernelDim_ImFrac, kernelSqrDims, kernelPosCnts)
% put the features in the noise image in the location with the highest
% local correlation value relative to the base image.. ie put it them the most
% similar location in the base image..

% convert to double..
baseIm=double(baseIm); 
noiseIm=double(noiseIm);

rng("shuffle"); % Make sure the random number generator is shuffled..

% get the number of kernels.. 
nKernels=length(kernelSqrDims);
kNumberz=linspace(1,nKernels,nKernels); % set kNumberz equal to list 1-->nKernels..

kernelFrames=zeros(size(baseIm,1),size(baseIm,2),nKernels); % preallocate stack for intermediate kernel frames..
krnlWindwStax=cell(2,nKernels);

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

        % rescale sample to range 0-1 and add img it to the stack..
        kStack(:,:,countr2)=rescale(kSubImg,0,1);

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

    % pull max cor value and position data
    MaxCoordzTmp=kStackCorData{aa,1};
    MaxCorValTmp={aa,3};
    
    % pull the corresonding region from the kernelFrame
    matchedRegionMax = kernelFrame(MaxCoordzTmp(1):MaxCoordzTmp(2), MaxCoordzTmp(3):MaxCoordzTmp(4));
    
    % weight the kernel window values by the max correlation value achieved
    % at this position
    kWinWeighted=((kStack(:,:,aa)).*MaxCorValTmp);

    % set pixel values equal to sum of kWinWeighted and what ever is
    % already there..
    kernelFrame(MaxCoordzTmp(1):MaxCoordzTmp(2), MaxCoordzTmp(3):MaxCoordzTmp(4)) = kWinWeighted+matchedRegionMax;
end

kernelFrames(:,:,ii)=kernelFrame; % save kernel frame..

countr=countr+1; % update countr..

end
% clear arbitrary counter variables..
clear countr
clear countr2
%==========================================================================


end