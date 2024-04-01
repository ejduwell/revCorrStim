function [kernelFrames, betterSSIMz, bestSSIMarray]= ssimKrnlNzReArrange(baseIm, noiseIm,kSampleArrayz, kernelDim_ImFrac, kernelSqrDims, kernelPosCnts,initSSIM)
% put the features in the noise image in the location with the highest ssim
% value relative to the base image.. ie put it them the most structurally
% similar location in the base image..

nScramblz=1000; % number of times each kernel is used to make a scrambled image to compare against BI via ssim

baseIm=double(baseIm); 
noiseIm=double(noiseIm);

rng("shuffle"); % Make sure the random number generator is shuffled..

% get the number of kernels.. 
nKernels=length(kernelSqrDims);
kNumberz=linspace(1,nKernels,nKernels); % set kNumberz equal to list 1-->nKernels..

kernelFrames=zeros(size(baseIm,1),size(baseIm,2),nKernels); % preallocate stack for intermediate kernel frames..
countr=1; % initialize countr..

bestSSIMarray=zeros(1,nKernels);

% Loop through the kernels and generate kernel subsample frames for each..
for ii = 1:nKernels

kNumber=kNumberz(ii);

% pull info for this kernel..
kDms=[kernelSqrDims(kNumber),kernelSqrDims(kNumber)];
kFrkn=kernelDim_ImFrac(kNumber);

imsz=size(baseIm,1);
nKwindowz=kFrkn*kFrkn;

kStack=zeros(kDms(1),kDms(2),nKwindowz); % preallocate for kernel window stack
countr2=1; % initialize counter number 2..
% build stack of kernel windows for this kernel size..
for yy = 1:kFrkn
    yStrt=1+((imsz/kFrkn)*(yy-1));
    for xx = 1:kFrkn

        xStrt=1+((imsz/kFrkn)*(xx-1));

        % select the next non-overlapping kSubImg tile for this kernel size...
        kSubImg=noiseIm(yStrt:(yStrt+kDms(1)-1),xStrt:(xStrt+kDms(1)-1));

        % add it to the stack..
        kStack(:,:,countr2)=kSubImg;
        countr2=countr2+1; % update countr2 for next pass..
    end
end

% Now generate nScramblz number of scrambled noise images using this kernel
% size as the "scrambled chunk size" (ie the size chunk of the image that 
% is getting rearranged/re-tiled).
%
% On each pass, compute the ssim of the scrambled version with the base
% image.. keep only the one with the highest ssim value..
%||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
bestSSIM=0; % initialize variable for tracking the best ssim value
stackInx=linspace(1,nKwindowz,nKwindowz); % create vector of index values from 1->the number of imgs in kStack

for ee = 1:nScramblz
    % preallocate array..
    nzScrmld=zeros(size(baseIm,1),size(baseIm,2)); 
    % generate randomized index array for scrambling the order in which the
    % windows in kStack are re-assembled into the noise frame (nzScrmld)
    perm = randperm(nKwindowz); % Generate a random permutation of indices
    ranStackInx = stackInx(perm); % Create the new vector with scrambled positions

    % Now loop through and assign the kernel-windowed images in kStack to
    % their randomized positions..
    countr2=1; % re-initialize counter number 2..
    % build stack of kernel windows for this kernel size..
    for yy = 1:kFrkn
        yStrt=1+((imsz/kFrkn)*(yy-1));
        for xx = 1:kFrkn
            xStrt=1+((imsz/kFrkn)*(xx-1));
            % get the random index for this pass
            idxTmp=ranStackInx(countr2);
            % assign the kernel windowed image at this index in kStack to
            % the next position in nzScrmld
            nzScrmld(yStrt:(yStrt+kDms(1)-1),xStrt:(xStrt+kDms(1)-1))=kStack(:,:,idxTmp);
            countr2=countr2+1; % update countr2 for next pass..
        end
    end

    % Apply random shifts to jitter the position of the tiles 
    vShftOptions={'up','down'};
    hShftOptions={'left','right'};
    shftCon1=randi(2); % vert shft random idx selector
    shftCon2=randi(2); % horz shft random idx selector
    vShftDir=vShftOptions{shftCon1};
    hShftDir=hShftOptions{shftCon2};
    vShftVal=randi(imsz);
    hShftVal=randi(imsz);
    nzScrmld = shiftImage(nzScrmld, vShftDir, vShftVal); % apply vertical shift
    nzScrmld = shiftImage(nzScrmld, hShftDir, hShftVal); % apply horizontal shift


    % compute ssim between nzScrmld and the base image for this pass
    ssimVal = ssim(nzScrmld,baseIm);

    % If this is the first pass: just assign ssimVal to bestSSIM
    % and nzScrmld to bestSSIM_img..
    if ee == 1
        bestSSIM=ssimVal;
        bestSSIM_img=nzScrmld;
    % If not: check whether this ssimVal beats the current best..
    % if so, update bestSSIM and save this nzScrmld as 
    % bestSSIM_img
    else
        if ssimVal>bestSSIM
        bestSSIM=ssimVal;
        bestSSIM_img=nzScrmld;  
        end
    end

    
    

end
%||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% save image with highest ssim relative to the base img (bestSSIM_img)
% into kernelFrames stack for this pass/kernel along with ssim value..
kernelFrames(:,:,countr)=bestSSIM_img;
bestSSIMarray(1,countr)=bestSSIM;

countr=countr+1; % update countr..

end
%==========================================================================

% label which kernel's best scrambled image resulted in ssim values better
% than the initial ssim between the noise and base image..
betterSSIMz=zeros(1,size(bestSSIMarray,2));
for xz=1:size(betterSSIMz,2)
if bestSSIMarray(xz)>initSSIM
    betterSSIMz(xz)=1;
else
    betterSSIMz(xz)=0;
end

end