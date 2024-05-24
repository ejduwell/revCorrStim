function [nzImz] = kernelNoise_v1_4(imgsPerKrnl,desiredSize,krnlNzBIz,nReps,krnlImgWgts,smoothTiles,gSmthK1sz,smoothFinal,gSmthK2sz)

rng("shuffle"); % Make sure the random number generator is shuffled..
meanLumIn=mean(krnlNzBIz,"all"); % compute mean luminance of input image..
maxLumIn=max(max(max(krnlNzBIz)));
minLumIn=min(min(min(krnlNzBIz)));

krnlNzBI=krnlNzBIz(:,:,1); % grab the top image as instance for computing kernel sizes..
%% Subsample the image with kernelSubSampler

tic
totKimz=sum(imgsPerKrnl); % count up total number of intermediate images..
nzImz=uint8(zeros(desiredSize(1),desiredSize(2),nReps));

parfor kk = 1:nReps
%for kk = 1:nReps

[adjImg, kSampleArrayz, kernelDim_ImFrac, kernelSqrDims, kernelPosCnts] = kernelSubSampler_v3(krnlNzBI);

%% Generate intermediate "Tile Images" by tiling the image space with kernel subsamples

% get the number of kernels.. 
nKernels=length(kernelSqrDims);
kNumberz=linspace(1,nKernels,nKernels); % set kNumberz equal to list 1-->nKernels..

kernelFrames=zeros(size(adjImg,1),size(adjImg,2),totKimz); % preallocate stack for intermediate kernel frames..
countr=1; % initialize countr..

% Loop through the kernels and generate kernel subsample frames for each..
for ii = 1:nKernels


kNumber=kNumberz(ii);

for nn = 1:imgsPerKrnl(ii)
% pull info for this kernel..
kDms=[kernelSqrDims(kNumber),kernelSqrDims(kNumber)];
kFrkn=kernelDim_ImFrac(kNumber);

tileImg=uint8(zeros(size(adjImg,1),size(adjImg,2))); % preallocate..
imsz=size(adjImg,1);

for yy = 1:kFrkn
    yStrt=1+((imsz/kFrkn)*(yy-1));
    for xx = 1:kFrkn

        xStrt=1+((imsz/kFrkn)*(xx-1));

        % select a random kernel window...
        kIndx=randi(kernelPosCnts(kNumber)); % grab a random kernel position index..
        
        % get coordinates for kernel at that index
        kCoords=kSampleArrayz{1,kNumber}(kIndx,:);

        % use coordinates to grab sub-image sample from base image
        %kSubImg=uint8(adjImg(kCoords(2):kCoords(4),kCoords(1):kCoords(3)));
        %kSubImg=rotateKrnlWindow(kCoords,adjImg);
        
        krotlimt=512;
        randImIdx=randi(size(krnlNzBIz,3));
        if (max(size(tileImg)) / max(kDms)) < krotlimt
        kSubImg=rotateKrnlWindow2(kCoords,krnlNzBIz(:,:,randImIdx));
        else
         kSubImg=uint8(krnlNzBIz(kCoords(2):kCoords(4),kCoords(1):kCoords(3),randImIdx));   
        end

        %kSubImg=kSampleArrayz{1,kNumber}(:,:,kIndx); % grab rand
        
        % select a random flip condition..
        flipCon1=randi(2); % horz flip
        flipCon2=randi(2); % vert flip
        flipCon3=randi(2); % invert

        if flipCon1==2
        % Flip the image horizontally
        kSubImg = fliplr(kSubImg);
        end
        
        if flipCon2==2
        % Flip the image vertically
        kSubImg = flipud(kSubImg);
        end

        if flipCon3==2
        % Flip the image both horizontally and vertically
        kSubImg = kSubImg';
        end

        tileImg(yStrt:(yStrt+kDms(1)-1),xStrt:(xStrt+kDms(1)-1))=kSubImg;
    end

end

% Convert to double and scale from 0-1
tileImg = double(tileImg);
maxVal=max(max(tileImg));
tileImg = tileImg./maxVal; % divide all by max value to rescale range from 0-1

% Apply random shifts to jitter the position of the tile boundaries 
vShftOptions={'up','down'};
hShftOptions={'left','right'};
shftCon1=randi(2); % vert shft random idx selector
shftCon2=randi(2); % horz shft random idx selector
vShftDir=vShftOptions{shftCon1};
hShftDir=hShftOptions{shftCon2};
vShftVal=randi(imsz);
hShftVal=randi(imsz);
tileImg = shiftImage(tileImg, vShftDir, vShftVal); % apply vertical shift
tileImg = shiftImage(tileImg, hShftDir, hShftVal); % apply horizontal shift

%Apply Random Circle Rotations if requested
klimt=0;
if (max(size(tileImg)) / max(kDms)) < klimt
randCircRotn = 1;
else
randCircRotn = 0;    
end

if randCircRotn==1
% nCircles=10;
% radIn=(size(tileImg,1)/4);
% tileImg = randomCircleRotation(tileImg, nCircles, radIn);

nPolygons=((max(size(tileImg))/max(kDms))^2)/4;
nSides=8;
%sizeIn=round(size(tileImg,1)/4);
sizeIn=round((max(kDms)/2)*1.5);
nRotTimes=1;

for ff=1:nRotTimes
tileImg = randomPolygonRotation5(tileImg, nPolygons, nSides, sizeIn);
end

end

% smooth tile image if requested..
if smoothTiles(ii)==1
tileImg = imgaussfilt(tileImg,gSmthK1sz(ii));
end

% apply requested weighting
tileImg = tileImg.*krnlImgWgts(ii);

% save tile image frame into kernelFrames stack
kernelFrames(:,:,countr)=tileImg;
countr=countr+1;

end

end
%==========================================================================

%% Combine the stack of kernel subsample images linearly to form noise image

nSubImgs=size(kernelFrames,3); % get total number of intermediate kernel images..
noiseImg=sum(kernelFrames,3)./nSubImgs; % take weighted average of the stack..
% smooth final noise image if requested..
if smoothFinal==1
noiseImg = imgaussfilt(noiseImg,gSmthK2sz);
end


%Adjust mean luminance to equal input base image.. 
noiseImg = adjustMeanLum(noiseImg,meanLumIn);
% scale values of final image to range from min to max input values and convert to uint8..
noiseImg=uint8(rescale(noiseImg,minLumIn,maxLumIn)); 

nzImz(:,:,kk)=noiseImg;
end
disp(" ")
disp("Making Kernel Noise Imz Took:")
toc
end