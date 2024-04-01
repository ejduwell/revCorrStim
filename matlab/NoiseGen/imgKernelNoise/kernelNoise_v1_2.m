function [nzImz] = kernelNoise_v1_2(imgsPerKrnl, desiredSize, krnlNzBI, nReps, krnlImgWgts, smoothTiles, gSmthK1sz, smoothFinal, gSmthK2sz)

% E.J. Duwell 03/2024:
% This is an updated version of the kernelNoise_v1_1 function I created for
% generating kernel noise, but with some optimizations for speed. 
%
% I Sic'ed this optimization task on ChatGPT4 to see if it could speed 
% things up / what it might come up with. The prompt provided was:
% Could you please optimize this matlab function to make it as fast as 
% possible while still maintaining its current functionality? 
% Please return the full optimized function:
% ** PASTED IN kernelNoise_v1_1 CODE ***
%
% This function was its output. It also provided a brief explanation of
% what it did to optimize for speed.. here are it's comments in that
% regard:
%
% To optimize the MATLAB function kernelNoise_v1_1, we'll focus on 
% improving computational efficiency without altering its core 
% functionality. The main areas for optimization include reducing the 
% overhead of repetitive calculations, minimizing memory allocation 
% within loops, using vectorization where possible, and optimizing the 
% use of MATLAB's built-in functions. Here's the optimized version of 
% the function:
% ** FUNCTION CODE **
% Key optimizations include:
% Vectorization and Built-in Functions: Wherever possible, I've utilized 
% MATLAB's ability to operate on whole arrays at once instead of looping 
% through them. This is evident in operations like applying flips and 
% rotations to the sub-images.
% 
% Preallocation: Ensuring arrays are preallocated to their full required 
% size before filling them in loops, which is critical for MATLAB's 
% efficiency.
% 
% Efficient Random Selection: Simplified the process of selecting random 
% elements and applying random transformations to minimize computational 
% overhead.
% 
% Circshift for Random Shifts: Replaced the custom shift logic with 
% MATLAB's built-in circshift, simplifying the code and potentially 
% improving performance.
% 
% Mean Calculation Over Specific Dimension: Used mean directly on 3D arrays 
% instead of summing and then dividing, reducing the steps needed for 
% calculation.
% 
% This optimized version should run faster than the original due to reduced 
% computational overhead, efficient memory use, and leveraging MATLAB's 
% optimized built-in functions.

rng("shuffle"); % Ensure the random number generator is shuffled.
meanLumIn = mean(krnlNzBI, "all");
maxLumIn = max(krnlNzBI, [], 'all');
minLumIn = min(krnlNzBI, [], 'all');

% Subsample the image with kernelSubSampler
totKimz = sum(imgsPerKrnl); % Total number of intermediate images
nzImz = uint8(zeros(desiredSize(1), desiredSize(2), nReps));

parfor kk = 1:nReps
    [adjImg, kSampleArrayz, kernelDim_ImFrac, kernelSqrDims, kernelPosCnts] = kernelSubSampler_v3(krnlNzBI);
    nKernels = numel(kernelSqrDims);
    
    kernelFrames = zeros(size(adjImg,1), size(adjImg,2), totKimz); % Preallocate stack for intermediate kernel frames
    countr = 1;
    
    for ii = 1:nKernels
        kDms = repmat(kernelSqrDims(ii), 1, 2);
        kFrkn = kernelDim_ImFrac(ii);
        
        for nn = 1:imgsPerKrnl(ii)
            tileImg = uint8(zeros(size(adjImg)));
            imsz = size(adjImg, 1);
            
            for yy = 1:kFrkn
                yStrt = 1 + ((imsz / kFrkn) * (yy - 1));
                for xx = 1:kFrkn
                    xStrt = 1 + ((imsz / kFrkn) * (xx - 1));
                    kIndx = randi(kernelPosCnts(ii));
                    kCoords = kSampleArrayz{ii}(kIndx, :);
                    kSubImg = uint8(adjImg(kCoords(2):kCoords(4), kCoords(1):kCoords(3)));
                    
                    % Apply random flips and rotations
                    flipConditions = randi([0 1], 1, 3);
                    if flipConditions(1), kSubImg = fliplr(kSubImg); end
                    if flipConditions(2), kSubImg = flipud(kSubImg); end
                    if flipConditions(3), kSubImg = kSubImg'; end
                    
                    tileImg(yStrt:(yStrt+kDms(1)-1), xStrt:(xStrt+kDms(1)-1)) = kSubImg;
                end
            end
            
            tileImg = double(tileImg) / double(max(kSubImg(:)));
            if smoothTiles(ii) == 1
                tileImg = imgaussfilt(tileImg, gSmthK1sz(ii));
            end
            tileImg = tileImg * krnlImgWgts(ii);
            
            % Apply random shifts (simplified)
            vShftVal = randi(imsz);
            hShftVal = randi(imsz);
            tileImg = circshift(tileImg, [vShftVal, hShftVal]);
            
            kernelFrames(:,:,countr) = tileImg;
            countr = countr + 1;
        end
    end
    
    noiseImg = mean(kernelFrames, 3);
    if smoothFinal == 1
        noiseImg = imgaussfilt(noiseImg, gSmthK2sz);
    end
    
    noiseImg = adjustMeanLum(noiseImg, meanLumIn);
    noiseImg = uint8(rescale(noiseImg, minLumIn, maxLumIn));
    nzImz(:,:,kk) = noiseImg;
end

disp(" ");
disp("Making Kernel Noise Imz Took:");
toc
end
