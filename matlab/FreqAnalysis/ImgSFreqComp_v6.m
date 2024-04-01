function meanError = ImgSFreqComp_v6(varargin)


%noise = ["noiseData.NoiseType = 'white';","noiseData.mu = 0;","noiseData.std = 1;","noiseData.imgSize   = 512;"];   % Size of noise image (enter 2 values if x & y dims are different, ex: [256,512] )
nReps = getGlobalNumSamplz;

err_mat = zeros(1,nReps);
for jj = 1:nReps
% Get base and noise images
[im1, noiseDataLocal] = noiseBaseImDescFile(0); % call with option 0 to get the global variables..
% What Type Noise?
noizType = noiseDataLocal.NoiseType; % check the type of noise and save locally as "noizType"..
inputs = varargin{1:end};


prz2opt = getGlobalOptPars;

% --- Required Input1 (base image) --- %
% im1_in = "/Users/eduwell/OneDrive - mcw.edu/Documents/SNAP/projects/Rev_Corr/code/RevCorr/working_dir/out_test/revcorrBI_occ_0_ori_z_CR1CRpw1CRpl75CRal1CRob2CRdm3214_T1spheresTr11_5Tr20_5T_al1_L1l130l230LWT0.7L_al0_7.png";
%im1 = imread(baseIm);
if length(size(im1))>2 % check if it has a color channel..
im1 = rgb2gray(im1); % convert to greyscale..
end
im1_size = size(im1);

% --- Required Input2 (noise image) --- %
if noiseDataLocal.imgSize == "match"
%Make the size match the base image..
noiseDataLocal.imgSize   = [im1_size(1), im1_size(2)];      % Size of noise image (enter 2 values if x & y dims are different, ex: [256,512] )
end

% replace parameters in struct with the input parameters to be optimized.. 
% .. based on the type of noise we're running..
switch  noizType
    % if NoiseType is "gabor" update the gaussian blobs stuff...
    case "gabor"
       for ii = 1: length(inputs)
           par2update = inputs(ii);
           parName = prz2opt(ii);
           switch parName
               case "sig1"
                   noiseDataLocal.blobs.sig1 = par2update;
               case "sig2"
                   noiseDataLocal.blobs.sig2 = par2update;
               case "sig3"
                   noiseDataLocal.blobs.sig3 = par2update;
               case "gaussNum1"
                   noiseDataLocal.blobs.gaussNum1 = par2update;
               case "gaussNum2"
                   noiseDataLocal.blobs.gaussNum2 = par2update;
               case "gaussNum3"
                   noiseDataLocal.blobs.gaussNum3 = par2update;
               case "n_gauss1"
                   noiseDataLocal.blobs.n_gauss1 = par2update;
               case "n_gauss2"
                   noiseDataLocal.blobs.n_gauss2 = par2update;
               case "n_gauss3"
                   noiseDataLocal.blobs.n_gauss3 = par2update;
               case "grn"
                   noiseDataLocal.blobs.grn = par2update;
           end
       end 

    case "white"
       for ii = 1: length(inputs)
           par2update = inputs(ii);
           parName = prz2opt(ii);
           switch parName
               case "binNoise"
                   par2update = round(par2update); % make sure its an integer..
                   noiseDataLocal.binNoise = par2update;
               case "range_abs"
                   par2update = round(par2update); % make sure its an integer..
                   %noiseDataLocal.range(1) = -par2update; % update min
                   noiseDataLocal.range(1) = 0; % update min
                   noiseDataLocal.range(2) = par2update; % update max
               case "effectiveRS"
                   par2update = round(par2update); % make sure its an integer..
                   noiseDataLocal.effectiveRS = par2update; % update min                
           end
       end 
       clear ii
end


% Call revCorrNoise function
[im2,~] = revCorrNoise_v2(noiseDataLocal); % myNoise returns noise field 

%% Generate Spatial Frequency Power Spectra for Images 1 and 2
% Image 1 (base image)
% #########################################################################
% get powerspectrum data


% Because the base image will be the same across comparisons, we only need
% to make the powerspectrum once. For speed considerations, check the
% global iteration count, and only create the base image powerspectrum if
% it is the first pass through. Then assign it to a global variable such
% that it can pulled from memory/the global workspace all the other times..

IterCount = getGlobalIterCount;
if IterCount < 2
        try firstPass;
            % if first firstPass exists, do nothing.. this means we've
            % already looped through/we're not on the first pass..
        catch
        firstPass = 1; % firstPass won't exist on the first pass, so it will get defined here initially to flag we're on the first pass
        end
else
    firstPass = 0;
end

if firstPass == 1
    [sf_im1, nps_im1] = powerspectrum(im1,[1,1],1000);
    setGlobalBaseImPSpec(nps_im1)
    firstPass = 0; % now switch firstPass to 0 so we know after this pass that we're no longer on the first pass..
else
    nps_im1 = getGlobalBaseImPspec;
end
% #########################################################################

% Image 2 (noise image)
% #########################################################################
% get powerspectrum data
[sf_im2, nps_im2] = powerspectrum(im2,[1,1],1000);
% #########################################################################

%% Compare the Power Spectra
% Calculate RMSE between the two powerspectra..
V1 = nps_im1;
V2 = nps_im2;
%RMSE = sqrt(mean((V1-V2).^2));
MSE = mean((V1-V2).^2);
err_mat(1,jj) = MSE;
end
IterCount = getGlobalIterCount;
prz2opt = getGlobalOptPars;
inputs = varargin{1:end};

meanError= mean(err_mat(1,1:end));
% if isnan(meanRMSE)
%  pause = "";
% end
iNum = getGlobalIterCount;
disp("====================================================================")
disp(strcat("Iteration Number:"," ",num2str(iNum)))
disp("Parameters Optimized:")
disp(prz2opt)
disp("Parameter Values:")
disp(num2str(inputs))

disp(strcat("Mean Square Error Between im1 and im2 Spatial Frequency Power Spectra:"," ",num2str(meanError)))

setGlobalIterCount(IterCount+1); % update the global iteration counter..
end