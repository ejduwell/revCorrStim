function [baseIm, noiseData] = noiseBaseImDescFile(varargin)
% To initialize noiseData and baseIm input 1. Otherwise, input 0.
% This means you are checking the existing baseIm and noiseData
% values saved in the global variables. If no input is provided,
% the default is to set input to 1 and initialize variables..

if length(varargin)>0
set = varargin{1};
else
set = 1; % default value is 1.. ie default is to define the variables, not call the globals..
end

% If "set" = 1, set/define the output variables..
if set == 1
%############################### Base Image ###############################
%##########################################################################
% Specify Path to Base Image..
% Note: baseIm is set as a global var here to be accessible across functions without calling/running everything else in this file...

%Genna's Gabor Stimuli Base images
%baseIm = "/Users/eduwell/Library/CloudStorage/OneDrive-mcw.edu/ReverseCorrelation/NoiseGen/standard_samples/GennaBaseImage1.jpg";
%baseIm = "/Users/eduwell/Library/CloudStorage/OneDrive-mcw.edu/ReverseCorrelation/NoiseGen/standard_samples/GennaBaseImage1.jpg";

% Standard whitenoise samples for validation testing..
%baseIm = "/Users/eduwell/Library/CloudStorage/OneDrive-mcw.edu/ReverseCorrelation/NoiseGen/standard_samples/whitenoise_r255_rs1_nobin.jpg";
%baseIm = "/Users/eduwell/Library/CloudStorage/OneDrive-mcw.edu/ReverseCorrelation/NoiseGen/standard_samples/whitenoise_r255_rs0pt5_nobin.jpg";
%baseIm = "/Users/eduwell/Library/CloudStorage/OneDrive-mcw.edu/ReverseCorrelation/NoiseGen/standard_samples/whitenoise_r255_rs0pt1_nobin.jpg";
%baseIm = "/Users/eduwell/Library/CloudStorage/OneDrive-mcw.edu/ReverseCorrelation/NoiseGen/standard_samples/whitenoise_r255_rs0pt05_nobin.jpg";
%baseIm = "/Users/eduwell/Library/CloudStorage/OneDrive-mcw.edu/ReverseCorrelation/NoiseGen/standard_samples/whitenoise_r255_rs0pt01_nobin.jpg";

%baseIm = "/Users/eduwell/Library/CloudStorage/SynologyDrive-RevCorr_Shared/NoiseGen/standard_samples/rand_gauss_liu5_opt_64_256_1024.jpg";
baseIm = "/Users/eduwell/Library/CloudStorage/SynologyDrive-RevCorr_Shared/NoiseGen/standard_samples/liu5_64_256_1024_128by128.jpg";

%##########################################################################

%############################## Noise Image ###############################
%##########################################################################
% The input to the function is a single struct with fields as options
% 
% =========================== Required fields =============================

% Specify the number of samples of noise at used to compute the error signal 
% at each parameter setting 

%nNoizeSamplzPerEval = 10;
nNoizeSamplzPerEval = 10;

% noiseData.NoiseType        - type of noise: 'white','gauss','grating','gabor'
noiseData.NoiseType = 'gabor';

% noiseData.imgSize          - dimensions of noise field ([xDim, yDim]), not required if base image is provided (will use size of base image)
noiseData.imgSize   = "match";      % Size of noise image (enter 2 values if x & y dims are different, ex: [256,512] )
                                        % Note: if noiseData.imgSize = "match",
                                        % then the size of the noise image will
                                        % be sized to match the base image..
%noiseData.imgSize   = 480; 
% ============================ Optional fields ============================
if noiseData.NoiseType == 'white'
%     noiseData.binNoise = 25;
    noiseData.range_abs = 255;
    noiseData.range =[0, noiseData.range_abs];            %- noise range [min, max] (for 'guass' & 'white')
    noiseData.effectiveRS = 0.0159;

    %optPars = ["binNoise", "range_abs", "effectiveRS"];
    optPars = ["effectiveRS"];
    
    setGlobalOptPars(optPars)
    %pars2opt_start = [noiseData.binNoise, noiseData.range_abs, noiseData.range_abs];
    
    rng("shuffle");
    randLL = 0.005;
    randUL = 1;
    randStart = randLL + (randUL-randLL) .* rand(1,1);
    %randStart = randi([1 50],1,1);
    
    pars2opt_start = [randStart];
    %pars2opt_start = [noiseData.effectiveRS];
    
    %pars2opt_lb = [2, 1, 1];
    %pars2opt_ub = [255, 255, 50];
    pars2opt_lb = [0.05];
    pars2opt_ub = [1];
    %pars2opt_stpSiz = [1, 1, 1];
    pars2opt_stpSiz = [1];
end

% noiseData.mu               - noise mean (for 'gauss' option)
%noiseData.mu = 0;

% noiseData.std              - noise standard deviation (for 'gauss' option)
%noiseData.std = 1;

% noiseData.range            - noise range [min, max] (for 'guass' & 'white')

% noiseData.baseImg          - Base image. If present, function will return BaseImg + noise

% noiseData.normBase         - normalize base

% noiseData.normNoise        - normalize noise

% noiseData.normBaseNoise    - normalize base + noise image

% noiseData.binNoise         - enter number of bins

% noiseData.weightedMean     - Merge base+noise using weighted mean instead of sum (default). [Input is weight 0-1]

% noiseData.smoothBase       - Smooth base image [input is sigma value]      

% noiseData.smoothNoise      - Smooth noise image [input is sigma value]

% noiseData.smoothBN         - Smooth noise + base [input is sigma value]

% noiseData.newFig           - open new figure, 1 default or 0

% noiseData.cMap             - Add colormap for displaying images (under construction)

% noiseData.pltBase          - plot base image

% noiseData.pltBaseHist      - plot base image histogram

% noiseData.pltNoise         - plot noise

% noiseData.pltNoiseHist     - plot noise histogram

% noiseData.pltBaseNoise     - plot base + noise

% noiseData.pltBaseNoiseHist - plot base + noise histogram 

% Gaussian-Blob-Specific-Pars (noiseData.NoiseType = 'gabor'):


% White-Noise-Specific-Pars (noiseData.NoiseType = 'white'):

% Specify which Gaussian-Blob-Specific-Pars you want to have optimized..
if noiseData.NoiseType == 'gabor'
    noiseData.blobs2.sig1 = 64;
    noiseData.blobs2.sig2 = 256;
    noiseData.blobs2.sig3 = 1024;
    
%     noiseData.blobs.sig1 = 64;
%     noiseData.blobs.sig2 = 256;
%     noiseData.blobs.sig3 = 1024;
    %noiseData.blobs.gaussNum1 = 1;
    %noiseData.blobs.gaussNum2 = 5;
    %noiseData.blobs.gaussNum3 = 25;
    %noiseData.blobs.n_gauss1 = 2500;
    %noiseData.blobs.n_gauss2 = 156;
    %noiseData.blobs.n_gauss3 = 10;
    %noiseData.blobs.grn = 480;
    %options = optimset('sig1',noiseData.blobs.sig1,'sig2',noiseData.blobs.sig2,'sig3',noiseData.blobs.sig3,'gaussNum1',noiseData.blobs.gaussNum1,'gaussNum2',noiseData.blobs.gaussNum2 ,'gaussNum3',noiseData.blobs.gaussNum3,'n_gauss1',noiseData.blobs.n_gauss1 ,'n_gauss2',noiseData.blobs.n_gauss2,'n_gauss3',noiseData.blobs.n_gauss3,'grn',noiseData.blobs.grn);
    
    %optimset('TolX',1e-6,'PlotFcns',@rand_gauss_liu4)
    %options = optimset(@rand_gauss_liu4);

    %optPars = ["sig1", "sig2", "sig3", "gaussNum1", "gaussNum2", "gaussNum3", "n_gauss1", "n_gauss2", "n_gauss3", "grn"];
    optPars = ["sig1", "sig2", "sig3"];
    setGlobalOptPars(optPars)
    
%     pars2opt_start = [noiseData.blobs.sig1, noiseData.blobs.sig2, noiseData.blobs.sig3, noiseData.blobs.gaussNum1, noiseData.blobs.gaussNum2, noiseData.blobs.gaussNum3, noiseData.blobs.n_gauss1, noiseData.blobs.n_gauss2, noiseData.blobs.n_gauss3, noiseData.blobs.grn];
%     pars2opt_lb = [1, 1, 1, 1, 1, 1, 1, 1, 1, 25];
%     pars2opt_ub = [100, 500, 2000, 5, 10, 50, 5000, 300, 50, 480];
%     pars2opt_stpSiz = [5, 10, 25, 1, 1, 2, 25, 10, 2, 10];

    pars2opt_start = [noiseData.blobs2.sig1, noiseData.blobs2.sig2, noiseData.blobs2.sig3];
    pars2opt_lb = [1, 1, 1];
    pars2opt_ub = [100, 500, 2000];
    pars2opt_stpSiz = [5, 10, 25];

end


% ========================== Order of Operations ==========================
%                         (In Creating Noise Image)
% (1) Resize    (noise; base)
% (2) Smooth    (noise; base)
% (3) Normalize (noise; base)
% (4) Bin       (noise)
% (5) Merge base & noise into combined image
% (6) Resize    (combined image)
% (7) Smooth    (combined image)
% (8) Normalize (combined image)
% =========================================================================

%##########################################################################
% ======================== Save pars2opt_start to noiseData output ========================
noiseData.pars2opt_start = pars2opt_start;
noiseData.pars2opt_lb = pars2opt_lb;
noiseData.pars2opt_ub = pars2opt_ub;
noiseData.pars2opt_stpSiz = pars2opt_stpSiz;
% =================== Read in and save baseIm to output ===================
baseIm = imread(baseIm);
if length(size(baseIm))>2 % check if it has a color channel..
baseIm = rgb2gray(baseIm); % convert to greyscale..
end

% ========== Save base image and noise data to global variables ===========
setGlobalBaseIm(baseIm);
setGlobalNoiseIm(noiseData)
setGlobalNumSamplz(nNoizeSamplzPerEval);
setGlobalIterCount(1); % initialize the GlobalIterCount (keeps track of the number of iterations through the fminsearch loop) as 1.

% If "set" is set to zero, return the values already saved in the global
% variables..
% NOTE: running with this option requires that this script has already been
% called previously with "set" set to 1. Otherwise there will be no
% variables to call and you will get an error..
elseif set == 0
    baseIm = getGlobalBaseIm;
    noiseData = getGlobalNoiseIm;
end

%##########################################################################

end