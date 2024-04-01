
%** RevCorrSimulations2.m Snippet for making gabor with a mask...
% Ethan cut and pasted this snippet out of RevCorrSimulations2.m
% which was written by Gena for the first "simulation" paper where we
% simulated reverse correlation experiments using a 2 alternative forced
% choice task descriminating between right vs. left angled gabors.
% Ethan cut this portion out for Lauren. It simply creates a circular mask
% for the central region of the image containing the gabor. It also creates
% the associated base images. "MyMask" is the output variable of interest 
% (ie the mask..) Parameters controlling the creation of the mask/gabor 
% images are located below under "Gabor functions".

% Note: this script/code depends on also having the following functions on
% your path:
%           makeGabor.m


% --- Gabor functions --- %
gabOpts.res      = 128;   % Resolution (of base and noise)
gabOpts.tilt     = 15;    % Tilt degrees
gabOpts.contrast = 1;     % Contrast
gabOpts.aspRatio = 1;     % Aspect ration
mySF             = [4];
gabOpts.SF       = [4];     % Spatial freq % 8
gabOpts.Deg      = 1;     % Size of grating (absolute)
gabOpts.sigSize  = .2;    % Size of gaussian
MakeMask         = 1;     % Make mask
% --- Generate Base Stimuli & Mask --- %
noiseStr = [];
for iSF = 1:numel(mySF)
    %gabOpts.SF = mySF(iSF);
    BaseStim{iSF,1} = makeGabor(gabOpts.tilt,     gabOpts.contrast, gabOpts.res, gabOpts.aspRatio, gabOpts.SF, gabOpts.sigSize, gabOpts.Deg, 0);
    BaseStim{iSF,2} = makeGabor(360-gabOpts.tilt, gabOpts.contrast, gabOpts.res, gabOpts.aspRatio, gabOpts.SF, gabOpts.sigSize, gabOpts.Deg, 0);
    if MakeMask
        [~,~,MyMask] = makeGabor(gabOpts.tilt, gabOpts.contrast, gabOpts.res, gabOpts.aspRatio, gabOpts.SF, gabOpts.sigSize, gabOpts.Deg, 1);
    end

    % -- String to load/save noise (this will be long but descriptive) --- %
    %EJD COMMENTED BC this is not required for making the mask...
    %myStruct{1} = noiseData;
    %myStruct{2} = gabOpts;
    %noiseStr{iSF}    = strcat('Ntrial',num2str(numTrials),structStringify(myStruct));
end



% MASK THE BASE IMAGE
BaseImg=BaseStim{1,1}; % select one of the two base images
MaskedBaseImg=BaseImg.*MyMask; % mask by multiplying by the mask image

% PLOT THE MASK
figure();
imagesc(MyMask);
colormap("gray");
axis image

% PLOT THE BASE IMAGE ALONE
figure();
imagesc(BaseStim{1,1});
colormap("gray");
axis image

% PLOT THE MASKED BASE IMAGE
figure();
imagesc(MaskedBaseImg);
colormap("gray");
axis image
