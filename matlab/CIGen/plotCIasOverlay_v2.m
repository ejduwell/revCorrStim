function plotCIasOverlay_v2(biIn1,biIn2,cisIn,ciThr,figDims, clrMap,alphVal,configStr,noiseConfigStr)

% figure size parameters..
figSizeX=figDims(1);
figSizeY=figDims(2);

% initialize structs for colored/thresholded cis 
colorCIs=struct;
colorCIs_thr=struct;

% Get list of fields in the input struct
cisIn_Fields=fieldnames(cisIn);

% Loop through fields and make colored/thresholded versions
% saved in struct copies initialized above
%==========================================================================
for ii = 1:size(cisIn_Fields,1)
    % Compute thresholded images/masks
    %(for ciIn)
    [upperMask, lowerMask, upperMskCI,lowerMskCI] = CI_Thrshld_Percent(cisIn.(cisIn_Fields{ii,1}),ciThr);
    combMask=upperMask+lowerMask;
    
    % make colorized verison of CI (real ci)
    colorCIs.(cisIn_Fields{ii,1}) = applyColormapToImage(cisIn.(cisIn_Fields{ii,1}), clrMap);
    % apply the combined upper/lower threshold mask (real ci)
    colorCIs_thr.(cisIn_Fields{ii,1}) = combMask.*colorCIs.(cisIn_Fields{ii,1});
end
%==========================================================================

% note that labeloverlay() uses transparency, not opacity/alpha
%alpha = 0.5;
%thrCIBIolay = labeloverlay(biIn,combMskCIAdj,'colormap',"cool",'transparency',1-alpha);

figure('Position', [0 0 figSizeX figSizeY]); % open figure window
% plot base image..
subplot(2,3,1);
% Display the grayscale base image
%alphVal=0.5; % specify alpha/transparency value..
imshow(biIn1);
hold on;
% Display the colorized/thresholded CI image on top
%h = imshow(colorCI_thr.CI_results1both);
h = imshow(colorCIs_thr.(cisIn_Fields{1,1}));
% Adjust the transparency of the color image
% Set alpha to a value between 0 and 1, where 0 is fully transparent and 1 is fully opaque
alpha(h, alphVal);
set(gca,'DataAspectRatio',[1 1 1]); % preserve the aspect ratio
titleStr=strcat("(upper and lower ",num2str(ciThr),"th percentile)");
title({'BI Overlaid with Thresholded CI1 both Occ and Nocc',titleStr});
hold off;


% plot CI
subplot(2,3,2);
% Display the grayscale base image
%alphVal=0.5; % specify alpha/transparency value..
imshow(biIn1);
hold on;
% Display the colorized/thresholded CI image on top
%h = imshow(colorCI_thr.CI_results1occ);
h = imshow(colorCIs_thr.(cisIn_Fields{3,1}));
% Adjust the transparency of the color image
% Set alpha to a value between 0 and 1, where 0 is fully transparent and 1 is fully opaque
alpha(h, alphVal);
set(gca,'DataAspectRatio',[1 1 1]); % preserve the aspect ratio
titleStr=strcat("(upper and lower ",num2str(ciThr),"th percentile)");
title({'BI Overlaid with Thresholded CI1 Occ Only',titleStr})
hold off;

% plot ideal CI
subplot(2,3,3);
% Display the grayscale base image
%alphVal=0.5; % specify alpha/transparency value..
imshow(biIn2);
hold on;
% Display the colorized/thresholded CI image on top
%h = imshow(colorCI_thr.CI_results1nocc);
h = imshow(colorCIs_thr.(cisIn_Fields{5,1}));
% Adjust the transparency of the color image
% Set alpha to a value between 0 and 1, where 0 is fully transparent and 1 is fully opaque
alpha(h, alphVal);
set(gca,'DataAspectRatio',[1 1 1]); % preserve the aspect ratio
titleStr=strcat("(upper and lower ",num2str(ciThr),"th percentile)");
title({'BI Overlaid with Thresholded CI1 Nocc Only',titleStr})
hold off;


% plot CI/BI Overlay..
subplot(2,3,4);
% Display the grayscale base image
%alphVal=0.5; % specify alpha/transparency value..
imshow(biIn1);
hold on;
% Display the colorized/thresholded CI image on top
%h = imshow(colorCI_thr.CI_results2both);
h = imshow(colorCIs_thr.(cisIn_Fields{2,1}));
% Adjust the transparency of the color image
% Set alpha to a value between 0 and 1, where 0 is fully transparent and 1 is fully opaque
alpha(h, alphVal);
set(gca,'DataAspectRatio',[1 1 1]); % preserve the aspect ratio
titleStr=strcat("(upper and lower ",num2str(ciThr),"th percentile)");
title({'BI Overlaid with Thresholded CI2 both Occ and Nocc',titleStr})
hold off;

% plot BI with thresholded CI Overlay..
subplot(2,3,5);
% Display the grayscale base image
%alphVal=0.5; % specify alpha/transparency value..
imshow(biIn1);
hold on;
% Display the colorized/thresholded CI image on top
%h = imshow(colorCI_thr.CI_results2occ);
h = imshow(colorCIs_thr.(cisIn_Fields{4,1}));
% Adjust the transparency of the color image
% Set alpha to a value between 0 and 1, where 0 is fully transparent and 1 is fully opaque
alpha(h, alphVal);
set(gca,'DataAspectRatio',[1 1 1]); % preserve the aspect ratio
titleStr=strcat("(upper and lower ",num2str(ciThr),"th percentile)");
title({'BI Overlaid with Thresholded CI2 Occ Only',titleStr})
hold off;

% plot BI with "ideal" CI as colorized overlay..
subplot(2,3,6);
% Display the grayscale base image
%alphVal=0.5; % specify alpha/transparency value..
imshow(biIn2);
hold on;
% Display the colorized/thresholded CI image on top
%h = imshow(colorCI_thr.CI_results2nocc);
h = imshow(colorCIs_thr.(cisIn_Fields{6,1}));
% Adjust the transparency of the color image
% Set alpha to a value between 0 and 1, where 0 is fully transparent and 1 is fully opaque
alpha(h, alphVal);
set(gca,'DataAspectRatio',[1 1 1]); % preserve the aspect ratio
titleStr=strcat("(upper and lower ",num2str(ciThr),"th percentile)");
title({'BI Overlaid with Thresholded CI2 Nocc Only',titleStr})
hold off;

% Split title in half if greater than threshold
titleLngth=length(char(noiseConfigStr));
lenThrld=100;
splitChar=" ";
if titleLngth>lenThrld
    hlfLen= round(titleLngth/2);
    noiseConfigStr=char(noiseConfigStr);
    [firstPart, secondPart] = splitTitle(noiseConfigStr, hlfLen, splitChar);
    % Add title for entire panel.
    sgt = sgtitle({firstPart,secondPart," "},'Color','red');
    sgt.FontSize = 20;
else
    % Add title for entire panel.
    sgt = sgtitle({noiseConfigStr,""},'Color','red');
    sgt.FontSize = 20;
end




end