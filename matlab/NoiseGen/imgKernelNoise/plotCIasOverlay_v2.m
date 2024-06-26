function plotCIasOverlay_v2(biIn,ciIn,ciWt,ciThr, figDims, clrMap,alphVal,idealCI,biStr,ciStr,IdlCIStr,configStr,noiseConfigStr)
% figure size parameters..
figSizeX=figDims(1);
figSizeY=figDims(2);

% Compute weighted mean image 
biWt=1-ciWt;
wtdMeanImg= ((double(rescale(biIn,0,1)).*biWt)+(rescale(ciIn,0,1).*ciWt))./2; %compute weighted mean (NOTE: this line rescales both to 0-1 to simplify the weighting/make it more intuitive..)


% Compute thresholded images/masks
%(for ciIn)
[upperMask, lowerMask, upperMskCI,lowerMskCI] = CI_Thrshld_Percent(ciIn,ciThr);
combMask=upperMask+lowerMask;
%(for ideal ci too)
[upperMaskIdl, lowerMaskIdl, upperMskCIIdl,lowerMskCIIdl] = CI_Thrshld_Percent(idealCI,ciThr);
combMaskIdl=upperMaskIdl+lowerMaskIdl;

% make colorized verison of CI (real ci)
colorCI = applyColormapToImage(ciIn, clrMap);
% apply the combined upper/lower threshold mask (real ci)
colorCI_thr = combMask.*colorCI;

% make colorized "ideal ci" too..
ClrdIdealCI= applyColormapToImage(idealCI, clrMap);
% apply the combined upper/lower threshold mask (ideal ci)
ClrdIdealCI_thr = combMaskIdl.*ClrdIdealCI;

% note that labeloverlay() uses transparency, not opacity/alpha
%alpha = 0.5;
%thrCIBIolay = labeloverlay(biIn,combMskCIAdj,'colormap',"cool",'transparency',1-alpha);

figure('Position', [10 10 figSizeX figSizeY]); % open figure window
% plot base image..
subplot(2,3,1);
imagesc(biIn);
colormap("gray");
set(gca,'DataAspectRatio',[1 1 1]); % preserve the aspect ratio
colorbar;
title({"BI:",strcat("(",biStr,")")});


% plot CI
subplot(2,3,2);
imagesc(ciIn);
colormap("gray");
set(gca,'DataAspectRatio',[1 1 1]); % preserve the aspect ratio
colorbar;
title({"Real CI:",strcat("(",ciStr,")")});

% plot ideal CI
subplot(2,3,3);
imagesc(idealCI);
colormap("gray");
set(gca,'DataAspectRatio',[1 1 1]); % preserve the aspect ratio
colorbar;
title({"Ideal CI:",strcat("(",IdlCIStr,")")});


% plot CI/BI Overlay..
subplot(2,3,4);
imagesc(wtdMeanImg);
colormap("gray");
set(gca,'DataAspectRatio',[1 1 1]); % preserve the aspect ratio
colorbar;
titleStr=strcat("(",num2str(ciWt),"/",num2str(biWt),")");
title({"RealCI/BI Overlay", strcat("weighted mean: ",titleStr)})

% plot BI with thresholded CI Overlay..
subplot(2,3,5);
% Display the grayscale base image
%alphVal=0.5; % specify alpha/transparency value..
imshow(biIn);
hold on;
% Display the colorized/thresholded CI image on top
h = imshow(colorCI_thr);
% Adjust the transparency of the color image
% Set alpha to a value between 0 and 1, where 0 is fully transparent and 1 is fully opaque
alpha(h, alphVal);
set(gca,'DataAspectRatio',[1 1 1]); % preserve the aspect ratio
colorbar;
titleStr=strcat("(upper and lower ",num2str(ciThr),"th percentile)");
title({'BI Overlaid with Thresholded CI',titleStr})
hold off;

% plot BI with "ideal" CI as colorized overlay..
subplot(2,3,6);
% Display the grayscale base image
%alphVal=0.5; % specify alpha/transparency value..
imshow(biIn);
hold on;
% Display the colorized/thresholded CI image on tClrdIdealCI_throp
h = imshow(ClrdIdealCI_thr);
% Adjust the transparency of the color image
% Set alpha to a value between 0 and 1, where 0 is fully transparent and 1 is fully opaque
alpha(h, alphVal);
set(gca,'DataAspectRatio',[1 1 1]); % preserve the aspect ratio
colorbar;
titleStr=strcat("(upper and lower ",num2str(ciThr),"th percentile)");
title({'BI Overlaid with Thresholded Ideal CI',titleStr})
hold off;

% Add title for entire panel.
sgt = sgtitle({strcat(configStr,": ",noiseConfigStr),""},'Color','red');
sgt.FontSize = 20;

end