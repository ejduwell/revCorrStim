function [imNoizeComb] = revCorrAddNoise_v1(imgIn,noiseImgIn,nWt,fixIndOlayIm,nonFxnIndMask)

imWt = 1-nWt; % set texture image weight to be 1-lum weight..
%take weighted average of lum-only and tex-only images...
%imNoizeComb = uint8(imWt*double(imgIn) + nWt*double(noiseImgIn));
imNoizeComb = uint8(imWt*double(imgIn) + nWt*double(noiseImgIn));

% combine fixIndOlayIm logically with imNoizeComb such that all 0-value
% pixels in fixIndOlayIm (locations without fixation point/inducers) are
% set equal to the corresponding index values in imNoizeComb, but the
% non-0-value pixels remain whatever they are in the fixIndOlayIm
imNoizeComb=nonFxnIndMask.*imNoizeComb;
imNoizeComb=imNoizeComb+fixIndOlayIm;
pause="";
end