function xCorOptNzImgOut = weightByTallyImg(combTallyImg,xCorAdjNzFrames)

% Combine outputs from each kernel
%--------------------------------------------------------------------------
% take weighted average..
%-------------------------------------------------
% First make adjusted version of tally image where 0s are replaced with 1s
% (This avoids divide by zero errors in the following weighted mean calc)
tallyImg2=combTallyImg;

ZeroValzMask = (tallyImg2 == 0);
% Replace all values in this mask with the mean computed above..
tallyImg2(ZeroValzMask) = 1; % Replace values

% Finally divide sum of frames in xCorOptNzImg by tallyImg2 to make
% weighted mean image that incorporates info regarding the number of totalindTallyImgz
% things (kernel image values) added at each location.
xCorOptNzImg=((sum(xCorAdjNzFrames,3))./tallyImg2); 
%-------------------------------------------------

% find all locations that equal zero in the tally image
% (places that never had a tile match as its max correlation position) 
% and replace these with half the min of the other non-zero values..
%
% First find the the non-zero values in the tally image:
nonZeroValzMask = (combTallyImg ~= 0);
% Get set of values withing this mask in xCorOptNzImg
NonZeroValz=xCorOptNzImg(nonZeroValzMask);
% Vectorize and get the mean..
NonZeroValz(:);
NonZeroValzMin=min(NonZeroValz);
% Now make mask of 0s in tally img..
ZeroValzMask = (combTallyImg == 0);
% Replace all values in this mask with the min computed above..
xCorOptNzImg(ZeroValzMask) = (NonZeroValzMin); % Replace values

% rescale xCorOptNzImg to range from 0 to 255
xCorOptNzImg=rescale(xCorOptNzImg,0,255); % scale from 0-255 
xCorOptNzImgOut=uint8(xCorOptNzImg);

end