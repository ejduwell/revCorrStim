function imgOut = imgReformater(im,desiredDimz)

%% Parameters

sizeOut=desiredDimz;

%% Detect if this is a color image, and if so, convert to grayscale
Dimz=size(im);
nDimz=length(Dimz); % check the number of dimensions/channels..

%if nDimz is > 2 this is a color image. convert to gray.
if nDimz>2
im=rgb2gray(im);
Dimz=size(im);
end

%% Resize if input image dimensions do not match desired dimensions

if Dimz(1)~=desiredDimz(1) || Dimz(2)~=desiredDimz(2)

%% Cut square out of center to ensure x and y dimensions are equal

% save height/width
imH=Dimz(1);
imW=Dimz(2);

% get the images smallest spatial dimension to set square size
squareDim=min([imH,imW]); 

% calculate position of centered square
squareXpos=((imW/2)-(squareDim/2)+1);
squareYpos=((imH/2)-(squareDim/2)+1);

sqrIm=im((squareYpos:(squareYpos+squareDim-1)),(squareXpos:(squareXpos+squareDim)-1));

%% Resize Image
imgOut=imresize(sqrIm,sizeOut);

else

% if sizes match.. do nothing..
imgOut=im;

end

end