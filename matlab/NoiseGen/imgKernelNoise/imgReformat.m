
%% Parameters

% images
%imgs=["107_08.jpg","107_03.jpg"];
imgs=["031_08.jpg","031_03.jpg"];
%imgs=["right.png","left.png"];


%output image dimensions..
sizeOut=[256,256];

%output image name tag
outTag="_rfmt";

for ii = 1: length(imgs)
%% Read in image
im=imread(imgs(ii));


%% Detect if this is a color image, and if so, convert to grayscale
Dimz=size(im);
nDimz=length(Dimz); % check the number of dimensions/channels..

%if nDimz is > 2 this is a color image. convert to gray.
if nDimz>2
im=rgb2gray(im);
end

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
sqrIm=imresize(sqrIm,sizeOut);

%% Save image

% build output file name
[fPath,fName,fExt]=fileparts(imgs(ii));
outName=strcat(fName,outTag,fExt);

% save file
imwrite(sqrIm,outName);

end