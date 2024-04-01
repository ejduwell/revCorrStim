function imOut = fftMaskFilter(imIn,fSpaceRad,filter)


%reading in the image
%A=imread('coins.png');
A=imIn;

% run 2dfft..
D = fft2(double(A));

%applying fourier transform shift to zero postion and fouire transform for 2D
D=fftshift(D);

%Create Mask
imgWidth = size(D,2);  % Width of the output image
imgHeight = size(D,1);  % Height of the output image

%circleRadius=10;
circleRadius=fSpaceRad;

binaryMask = createCircleMask(imgWidth, imgHeight, circleRadius);
binaryMaskInv=imcomplement(binaryMask);

% apply masks to frequency space image..
D_LF=binaryMask.*D; % pass low frquencies/block high frequencies
D_HF=binaryMaskInv.*D; % pass high frquencies/block low frequencies


% Inverse fftshift..
D_LF=ifftshift(D_LF);
D_HF=ifftshift(D_HF);


% return to normal domain space  
D_LF=uint8(ifft2(D_LF));
D_HF=uint8(ifft2(D_HF));


D_LFr=real(D_LF);
D_LFi=imag(D_LF);
D_HFr=real(D_HF);
D_HFi=imag(D_HF);


% figure;
% imshow(D_LFr);
% 
% figure;
% imshow(D_HFr);

% D_LFr_rescale=uint8(rescale(D_LFr,0,255));
% D_HFr_rescale=uint8(rescale(D_HFr,0,255));

% figure;
% imshow(D_LFr_rescale);
% 
% figure
% imshow(D_HFr_rescale);

if filter == "hf"
imOut=D_HFr;
end

if filter == "lf"
imOut=D_LFr;
end

if filter == "both"
imOut=cell(1,2);
imOut{1,1}=D_HFr;
imOut{1,D_HFr2}=D_LFr;
end


end