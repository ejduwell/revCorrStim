%% Power spectrum playground %%

%% Load Image
im = imread('ngc6543a.jpg');
im = rgb2gray(im);
im = imresize(im,[256,256]);
%im = mynoise;
figure;
imshow(im);

%%

F=fft2(im);
F=fftshift(F);
p=(abs(F).^2);
p=log10(p);
fs=1;
f=fs/2*logspace(-10,10,size(im,1));
figure;
plot(f,p);
xlabel('frequency(Hz)');
ylabel('power');
title('{\bf LogLogPlot}');


%%

 F=fftshift(fft2(im));   
 figure
 subplot(1,2,1)
 plot(abs(F))        

 IF=ifft2(fftshift(F));
 subplot(1,2,2)
 imshow(IF)

 %%

fftImage         = fft2(im);
fftImageRealPart = real(fftImage);
[rows, columns]  = size(fftImageRealPart);
lineNumber       = floor(rows/2); % Middle row
oneLine          = fftImageRealPart(lineNumber, :);
figure
plot(oneLine, 'r-', 'LineWidth', 2);
grid on;

%%
  mag=zeros(100,100);  % magnitude (all zeros for now)
  ph=zeros(100,100);   % phase (all zeros for now)

  mag(2,3)=1;          % 1 cycle in x, 2 cycles in y

  y=mag.*exp(1i*ph);   % build fft combining mag and phase
  %x=real(ifft2(y));    % inverse fft (then take real part)
  x=real(ifft2(im));    % inverse fft (then take real part)
%   figure
%   imagesc(x)           % plot

  figure,hold on
  plot(x(:,1));
  plot(x(1,:),'r');
  legend({'x axis','y axis'});

  %%
    y = fft2(im);
    clim=quantile(abs(y(:)),[.01 .99]);
    figure
    subplot(1,2,1)
    q = fftshift(abs(y));
    imagesc(fftshift(abs(y)),clim);colormap gray
    title('magnitude');
    clim=quantile(angle(y(:)),[.01 .99]);
    subplot(1,2,2)
    imagesc(fftshift(angle(y)),clim);colormap gray
    title('phase')

    %%
%     D = load('matlab.mat');
% mat = D.mat; mat = im;

t = 0:numel(im)-1;                                         % Time Vector (Units: Time Units)
Ts = mean(diff(t));                                         % Sampling Time Interval
Fs = 1/Ts;                                                  % Sampling Frequency
Fn = Fs/2;                                                  % Nyquist Frequency
L = numel(t);                                               % Signal Length
FTmat = fft(im)/L;                                         % Fourier Transform
Fv1 = linspace(0, 1, fix(L/2)+1)*Fn;                       % Frequency Vector - One-Sided Fourier Transform (Units: Cycles/Time Unit)
Iv = 1:numel(Fv1);                                           % Index Vector
figure
plot(Fv1, abs(FTmat(Iv))*2)
grid
title('One-sided Fourier Transform')
xlabel('Frequency')
ylabel('Amplitude')
Fv2 = linspace(-Fn, +Fn, L);                               % Frequency Vector - Two-Sided Fourier Transform (Units: Cycles/Time Unit)

% figure
% plot(Fv2, fftshift(abs(FTmat)))
% grid
% title('Two-Sided Fourier Transform')
% xlabel('Frequency')
% ylabel('Amplitude')

%%
figure()
[pxx,w] = periodogram(double(im));
plot(w,10*log10(pxx))
title('periodogram')

%%
[sf, nps] = powerspectrum(im);
%[sf, nps] = powerspectrum(grayImage)
figure;
plot(sf,log10(nps)); %semilogy(sf, nps);
title('powerspectrum')

%%
Fs = 1000;
t = (0:1/Fs:0.296)';
x = cos(2*pi*t*200)+0.1*randn(size(t));
xTable = timetable(seconds(t),x);

[pxx,f] = pspectrum(xTable);
figure
plot(f,pow2db(pxx))
grid on
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')
title('Default Frequency Resolution')

%%
% im=double(imread('ngc6543a.jpg'));
% im=rgb2gray(im);
% im=imresize(im,[256,256]);
% N=256;
% imf=fftshift(fft2(im));
% impf=abs(imf).^2;
% f=-N/2:N/2-1;
% imagesc(f,f,log10(impf)), axis xy
% Pf=rotavg(impf);
% f1=0:N/2;
% loglog(f1,Pf)

%%

% 2D FFT Demo to get average radial profile of the spectrum of an image.
fontSize = 20;

  
% Check that user has the Image Processing Toolbox installed.
hasIPT = license('test', 'image_toolbox');
if ~hasIPT
  % User does not have the toolbox installed.
  message = sprintf('Sorry, but you do not seem to have the Image Processing Toolbox.\nDo you want to try to continue anyway?');
  reply = questdlg(message, 'Toolbox missing', 'Yes', 'No', 'Yes');
  if strcmpi(reply, 'No')
    % User said No, so exit.
    return;
  end
end
% Read in a standard MATLAB gray scale demo image.
folder = fileparts(which('cameraman.tif')); % Determine where demo folder is (works with all versions).
baseFileName = 'cameraman.tif';
% Get the full filename, with path prepended.
fullFileName = fullfile(folder, baseFileName);
% Check if file exists.
if ~exist(fullFileName, 'file')
  % File doesn't exist -- didn't find it there.  Check the search path for it.
  fullFileName = baseFileName; % No path this time.
  if ~exist(fullFileName, 'file')
    % Still didn't find it.  Alert user.
    errorMessage = sprintf('Error: %s does not exist in the search path folders.', fullFileName);
    uiwait(warndlg(errorMessage));
    return;
  end
end
% Read in image.
grayImage = imread('cameraman.tif');
[rows, columns, numberOfColorChannels] = size(grayImage);
if numberOfColorChannels > 1
  grayImage = rgb2gray(grayImage);
end

%grayImage = im;
% Display original grayscale image.
subplot(2, 3, 1);
imshow(grayImage)
axis on;
title('Original Gray Scale Image', 'FontSize', fontSize)
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% Perform 2D FFTs
fftOriginal = fft2(double(grayImage));
% Move center from (1,1) to (129, 129) (the middle of the matrix).
shiftedFFT = fftshift(fftOriginal);
subplot(2, 3, 2);
scaledFFTr = 255 * mat2gray(real(shiftedFFT));
imshow(log(scaledFFTr), []);
title('Log of Real Part of Spectrum', 'FontSize', fontSize)
subplot(2, 3, 3);
scaledFFTi = mat2gray(imag(shiftedFFT));
imshow(log(scaledFFTi), []);
axis on;
title('Log of Imaginary Part of Spectrum', 'FontSize', fontSize)
% Display magnitude and phase of 2D FFTs
subplot(2, 3, 4);
shiftedFFTMagnitude = abs(shiftedFFT);
imshow(log(abs(shiftedFFTMagnitude)),[]);
axis on;
colormap gray
title('Log Magnitude of Spectrum', 'FontSize', fontSize)
% Get the average radial profile
midRow = rows/2+1
midCol = columns/2+1
maxRadius = ceil(sqrt(129^2 + 129^2))
radialProfile = zeros(maxRadius, 1);
count = zeros(maxRadius, 1);
for col = 1 : columns
  for row = 1 : rows
    radius = sqrt((row - midRow) ^ 2 + (col - midCol) ^ 2);
    thisIndex = ceil(radius) + 1;
    radialProfile(thisIndex) = radialProfile(thisIndex) + shiftedFFTMagnitude(row, col);
    count(thisIndex) = count(thisIndex) + 1;
  end
end
% Get average
radialProfile = radialProfile ./ count;
subplot(2, 3, 5:6);
plot(radialProfile, 'b-', 'LineWidth', 2);
grid on;
title('Average Radial Profile of Spectrum', 'FontSize', fontSize)

%%

figure;
imdata = im;
subplot(2,3,1)
imshow(imdata); title('Gray Image');
%Get Fourier Transform of an image
F = fft2(imdata);
% Fourier transform of an image
S = abs(F);
subplot(2,3,2);imshow(S,[]);title('Fourier transform of an image');
%get the centered spectrum
Fsh = fftshift(F);
subplot(2,3,3);imshow(abs(Fsh),[]);title('Centered fourier transform of Image')
%apply log transform
S2 = log(1+abs(Fsh));
subplot(2,3,4);imshow(S2,[]);title('log transformed Image')
%reconstruct the Image
F = ifftshift(Fsh);
f = ifft2(F);
subplot(2,3,5);imshow(f,[]),title('reconstructed Image')

%%

noiseData.NoiseType = 'grating';  % Noise type (grating)
noiseData.imgSize   = 512;        % Size of noise image

for i = 1:100

disp(i)
% --- Call Function --- %
[myNoise,~] = revCorrNoise(noiseData);

[sf, nps] = powerspectrum(myNoise,[1,1],400);
storeData(i,:) = log10(nps);

%plot(sf,log10(nps)); %semilogy(sf, nps);

end

plot(mean(storeData))