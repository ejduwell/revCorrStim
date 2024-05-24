function kWindowOut = rotateKrnlWindow2(kCoords, imgIn)

% Calculate centers
imgCenter = round(size(imgIn) / 2);
windowCenter = round([(kCoords(2) + kCoords(4)) / 2, (kCoords(1) + kCoords(3)) / 2]);

% Calculate shifts
vShftVal = imgCenter(1) - windowCenter(1);
hShftVal = imgCenter(2) - windowCenter(2);

% Shift the image
imgIn = circshift(imgIn, [vShftVal, hShftVal]);

% Calculate centered window coordinates
ctrWindowH = abs(kCoords(2) - kCoords(4));
ctrWindowW = abs(kCoords(1) - kCoords(3));
ctrWindow = round([imgCenter(2) - ctrWindowW / 2, imgCenter(1) - ctrWindowH / 2, imgCenter(2) + ctrWindowW / 2, imgCenter(1) + ctrWindowH / 2]);

% Rotate the image
ranAngle = rand() * 360;
imgInRotated = imrotate(imgIn, ranAngle, 'nearest', 'crop');

% Extract the window
kWindowOut = imgInRotated(ctrWindow(2):ctrWindow(4), ctrWindow(1):ctrWindow(3));

end
