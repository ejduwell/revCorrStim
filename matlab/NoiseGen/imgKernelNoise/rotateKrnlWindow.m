function kWindowOut = rotateKrnlWindow(kCoords,imgIn)


% Apply shifts to position the window at center of imgIn

% figure;
% imshow(imgIn);

imgCenter=[round(size(imgIn,1)/2),round(size(imgIn,2)/2)];
windowCenter=[round((kCoords(2)+kCoords(4))/2),round((kCoords(1)+kCoords(3))/2)];

%kSubImg=(imgIn(kCoords(2):kCoords(4),kCoords(1):kCoords(3)));

% figure;
% imshow(kSubImg);

vShftVal=imgCenter(1)-windowCenter(1);
hShftVal=imgCenter(2)-windowCenter(2);

vShftOptions={'up','down'};
hShftOptions={'left','right'};

if vShftVal<0
vShftDir=vShftOptions{1};
else
vShftDir=vShftOptions{2};
end

if hShftVal<0
hShftDir=hShftOptions{1};
else
hShftDir=hShftOptions{2};
end

vShftVal=round(abs(vShftVal));
hShftVal=round(abs(hShftVal));

imgIn = shiftImage(imgIn, vShftDir, vShftVal); % apply vertical shift
imgIn = shiftImage(imgIn, hShftDir, hShftVal); % apply horizontal shift

% figure;
% imshow(imgIn);

ctrWindowH=abs(kCoords(2)-kCoords(4));
ctrWindowW=abs(kCoords(1)-kCoords(3));
ctrWindow=[round(imgCenter(2)-ctrWindowW/2),round(imgCenter(1)-ctrWindowH/2),round(imgCenter(2)+ctrWindowW/2),round(imgCenter(1)+ctrWindowH/2)];

%centerTest=imgIn(ctrWindow(2):ctrWindow(4),ctrWindow(1):ctrWindow(3));
% figure;
% imshow(centerTest);

ranAngle=(rand(1)*360);
imgInRotated=imrotate(imgIn,ranAngle,"nearest","crop");

% figure;
% imshow(imgInRotated);

kWindowOut=imgInRotated(ctrWindow(2):ctrWindow(4),ctrWindow(1):ctrWindow(3));

% figure;
% imshow(kWindowOut);

end