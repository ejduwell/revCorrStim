function imgout = revCorrMakeExmplImgs(parzIn)

%% Unpack Input Parameters

occImgR=parzIn.occImgR;
occImgL=parzIn.occImgL;
noccImgR=parzIn.noccImgR;
noccImgL=parzIn.noccImgL;
noiseDir=parzIn.noiseDir;
useNoise=parzIn.useNoise;
fix_im=parzIn.fix_im;
nWt=parzIn.nWt;

%% If noise requested, combine noise with base images

if useNoise == 1

    % Get noise images..
    strtDir=pwd; % save start path
    cd(noiseDir); % go to noise dir
    % take the top 4 images
    noiseImg1=imread('noiseSample00001.png');
    noiseImg2=imread('noiseSample00002.png');
    noiseImg3=imread('noiseSample00003.png');
    noiseImg4=imread('noiseSample00004.png');
    cd(strtDir); % go back

    % Using the fixation only image as a starting point,
    % generate an image with 0s everywhere but where the fixation point and
    % inducer lines are (ie everything with a value darker than 127..
    % everything that isn't gray..) This will be combined logically with
    % each BI+Noise image so the appearance of the fixation point/inducers
    % stays constant regardless of the noise..
    fixIndOlayIm=fix_im<127;
    fixIndOlayIm=uint8(double(fixIndOlayIm).*double(fix_im));
    nonFxnIndMask=uint8(fixIndOlayIm==0); %mask with 1s where fixIndOlayIm is 0 (no fixation marker/inducers) and 0s where fixIndOlayIm is non-zero.

    % Combine Noise with base images
    [occImgR] = revCorrAddNoise_v1(occImgR,noiseImg1,nWt,fixIndOlayIm,nonFxnIndMask);
    [occImgL] = revCorrAddNoise_v1(occImgL,noiseImg2,nWt,fixIndOlayIm,nonFxnIndMask);
    [noccImgR] = revCorrAddNoise_v1(noccImgR,noiseImg3,nWt,fixIndOlayIm,nonFxnIndMask);
    [noccImgL] = revCorrAddNoise_v1(noccImgL,noiseImg4,nWt,fixIndOlayIm,nonFxnIndMask);


end

%% Create example image array figure

occImgR=cat(3,occImgR,occImgR,occImgR);
occImgL=cat(3,occImgL,occImgL,occImgL);
noccImgR=cat(3,noccImgR,noccImgR,noccImgR);
noccImgL=cat(3,noccImgL,noccImgL,noccImgL);

f = figure('visible','off');
f.Position = [100 100 (size(occImgR,1)*2)+14 (size(occImgR,2)*2)+14];

if useNoise==1
panelTitle="Objects With Noise Overlaid";
elseif useNoise==0
panelTitle="Objects Without Noise Overlaid";
end

set(gcf,'color',[0.5,0.5,0.5]);
hold on;

%check if the tiledlayout function exists on this machine
if exist('tiledlayout')>0
tL=1;
else
tL=0;
end

%tiledlayout(2,2,"TileSpacing","tight","Padding","tight");

if tL==1
t=tiledlayout("flow","TileSpacing","tight","Padding","tight");
title(t,panelTitle,'FontSize',35,'FontWeight','bold','Color','r');
nexttile
else
  subplot(2,2,1);
  %subplot_tight(2,2,1,[0.01]); 
end

imagesc(noccImgL);
colormap("gray");
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
axis image;
ylabel("(Not Occluded)",'FontSize',30,'FontWeight','bold','Color','r');
title("LEFT",'FontSize',30,'FontWeight','bold','Color','r');

if tL==1
nexttile
else
subplot(2,2,2);
%subplot_tight(2,2,2,[0.01]);
end
imagesc(noccImgR);
colormap("gray");
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
axis image;
title("RIGHT",'FontSize',30,'FontWeight','bold','Color','r');

if tL==1
nexttile
else
subplot(2,2,3);
%subplot_tight(2,2,3,[0.01]);
end
%subplot(2,2,4);
imagesc(occImgL);
colormap("gray");
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
ylabel("(Occluded)",'FontSize',30,'FontWeight','bold','Color','r');
axis image;

if tL==1
nexttile
else
subplot(2,2,4);
%subplot_tight(2,2,4,[0.01]);
end
%subplot(2,2,3);
imagesc(occImgR);
colormap("gray");
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
axis image;

if tL==0
sgt = sgtitle(panelTitle,'Color','red','FontWeight','bold');
sgt.FontSize = 35;
end
hold off;

%% Convert figure to image
F = getframe(gcf);
[imgout, ~] = frame2im(F);

end