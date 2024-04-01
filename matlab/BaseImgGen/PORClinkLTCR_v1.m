function [backwindow] = PORClinkLTCR_v1(mainsrc, window, backwindow,vertOcc,horzOcc,boxUL,boxUR,boxLL,boxLR,gray,gray1,gray2,trialprop,occ,CR,CR_penWidth,CR_penLum,CR_al,CR_obj,scrnWidHei,rscale_f1, rscale_f2, image_in, db_mode,normalize,norm_range,maindim,T,T_al,L,L_al,lumWt,cr_long,cr_short,cmdLineOut)
% NOTE: THE ORIENTATION CONDITIONS (ORIENT) CURRENTLY ASSUME THAT THE
% OUTPUT IMAGE WILL BE ROTATED +45 DEGREES WHEN DISPLAYED..


% %For debugging/Dev.. run this to blank out image
%Screen('FillRect', backwindow, gray);  % blank out the offscreen window
%
% PORClinkL4a.m
%
% function for Perceptual Organization Reverse Correlation experiments
% generates linking display for Luminance blocks
% Called by: PORCpractice4a.m PORCquest4a.m
%
% Variables:
%
%


% EJD NOTE: 'm'= vert here (prior to rotation)
%           'z' = horiz here (prior to rotation)            

% Written by Adam Greenberg, CMU/Psych
% Sept, 2010
%
%
%% First some basic setup
Screen('FillRect', window, gray);  % blank out the onscreen window
Screen('FillRect', backwindow, gray);  % blank out the offscreen window
% Get the screen numbers
screens = Screen('Screens');

% Draw to the external screen if avaliable
screenNumber = max(screens);

% Define black and white
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
grey = white / 2;
inc = white - grey;

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', backwindow, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
% get screen center..
scrn_center = [(scrnWidHei(1)/2),(scrnWidHei(2)/2)];

% get normalized units of width/height in expressed as % of screen dimension..
w_perc = (scrnWidHei(1)/100);
h_perc = (scrnWidHei(2)/100);

% thickness of skinny rects
skinny = round(RectWidth(boxUL)*.2); %orig
%skinny = round(RectWidth(boxUL)*0.5); %ejd test..

% a little logic to figure out which rect gets canonical color
orient = trialprop{1};
%orient = trialprop(1);
%% Common Region Parameters..
%set dims of long and short lengths for the Common Region % (NOTE: total widths/heights will be 2x these
% values..)
% cr_long = 22; % length of the "long" dimension of the Common Region.. (in units of % of screen height)
% cr_short = 10; % length of the "short" dimension of the Common Region.. (in units of % of screen height)

if orient=='m'
cr_wid = cr_short*h_perc; % width
cr_hei = cr_long*h_perc; % height
cr_wid_na = cr_long*h_perc; % width
cr_hei_na = cr_short*h_perc; % height
end
if orient=='z'
cr_wid = cr_long*h_perc; % width
cr_hei = cr_short*h_perc; % height
cr_wid_na = cr_short*h_perc; % width
cr_hei_na = cr_long*h_perc; % height
end

%% Build Shapes/Objects..
    if L_al == 1
        L_orient = orient;
    elseif L_al == 0
    if orient == 'm'
        L_orient = 'z';
    end
    if orient == 'z'
        L_orient = 'm';
    end
    end
if L_orient=='m'
    if trialprop{2}==1
        ULcolr = gray1;
        LLcolr = gray1;
        URcolr = gray2;
        LRcolr = gray2;
        LTcolr = gray1;
        RBcolr = gray2;
        LTskinnycolr = gray1;
        RBskinnycolr = gray2;
    elseif trialprop{2}==2
        ULcolr = gray2;
        LLcolr = gray2;
        URcolr = gray1;
        LRcolr = gray1;
        LTcolr = gray2;
        RBcolr = gray1;
        LTskinnycolr = gray2;
        RBskinnycolr = gray1;
    end
elseif L_orient=='z'
    if trialprop{2}==1
        ULcolr = gray1;
        LLcolr = gray2;
        URcolr = gray1;
        LRcolr = gray2;
        LTcolr = gray1;
        RBcolr = gray2;
        LTskinnycolr = gray1;
        RBskinnycolr = gray2;
    elseif trialprop{2}==2
        ULcolr = gray2;
        LLcolr = gray1;
        URcolr = gray2;
        LRcolr = gray1;
        LTcolr = gray2;
        RBcolr = gray1;
        LTskinnycolr = gray2;
        RBskinnycolr = gray1;
    end
end

%%% draw background shapes
% first, for occluded condition (occ==0)
pointListUL = [boxUL(1),boxUL(4);boxUL(1),boxUL(2)+RectHeight(boxUL)/2;...
    boxUL(1)+RectWidth(boxUL)/2,boxUL(2);boxUL(3),boxUL(2);boxUL(3),boxUL(4)];
pointListUR = [boxUR(1),boxUR(4);boxUR(1),boxUR(2);boxUR(1)+RectWidth(boxUR)/2,boxUR(2);...
    boxUR(3),boxUR(2)+RectHeight(boxUR)/2;boxUR(3),boxUR(4)];
pointListLL = [boxLL(1),boxLL(2)+RectHeight(boxLL)/2;boxLL(1),boxLL(2);...
    boxLL(3),boxLL(2);boxLL(3),boxLL(4);boxLL(1)+RectWidth(boxLL)/2,boxLL(4)];
pointListLR = [boxLR(1),boxLR(4);boxLR(1),boxLR(2);boxLR(3),boxLR(2);...
    boxLR(3),boxLR(2)+RectHeight(boxLR)/2;boxLR(1)+RectWidth(boxLR)/2,boxLR(4)];

% now, for non-occluded condition (occ==1)
if orient=='m'
    LTrect = [boxUL(1) boxUL(4) boxUL(3) boxLL(2)];
    RBrect = [boxUR(1) boxUR(4) boxUR(3) boxLR(2)];
    LTskinny = [boxUL(3) boxUL(2) boxUL(3)+skinny boxLL(4)]; 
    RBskinny = [boxUR(1)-skinny boxUR(2) boxUR(1) boxLR(4)];
    cr_trans = ((abs(LTrect(1)-scrn_center(1))+abs(LTrect(3)-scrn_center(1)))/2)-1.25*h_perc; % translation from center along x axis (L/R) or y axis (Up/down).. (distance of the object offset from center)

elseif orient=='z'
    LTrect = [boxUL(3) boxUL(2) boxUR(1) boxUL(4)];
    RBrect = [boxLL(3) boxLL(2) boxLR(1) boxLL(4)];
    LTskinny = [boxUL(1) boxUL(4) boxUR(3) boxUL(4)+skinny]; 
    RBskinny = [boxLL(1) boxLL(2)-skinny boxLR(3) boxLL(2)];
    cr_trans = ((abs(LTrect(2)-scrn_center(2))+abs(LTrect(4)-scrn_center(2)))/2)-1.25*h_perc; % translation from center along x axis (L/R) or y axis (Up/down).. (distance of the object offset from center)
    
end

%%% draw appropriate shapes, in correct order
% draw surface patches
Screen('FillPoly', backwindow, ULcolr, pointListUL);
Screen('FillPoly', backwindow, URcolr, pointListUR);
Screen('FillPoly', backwindow, LRcolr, pointListLR);
Screen('FillPoly', backwindow, LLcolr, pointListLL);

% draw occluders on offscreen window
% Screen('FillRect', backwindow, gray, vertOcc);
% Screen('FillRect', backwindow, gray, horzOcc);

if occ==0 % occluded condition
    % do nothing
%     % draw surface patches
%     Screen('FillPoly', backwindow, ULcolr, pointListUL);
%     Screen('FillPoly', backwindow, URcolr, pointListUR);
%     Screen('FillPoly', backwindow, LRcolr, pointListLR);
%     Screen('FillPoly', backwindow, LLcolr, pointListLL);
%     % draw occluders on offscreen window
%     Screen('FillRect', backwindow, gray, vertOcc);
%     Screen('FillRect', backwindow, gray, horzOcc);
elseif occ==1 % non-occluded condition
    %draw non-occluded portions to offscreen window
    Screen('FillRect', backwindow, LTcolr, LTrect);
    Screen('FillRect', backwindow, RBcolr, RBrect);
    Screen('FillRect', backwindow, LTskinnycolr, LTskinny);
    Screen('FillRect', backwindow, RBskinnycolr, RBskinny);
end;


% test_im=Screen('GetImage', backwindow);
% pause = "";

% if luminance modulation is requested.. save a copy of the lum only
% objects.. configured/printed on screen..
if L == 1
  
    lumOnly_im =Screen('GetImage', backwindow); %capture image of backwindow as it stands now...
 
%     if normalize == 1
%         lumOnly_im = rescale(lumOnly_im,norm_range(1),norm_range(2)); 
%         image_tex = Screen('MakeTexture', window, lumOnly_im);
%         Screen('DrawTexture', backwindow, image_tex); 
%         lumOnly_im =Screen('GetImage', backwindow); %capture image of backwindow as it stands now...
%     end
end

if CR == 1 
%...TEST ZONE...
if CR_al == 1 % common region boxes align with objects
if orient =="z"

% Reminder: framerects rectangle coordinate values are: [Left, Top, Right, Bottom]
if CR_obj == 1
crBoxL = (scrn_center(1)-cr_wid);
crBoxT = (scrn_center(2)-cr_hei)+cr_trans;
crBoxR = (scrn_center(1)+cr_wid);
crBoxB = (scrn_center(2)+cr_hei)+cr_trans;
end
if CR_obj == 2
crBoxL = (scrn_center(1)-cr_wid);
crBoxT = (scrn_center(2)-cr_hei)-cr_trans;
crBoxR = (scrn_center(1)+cr_wid);
crBoxB = (scrn_center(2)+cr_hei)-cr_trans;
end


elseif orient=="m"

% Reminder: framerects rectagle coordinate values are: [Left, Top, Right, Bottom]
if CR_obj == 1
crBoxL = (scrn_center(1)-cr_wid)+cr_trans;
crBoxT = (scrn_center(2)-cr_hei);
crBoxR = (scrn_center(1)+cr_wid)+cr_trans;
crBoxB = (scrn_center(2)+cr_hei);
end
if CR_obj == 2
crBoxL = (scrn_center(1)-cr_wid)-cr_trans;
crBoxT = (scrn_center(2)-cr_hei);
crBoxR = (scrn_center(1)+cr_wid)-cr_trans;
crBoxB = (scrn_center(2)+cr_hei);
end

end
end

if CR_al == 0 % common region boxes DON'T align with objects
    if orient =="m"

% Reminder: framerects rectagle coordinate values are: [Left, Top, Right, Bottom]
if CR_obj == 1
crBoxL = (scrn_center(1)-cr_wid_na);
crBoxT = (scrn_center(2)-cr_hei_na)+cr_trans;
crBoxR = (scrn_center(1)+cr_wid_na);
crBoxB = (scrn_center(2)+cr_hei_na)+cr_trans;
end
if CR_obj == 2
crBoxL = (scrn_center(1)-cr_wid_na);
crBoxT = (scrn_center(2)-cr_hei_na)-cr_trans;
crBoxR = (scrn_center(1)+cr_wid_na);
crBoxB = (scrn_center(2)+cr_hei_na)-cr_trans;
end

elseif orient=="z"

% Reminder: framerects rectagle coordinate values are: [Left, Top, Right, Bottom]
if CR_obj == 1
crBoxL = (scrn_center(1)-cr_wid_na)+cr_trans;
crBoxT = (scrn_center(2)-cr_hei_na);
crBoxR = (scrn_center(1)+cr_wid_na)+cr_trans;
crBoxB = (scrn_center(2)+cr_hei_na);
end
if CR_obj == 2
crBoxL = (scrn_center(1)-cr_wid_na)-cr_trans;
crBoxT = (scrn_center(2)-cr_hei_na);
crBoxR = (scrn_center(1)+cr_wid_na)-cr_trans;
crBoxB = (scrn_center(2)+cr_hei_na);
end

    end
end
penWidth = CR_penWidth;

% % Non-Occluded Shapes
% Screen('FillRect', backwindow, LTcolr, LTrect);
% Screen('FillRect', backwindow, RBcolr, RBrect);
% Screen('FillRect', backwindow, LTskinnycolr, LTskinny);
% Screen('FillRect', backwindow, RBskinnycolr, RBskinny);
% 
% % Occluded Shapes
% Screen('FillPoly', backwindow, ULcolr, pointListUL);
% Screen('FillPoly', backwindow, URcolr, pointListUR);
% Screen('FillPoly', backwindow, LRcolr, pointListLR);
% Screen('FillPoly', backwindow, LLcolr, pointListLL);


% Reminder: framerects rectagle coordinate values are: [Left, Top, Right, Bottom]
CRbox = [crBoxL, crBoxT, crBoxR, crBoxB];

%Screen('FrameRect', backwindow, 0, CRbox,penWidth);

%Screen('FrameRect', backwindow, RBskinnycolr, CRbox2,penWidth);
%...TEST ZONE...
end
%% If Running Texture Version: Load in the image and scale it to the 2 different scalings
if T == 1
theImage = imread(image_in); % load image in from file..
if length(size(theImage)) > 2
theImage = rgb2gray(theImage); % convert the image to be grey scale if it isn't already...
end
theImage = cast(theImage,'uint8'); % For ease of coding though, make it uint8 again 
theImage = cat(3, theImage, theImage, theImage); % and concatenate 3 together so it has the proper "RGB" dimensions for later operations..

% Apply the 2 size scalings
inputIm=theImage;
theImage = imresize(inputIm, rscale_f1);
theImage_2 = imresize(inputIm, rscale_f2);

% EJD INSERTED FOR DEBUGGING
%===================================
% sizestr1=strcat("theImage dimz before crop: ", mat2str(size(theImage)));
% sizestr2=strcat("theImage_2 dimz before crop: ", mat2str(size(theImage_2)));
% sizestr3=strcat("inputImCpy dimz before crop: ", mat2str(size(inputIm)));
%===================================

% Get the dimensions of the offscreen window as an image
Screen('FillRect', backwindow, grey);  % black out the offscreen window

im_4_WinDims = Screen('GetImage', backwindow, mainsrc); %Save the window as an image
[s1, s2, s3] = size(im_4_WinDims); %save it's dimensions..

% crop the central square out of the two images corresponding the offscreen
% window dimensions
targetSize = [s1 s2];
r = centerCropWindow2d(size(theImage),targetSize); % QUESTION: Do we want to randomize where in the image the 54X54 comes from? Or keep them both centered on the same region? (as we currently are)
r2 = centerCropWindow2d(size(theImage_2),targetSize);

% EJD INSERTED FOR DEBUGGING
%===================================
% r3 = centerCropWindow2d(size(inputIm),targetSize);
% inputIm = imcrop(inputIm,r3);
%===================================

theImage = imcrop(theImage,r);
theImage_2 = imcrop(theImage_2,r2);

% EJD INSERTED FOR DEBUGGING
%===================================
% prdctdDiff=(rscale_f1/rscale_f2);
% disp(" ");
% disp(sizestr3);
% disp(sizestr1);
% disp(sizestr2);
% disp(strcat("rscale_f1: ",num2str(rscale_f1)));
% disp(strcat("rscale_f2: ",num2str(rscale_f2)));
% disp(strcat("If assumptions are correct, image scaling diff should be: ",num2str(prdctdDiff)));
% 
% figure;
% hold on
% imshow(theImage);
% title("ScaledIm1");
% axis image;
% truesize;
% hold off
% 
% figure;
% hold on
% imshow(theImage_2);
% title("ScaledIm2");
% axis image;
% truesize;
% hold off
% 
% figure;
% hold on
% imshow(inputIm);
% title("NoScale");
% axis image;
% truesize;
% hold off
% 
% pause="";
% close all
%===================================

end

%% (If Texture..) Determine which shape gets theImage and which gets theImage2


if T == 1
orient = trialprop{1};
% check if texture is supposed to be aligned or not.. adjust accordingly..
if T_al == 0
    if orient == 'm'
        tex_orient = 'z';
    end
    if orient == 'z'
        tex_orient = 'm';
    end
elseif T_al == 1
    tex_orient = orient;
end
%if orient=='z' %orig
if tex_orient=='m' %ejd switched with z to try to flip..
    if trialprop{2}==1
        UL_im = theImage;
        LL_im = theImage;
        UR_im = theImage_2;
        LR_im = theImage_2;
        LT_im = theImage;
        RB_im= theImage_2;
        LTsk_im = theImage;
        RBsk_im= theImage_2;
    elseif trialprop{2}==2
        UL_im = theImage_2;
        LL_im = theImage_2;
        UR_im = theImage;
        LR_im = theImage;
        LT_im = theImage_2;
        RB_im = theImage;
        LTsk_im = theImage_2;
        RBsk_im = theImage;
    end
%elseif orient=='m' %orig
elseif tex_orient=='z'  %ejd switched with m to try to flip..
    if trialprop{2}==1
        UL_im = theImage;
        LL_im = theImage_2;
        UR_im = theImage;
        LR_im = theImage_2;
        LT_im = theImage;
        RB_im = theImage_2;
        LTsk_im = theImage;
        RBsk_im = theImage_2;
    elseif trialprop{2}==2
        UL_im = theImage_2;
        LL_im = theImage;
        UR_im = theImage_2;
        LR_im = theImage;
        LT_im = theImage_2;
        RB_im = theImage;
        LTsk_im = theImage_2;
        RBsk_im = theImage;
    end
end
end
%% (If Texture..) Make Masks of the Shapes
if T == 1
% Begin making masks for occluded condition

Screen('FillRect', backwindow, black);  % black out the offscreen window
Screen('FillPoly', backwindow, 1, pointListUL);
blk_shapes_ul = Screen('GetImage', backwindow, mainsrc);
Screen('FillRect', backwindow, black);  % black out the offscreen window

Screen('FillRect', backwindow, black);  % black out the offscreen window
Screen('FillPoly', backwindow, 1, pointListUR);
blk_shapes_ur = Screen('GetImage', backwindow, mainsrc);
Screen('FillRect', backwindow, black);  % black out the offscreen window

Screen('FillRect', backwindow, black);  % black out the offscreen window
Screen('FillPoly', backwindow, 1, pointListLR);
blk_shapes_lr = Screen('GetImage', backwindow, mainsrc);
Screen('FillRect', backwindow, black);  % black out the offscreen window

Screen('FillRect', backwindow, black);  % black out the offscreen window
Screen('FillPoly', backwindow, 1, pointListLL);
blk_shapes_ll = Screen('GetImage', backwindow, mainsrc);
Screen('FillRect', backwindow, black);  % black out the offscreen window

Screen('FillPoly', backwindow, 1, pointListUL);
Screen('FillPoly', backwindow, 1, pointListUR);
Screen('FillPoly', backwindow, 1, pointListLR);
Screen('FillPoly', backwindow, 1, pointListLL);
blk_shapes_all = Screen('GetImage', backwindow, mainsrc);

% begin making masks for un-occluded conditions too...
Screen('FillRect', backwindow, black);  % black out the offscreen window
Screen('FillRect', backwindow, 1, LTrect);
blk_shapes_LT = Screen('GetImage', backwindow, mainsrc);

Screen('FillRect', backwindow, black);  % black out the offscreen window
Screen('FillRect', backwindow, 1, RBrect);
blk_shapes_RB = Screen('GetImage', backwindow, mainsrc);

Screen('FillRect', backwindow, black);  % black out the offscreen window
Screen('FillRect', backwindow, 1, LTskinny);
blk_shapes_LTsk = Screen('GetImage', backwindow, mainsrc);

Screen('FillRect', backwindow, black);  % black out the offscreen window
Screen('FillRect', backwindow, 1, RBskinny);
blk_shapes_RBsk = Screen('GetImage', backwindow, mainsrc);

Screen('FillRect', backwindow, black);  % black out the offscreen window
Screen('FillPoly', backwindow, 1, pointListUL);
Screen('FillPoly', backwindow, 1, pointListUR);
Screen('FillPoly', backwindow, 1, pointListLR);
Screen('FillPoly', backwindow, 1, pointListLL);
Screen('FillRect', backwindow, 1, LTrect);
Screen('FillRect', backwindow, 1, RBrect);
Screen('FillRect', backwindow, 1, LTskinny);
Screen('FillRect', backwindow, 1, RBskinny);
blk_shapes_all_nocc = Screen('GetImage', backwindow, mainsrc);
Screen('FillRect', backwindow, black);  % black out the offscreen window

%Make an all grey, an all white, and an all black (0s) base image of the
%proper dimensions for constructing masks..

%zero_base = zeros(54,54,3);

zero_base = zeros(maindim+4,maindim+4,3); % orig..w/+4s..
%zero_base = zeros(maindim,maindim,3);


zero_base = cast(zero_base,'uint8');

%one_base = ones(54,54,3);

one_base = ones(maindim+4,maindim+4,3); % orig..w/+4s..
%one_base = ones(maindim,maindim,3);

one_base = cast(one_base,'uint8');
grey_base = one_base * grey;
grey_base = cast(grey_base,'uint8');
Screen('FillRect', backwindow, grey);  % grey out the offscreen window

%Make 2 binary masks (inverted versions of each other) for each shape
%individually and for all shapes combined..
% ** For all masks: in the "_1" version, the shapes are 1s and background is 0s
% ** and in the "_2" version, the shapes are 0s and the background is 1s..
% ** masks with a "_gry" tag have a grey background with shapes cut out of
% it as patchs of 0s... (like a grey cookie tray to place the image pieces in...)

mask_all_1 = zero_base + blk_shapes_all;  

%mask_all_2 = ~mask_all_1(1:end,1:54);
mask_all_2 = ~mask_all_1(1:end,1:maindim+4);

mask_all_2 = cast(mask_all_2,'uint8');

mask_ul_1 = zero_base + blk_shapes_ul;

%mask_ul_2 = ~mask_ul_1(1:end,1:54);
mask_ul_2 = ~mask_ul_1(1:end,1:maindim+4);

mask_ul_2 = cast(mask_ul_2,'uint8');

mask_ur_1 = zero_base + blk_shapes_ur;

%mask_ur_2 = ~mask_ur_1(1:end,1:54);
mask_ur_2 = ~mask_ur_1(1:end,1:maindim+4);

mask_ur_2 = cast(mask_ur_2,'uint8');

mask_ll_1 = zero_base + blk_shapes_ll;

%mask_ll_2 = ~mask_ll_1(1:end,1:54);
mask_ll_2 = ~mask_ll_1(1:end,1:maindim+4);

mask_ll_2 = cast(mask_ll_2,'uint8');

mask_lr_1 = zero_base + blk_shapes_lr;

%mask_lr_2 = ~mask_lr_1(1:end,1:54);
mask_lr_2 = ~mask_lr_1(1:end,1:maindim+4);

mask_lr_2 = cast(mask_lr_2,'uint8');

% Do the same for the un-occluded versions..
mask_all_1_nocc = zero_base + blk_shapes_all_nocc;  

%mask_all_2_nocc = ~mask_all_1_nocc(1:end,1:54);
mask_all_2_nocc = ~mask_all_1_nocc(1:end,1:maindim+4);

mask_all_2_nocc = cast(mask_all_2_nocc,'uint8');

mask_LT_1 = zero_base + blk_shapes_LT;

%mask_LT_2 = ~mask_LT_1(1:end,1:54);
mask_LT_2 = ~mask_LT_1(1:end,1:maindim+4);

mask_LT_2 = cast(mask_LT_2,'uint8');

mask_RB_1 = zero_base + blk_shapes_RB;

%mask_RB_2 = ~mask_RB_1(1:end,1:54);
mask_RB_2 = ~mask_RB_1(1:end,1:maindim+4);

mask_RB_2 = cast(mask_RB_2,'uint8');

mask_LTsk_1 = zero_base + blk_shapes_LTsk;

%mask_LTsk_2 = ~mask_LTsk_1(1:end,1:54);
mask_LTsk_2 = ~mask_LTsk_1(1:end,1:maindim+4);

mask_LTsk_2 = cast(mask_LTsk_2,'uint8');

mask_RBsk_1 = zero_base + blk_shapes_RBsk;

%mask_RBsk_2 = ~mask_RBsk_1(1:end,1:54);
mask_RBsk_2 = ~mask_RBsk_1(1:end,1:maindim+4);

mask_RBsk_2 = cast(mask_RBsk_2,'uint8');

%% (If Texture..) Use the masks created above to mask the images into the shapes
%Use the masks above to generate masked versions of the input images. Also
%generate grey background masks for the masked image shapes to fit in (ie zeros where the shapes are)...
%Then put the masked image shapes into the "cookie trays" (grey background
%masks)

%occluded
im_mask_ul_1 = mask_ul_1 .* UL_im;
im_mask_ul_gry = mask_ul_2 .* grey_base;
im_mask_ul_2 = im_mask_ul_1 + im_mask_ul_gry;

im_mask_ur_1 = mask_ur_1 .* UR_im;
im_mask_ur_gry = mask_ur_2 .* grey_base;
im_mask_ur_2 = im_mask_ur_1 + im_mask_ur_gry;

im_mask_ll_1 = mask_ll_1 .* LL_im;
im_mask_ll_gry = mask_ll_2 .* grey_base;
im_mask_ll_2 = im_mask_ll_1 + im_mask_ll_gry;

im_mask_lr_1 = mask_lr_1 .* LR_im;
im_mask_lr_gry = mask_lr_2 .* grey_base;
im_mask_lr_2 = im_mask_lr_1 + im_mask_lr_gry;


im_mask_all_gry = mask_all_2 .* grey_base;
im_masked_all_occ = im_mask_all_gry + im_mask_lr_1 + im_mask_ll_1 + im_mask_ur_1 + im_mask_ul_1; % This is output image if this is an occluded condition

%non-occluded
im_mask_LT_1 = mask_LT_1 .* LT_im;
im_mask_LT_gry = mask_LT_2 .* grey_base;
im_mask_LT_2 = im_mask_LT_1 + im_mask_LT_gry;

im_mask_RB_1 = mask_RB_1 .* RB_im;
im_mask_RB_gry = mask_RB_2 .* grey_base;
im_mask_RB_2 = im_mask_RB_1 + im_mask_RB_gry;

im_mask_LTsk_1 = mask_LTsk_1 .* LTsk_im;
im_mask_LTsk_gry = mask_LTsk_2 .* grey_base;
im_mask_LTsk_2 = im_mask_LTsk_1 + im_mask_LTsk_gry;

im_mask_RBsk_1 = mask_RBsk_1 .* RBsk_im;
im_mask_RBsk_gry = mask_RBsk_2 .* grey_base;
im_mask_RBsk_2 = im_mask_RBsk_1 + im_mask_RBsk_gry;

im_mask_all_nocc_gry = mask_all_2_nocc .* grey_base;
im_masked_all_nocc = im_mask_all_nocc_gry + im_mask_lr_1 + im_mask_ll_1 + im_mask_ur_1 + im_mask_ul_1 + im_mask_LT_1 + im_mask_RB_1 + im_mask_LTsk_1 + im_mask_RBsk_1; % This is output image if this is a non-occluded condition
end

%% Generate the final texture-image-only image if texture is requested..
if T == 1
% grey out the offscreen window
Screen('FillRect', backwindow, grey); 

% draw occluders on offscreen window
% Screen('FillRect', backwindow, gray, vertOcc);
% Screen('FillRect', backwindow, gray, horzOcc);

if occ==0 % occluded condition
    % draw occluded version of the masked images to the offscreen window
    
    if normalize == 0
        image_tex = Screen('MakeTexture', window, im_masked_all_occ);
        Screen('DrawTexture', backwindow, image_tex);
        texOnly_im =Screen('GetImage', backwindow); %capture image of backwindow as it stands now...
    end


    if normalize == 1
        im_masked_all_occ = rescale(im_masked_all_occ,norm_range(1),norm_range(2)); 
        image_tex = Screen('MakeTexture', window, im_masked_all_occ);
        Screen('DrawTexture', backwindow, image_tex); 
        texOnly_im =Screen('GetImage', backwindow); %capture image of backwindow as it stands now...
    end
elseif occ==1 % non-occluded condition
    % draw non-occluded version of the masked images to the offscreen window
    
    if normalize == 0
        image_tex = Screen('MakeTexture', window, im_masked_all_nocc);
        Screen('DrawTexture', backwindow, image_tex)
        texOnly_im =Screen('GetImage', backwindow); %capture image of backwindow as it stands now...
    end

    if normalize == 1
        im_masked_all_nocc = rescale(im_masked_all_nocc,norm_range(1),norm_range(2));
        image_tex = Screen('MakeTexture', window, im_masked_all_nocc);
        Screen('DrawTexture', backwindow, image_tex); 
        texOnly_im =Screen('GetImage', backwindow); %capture image of backwindow as it stands now...
    end
end
end


%% Combine the Luminance, Texture, and Common Region Versions as Requested into Final Output in Offscreen Window..

%First, clear the windows..
Screen('FillRect', backwindow, grey); 
Screen('FillRect', window, grey); 

% Check if both texture AND luminance cues are requested..
% if so, take the mean of the the texture- and the luminance-only
% images to apply the luminance modulations to the texture patches..
if (L == 1) && (T == 1)
    LumTexCombo = 1;    
    texWt = 1-lumWt; % set texture image weight to be 1-lum weight..
    %take weighted average of lum-only and tex-only images...
    LumTexture_im = uint8(texWt*double(texOnly_im) + lumWt*double(lumOnly_im));
    %LumTexture_im = LumTexture_im./(1+lumWt);
    image_tex = Screen('MakeTexture', window, LumTexture_im);
    Screen('DrawTexture', backwindow, image_tex);    
else
    LumTexCombo = 0;
end

% Now, if (and only if) it was NOT the case that BOTH luminance AND texture
% were requested, check if either was chosen individually
if LumTexCombo == 0
    if L == 1
        image_tex = Screen('MakeTexture', window, lumOnly_im);
        Screen('DrawTexture', backwindow, image_tex); 
    end

    if T == 1
        image_tex = Screen('MakeTexture', window, texOnly_im);
        Screen('DrawTexture', backwindow, image_tex); 
    end
    
end

% write CRbox to screen if common region was specified .. 
if CR == 1
Screen('FrameRect', backwindow, CR_penLum, CRbox,penWidth);
end

%% Debugging Sanity-Check Report Out Command Line Messages/Figures
% %For debugging/Dev.. run this to show backwindow as image
% test_im=Screen('GetImage', backwindow); %capture image of backwindow as it stands now...
% imshow(test_im,'InitialMagnification',200)
if cmdLineOut == 1 
disp("----------------------------------")
disp("Condition-General Parameters:")
disp("----------------------------------")
disp(strcat("'orient' is: ", orient))
disp(strcat("'occ' is: ", num2str(occ)))
disp("----------------------------------")
disp("Common-Region-Specific Parameters:")
disp("----------------------------------")
disp(strcat("'CR' is: ", num2str(CR)))
disp(strcat("'CR_obj' is: ", num2str(CR_obj)))
disp(strcat("'CR_al' is: ", num2str(CR_al)))
disp(strcat("'CR_penLum' is: ", num2str(CR_penLum)))
disp(strcat("'CR_penWidth' is: ", num2str(CR_penWidth)))
disp("----------------------------------")
disp("Texture-Specific Parameters:")
disp("----------------------------------")
disp(strcat("'T' is: ", num2str(T)))
disp(strcat("'T_al' is: ", num2str(T_al)))
disp(strcat("'rscale_f1' is: ", num2str(rscale_f1)))
disp(strcat("'rscale_f2' is: ", num2str(rscale_f2)))
disp(strcat("'normalize' is: ", num2str(normalize)))
disp(strcat("'norm_range' is: ", num2str(norm_range)))
disp("----------------------------------")
disp("Luminance-Specific Parameters:")
disp("----------------------------------")
disp(strcat("'L' is: ", num2str(L)))
disp(strcat("'L_al' is: ", num2str(L_al)))
disp(strcat("'luminance 1' is: ", num2str(gray1)))
disp(strcat("'luminance 2' is: ", num2str(gray2)))
disp("----------------------------------")
end



% % draw fixation box on offscreen window
% Screen('DrawLines', backwindow, fixlines, penW, black); % fixation box