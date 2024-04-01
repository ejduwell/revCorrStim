function imOutMat = RevCorrBaseImGen_v2(parsIn,window,change_screensize,scale_f,imOutMat,method,baseOutDir,norm,fxnParz)
%function img_snap = RevCorrBaseImGen_v1(orient,objcon,occluders,change_screensize,scale_f,CR,CR_penWidth,CR_penLum,CR_al,CR_obj,cr_long,cr_short,T,T_al,image_in,imName,rscale_f1,rscale_f2,L,lum1,lum2,lumWt,L_al,window,outdir,outdir_base,basename,out_fmt,cmdLineOut,crop,cropsz)

%% Some standard things that we've got to set up
%--------------------------------------------------------------
%Screen('Preference','Tick0Secs',nan);
%ListenChar(2);  % stop throwing characters to matlab windows     %EJD commented
%HideCursor;                                                      %EJD commented
%GetSecs;						% Pre-load and pre-allocate
%CharAvail;						% Pre-load and pre-allocate
%KbName('UnifyKeyNames');
%KbCheck(-1);

%%% choose appropriate resolution
%critstimdisplay = 100;
%rwidth = 1024;	% requested resolution width
%rheight = 768;	% requested resolution height
% res1=Screen('Resolutions', 0); %list of all resolutions possible
% res2=find([res1(:).width]==rwidth & [res1(:).height]==rheight);
% res3=unique([res1(res2).hz]);
% for r=1:length(res3)
%     rmod(r)=mod(critstimdisplay,1000*(1/res3(r)));  % match frame rate to value stored in critstimdisplay
% end
% [x,rmin]=min(rmod); % best frequency match
% res4=find([res1(:).width]==rwidth & [res1(:).height]==rheight & [res1(:).hz]==res3(rmin)); % best freq. & size
% [x,rmax]=max([res1(res4).pixelSize]);   % largest pixelSize value of good resolutions
% newres = res1(res4(rmax));
%oldResolution=Screen('Resolution', 0, newres.width, newres.height,
%newres.hz, newres.pixelSize); %set it %EJD COMMENTED BC of ERROR..
%HEADS UP THIS MAY BREAK STUFF

%%% get some info. on computer running this script
comp = Screen('Computer');
[ret, cname] = system('hostname');
if comp.windows
    kboard = -1;    % Windows PsychToolbox doesn't support PsychHID functions yet!
    if strcmp(cname(1:12),'behrmann-255')   %far-back room(new)
        room = 'farback';
        directory = 'C:\Documents and Settings\Adam\My Documents\My Dropbox\OBA\PORC4\Data\';
    elseif strcmp(cname(1:12),'behrmann-BL2')   % new back-left
        room = 'backleft';
        directory = 'C:\Documents and Settings\Adam\My Documents\My Dropbox\OBA\PORC4\Data\';
    elseif strcmp(cname(1:12),'behrmann-FL1')   % new front-left
        room = 'frontleft';
        directory = 'C:\Documents and Settings\Adam\My Documents\My Dropbox\OBA\PORC4\Data\';
    elseif strcmp(cname(1:12),'desk-9e4cbd7')        %far-back room - old
        directory = 'C:\My Dropbox\OBA\PORC4\Data';
    elseif strcmp(cname(1:12),'humphreys_35')    %back-left room - old
        directory = 'C:\Documents and Settings\Administrator\My Documents\My Dropbox\OBA\PORC4\Data\';
    end
    %         subjinfo = strcat(directory,'\',subjinfo1);
elseif strcmp(comp.machineName,'AdamsAppleMBP2')	% this is the new (2009) Intel MacBook Pro
    %         directory = '/Users/agreenb/Documents/Files/CMU/PORC1/Data/';comp.machineName
    directory = '/Users/agreenb/Dropbox/OBA/PORC4/Data/';
    kboard = GetKeyboardIndices([],[],73400320);    % internal laptop keyboard
    room = 'office';
elseif strcmp(comp.machineName,'AdamsAppleMB')	% this is the older (2007) Intel MacBook Pro
    directory = '/Users/agreenb/Dropbox/OBA/PORC4/Data/';
    kboard = GetKeyboardIndices([],[],487587840);    % external laptop keyboard (right USB location)
    %       kboard = GetKeyboardIndices([],[],????);    % external laptop keyboard (left USB location)
    room = 'office';
elseif strcmp(comp.machineName,'kern-mbp2-2021')	% Ethan's Macbook from Kern...
    directory = '/Users/eduwell/OneDrive - mcw.edu/Documents/SNAP/projects/Rev_Corr/code/RevCorr/working_dir/data/';
    %kboard = GetKeyboardIndices();    % external laptop keyboard (right USB location)
    %       kboard = GetKeyboardIndices([],[],????);    % external laptop keyboard (left USB location)
    room = 'office';
    Screen('Preference', 'SkipSyncTests', 1); %EJD: added to skip sync tests
    %oldEnableFlag=Screen('Preference', 'EmulateOldPTB'); %EJD added to try to deal with compatability issues with old psychtoolbox...
elseif strcmp(comp.machineName,'BME-TBRC-12341')
     Screen('Preference', 'SkipSyncTests', 1);
     room = 'office';
     kboard = GetKeyboardIndices(); 

elseif strcmp(comp.machineName,'tron')	% dell t7500 tower..
    directory = '/Users/eduwell/OneDrive - mcw.edu/Documents/SNAP/projects/Rev_Corr/code/RevCorr/working_dir/data/';
    %kboard = GetKeyboardIndices();    % external laptop keyboard (right USB location)
    %       kboard = GetKeyboardIndices([],[],????);    % external laptop keyboard (left USB location)
    room = 'office';
    Screen('Preference', 'SkipSyncTests', 1); %EJD: added to skip sync tests
    %oldEnableFlag=Screen('Preference', 'EmulateOldPTB'); %EJD added to try to deal with compatability issues with old psychtoolbox...
else

    fprintf('******************************\nYou''re not on a known machine!\n******************************\n');

    junk = input('However... If you would like to proceed, type Y, if not press N, and the program will terminate.. ','s');

    if junk == 'Y'
        fprintf('******************************\nOKIE DOKIE HERE WE GO....KIE!\n******************************\n');
    else
        return;
    end

end


%Screen('Preference','Tick0Secs',nan);

ListenChar(2);  % stop throwing characters to matlab windows
%GetSecs;						% Pre-load and pre-allocate
%CharAvail;						% Pre-load and pre-allocate
%KbName('UnifyKeyNames');
%KbCheck(-1);

screenrect = Screen('Rect',window);	% get the size of the display screen

%EJD COMMENTED TO SEE IF THIS FIXES BUG WITH STIM WINDOWING/NOT IN CENTER..
% This fixed it: Lesson-->If window is already resized,we don't need to do
% it here..
% if change_screensize == 1
% 
%     if strcmp(comp.machineName,'tron')
%         scale_f = 0.5;
%         % ejd attempt to rescale whole double screen into a single screen..
%         %screenrect = screenrect * scale_f; %ORIG
%         screenrect(3) = screenrect(3) * scale_f;
%     else
%         scale_f = 0.5;
%         screenrect = screenrect * scale_f;
%     end
% 
% end

scrn_wid = screenrect(3);
scrn_hei = screenrect(4);
scrnWidHei = [scrn_wid,scrn_hei];

%%% choose appropriate resolution
critstimdisplay = 100;
rwidth = screenrect(3);	% requested resolution width
rheight = screenrect(4);	% requested resolution height
zoom_factor = screenrect(4)/rheight;
% res1=Screen('Resolutions', 0); %list of all resolutions possible
% res2=find([res1(:).width]==rwidth & [res1(:).height]==rheight);
% res3=unique([res1(res2).hz]);
% 
% for r=1:length(res3)
%     rmod(r)=mod(critstimdisplay,1000*(1/res3(r)));  % match frame rate to value stored in critstimdisplay
% end
% [x,rmin]=min(rmod); % best frequency match
% res4=find([res1(:).width]==rwidth & [res1(:).height]==rheight & [res1(:).hz]==res3(rmin)); % best freq. & size
% [x,rmax]=max([res1(res4).pixelSize]);   % largest pixelSize value of good resolutions
% newres = res1(res4(rmax));

PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'UseFastOffscreenWindows');

% window = PsychImaging('OpenWindow', 0, GrayIndex(0,0.5), screenrect);	% Open generic on-screen window
% Screen(window,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);    % allow for transparency!

% main parameters (could prompt for/calculate)
% Make the height of the screen be the hypotenuse of the right angles
% forming the main dimension of the stim window such that when rotated 45 deg,
% it will fit the screen exactly..
maindim = round(sqrt((rheight^2)/2));	% dimension of mainrect size;max=750 (must jive with surfpatch=>5*surfpatch=maindim)
maindim = round(maindim*1.4); % now make it slightly bigger..

%surfpatch = round(maindim);	% dimension of stimulus surface patch squares
surfpatch = round(maindim/5);	% dimension of stimulus surface patch squares

occlength = round(surfpatch*0.7); % length of occluder outlines
magzoom = zoom_factor;   % zooming/magnification multiplier for display of stimulus

%penW = 2;
penW=fxnParz(1);
penH = 2;

trial_con = ['z','m'];  %!! : ETHAN GENERALIZED..  this had been "orientation"
% is now "trial_con" (for trial condition).. characters specify both the
% button press and which type of trial all in one.. 
% Now (for face version): trial_con(1)=APE ONLY, trial_con(2)=HUMAN PRESENT.. 

key1 = KbName(trial_con(1)); % NOW: APE ONLY WAS:Vertical (LL-UR) orientation
key2 = KbName(trial_con(2)); % NOW: APE+HUMAN WAS:Horizontal (UL-LR) orientation

%fixdim = round(maindim/50);	% dimension of fixation box
fixdim = round(maindim/fxnParz(2));	% dimension of fixation box


% determine rect sizes for stimuli
mainrect = CenterRect([0 -maindim maindim 0],screenrect); % rect defining outer corners of occluded objects
mainsrc = CenterRect([0 -(maindim+4) maindim+4 0],screenrect); % rect defining outer corners of non-zoomed display
maindest = CenterRect([0 -magzoom*(maindim+4) magzoom*(maindim+4) 0],screenrect); % rect defining outer corners of zoomed display

[x,y]=RectCenter(screenrect);
bigbox = CenterRectOnPoint([0 0 3*surfpatch 3*surfpatch],x,y); % box that just includes the four surface patches
boxUL = [bigbox(RectLeft) bigbox(RectTop) bigbox(RectLeft)+surfpatch bigbox(RectTop)+surfpatch];
boxUR = [bigbox(RectRight)-surfpatch bigbox(RectTop) bigbox(RectRight) bigbox(RectTop)+surfpatch];
boxLL = [bigbox(RectLeft) bigbox(RectBottom)-surfpatch bigbox(RectLeft)+surfpatch bigbox(RectBottom)];
boxLR = [bigbox(RectRight)-surfpatch bigbox(RectBottom)-surfpatch bigbox(RectRight) bigbox(RectBottom)];
vertOcc = [boxUL(RectRight) mainrect(RectTop) boxUR(RectLeft) mainrect(RectBottom)]; % rect for vertical occluder
horzOcc = [mainrect(RectLeft) boxUL(RectBottom) mainrect(RectRight) boxLL(RectTop)]; % rect for horizontal occluder
fixrect = floor(CenterRectOnPoint([0 -fixdim fixdim 0],x,y));
boxrect = [0 0 RectWidth(boxUL) RectHeight(boxUL)];	% size of any of the surface patches (for use when drawing cue around any patch)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% create fixation box for putting on the screen
minfixrect = [0 0 RectWidth(fixrect) RectHeight(fixrect)];
fixlines(:,1:2) = [fixrect(1)-round(penW/2) fixrect(3)+round(penW/2);fixrect(2) fixrect(2)];   % top of fixation square
fixlines(:,3:4) = [fixrect(1)-round(penW/2) fixrect(3)+round(penW/2);fixrect(4) fixrect(4)];   % bottom of fixation square
fixlines(:,5:6) = [fixrect(1) fixrect(1);fixrect(2) fixrect(4)];   % left of fixation square
fixlines(:,7:8) = [fixrect(3) fixrect(3);fixrect(2) fixrect(4)];   % right of fixation square
fixlines(:,9:10) = [fixrect(1) fixrect(3);fixrect(2) fixrect(4)];   % upper left to lower right of fixation square
fixlines(:,11:12) = [fixrect(1) fixrect(3);fixrect(4) fixrect(2)];   % lower left to upper right of fixation square
%%% create occluder outlines
occlines(:,1:2) = [vertOcc(1) vertOcc(1); vertOcc(2)-1+occlength vertOcc(2)-1];
occlines(:,3:4) = [vertOcc(1) vertOcc(3)+1; vertOcc(2)-1 vertOcc(2)-1];
occlines(:,5:6) = [vertOcc(3)+1 vertOcc(3)+1; vertOcc(2)-1+occlength vertOcc(2)-1];
occlines(:,7:8) = [vertOcc(1) vertOcc(1); vertOcc(4)+1-occlength vertOcc(4)+1];
occlines(:,9:10) = [vertOcc(1) vertOcc(3)+1; vertOcc(4) vertOcc(4)];
occlines(:,11:12) = [vertOcc(3)+1 vertOcc(3)+1; vertOcc(4)+1 vertOcc(4)+1-occlength];

occlines(:,13:14) = [horzOcc(1)-1+occlength horzOcc(1)-1; horzOcc(4) horzOcc(4)];
% occlines(:,15:16) = [horzOcc(1)-1 horzOcc(1)-1; horzOcc(4)+1 horzOcc(2)-1]; % this line for WinTels
% occlines(:,15:16) = [horzOcc(1) horzOcc(1); horzOcc(4)+1 horzOcc(2)-1]; % this line for MacBookPro
occlines(:,17:18) = [horzOcc(1)-1 horzOcc(1)-1+occlength; horzOcc(2)-1 horzOcc(2)-1];
occlines(:,19:20) = [horzOcc(3)+1-occlength horzOcc(3)+1; horzOcc(4) horzOcc(4)];
% occlines(:,21:22) = [horzOcc(3) horzOcc(3); horzOcc(4)+1 horzOcc(2)-1]; % this line for WinTels
% occlines(:,21:22) = [horzOcc(3)+1 horzOcc(3)+1; horzOcc(4)+1 horzOcc(2)-1];   % this line for MacBookPro
occlines(:,23:24) = [horzOcc(3)+1 horzOcc(3)+1-occlength; horzOcc(2)-1 horzOcc(2)-1];
switch room
    case 'frontleft'
        occlines(:,15:16) = [horzOcc(1) horzOcc(1); horzOcc(4)+1 horzOcc(2)-1]; % this line for MacBookPro
        occlines(:,21:22) = [horzOcc(3)+1 horzOcc(3)+1; horzOcc(4)+1 horzOcc(2)-1];   % this line for MacBookPro
    case 'backleft'
        occlines(:,15:16) = [horzOcc(1)-1 horzOcc(1)-1; horzOcc(4)+1 horzOcc(2)-1]; % this line for WinTels
        occlines(:,21:22) = [horzOcc(3) horzOcc(3); horzOcc(4)+1 horzOcc(2)-1]; % this line for WinTels
    case 'farback'
        occlines(:,15:16) = [horzOcc(1)-1 horzOcc(1)-1; horzOcc(4)+1 horzOcc(2)-1]; % this line for WinTels
        occlines(:,21:22) = [horzOcc(3) horzOcc(3); horzOcc(4)+1 horzOcc(2)-1]; % this line for WinTels
    case 'office'
        occlines(:,15:16) = [horzOcc(1) horzOcc(1); horzOcc(4)+1 horzOcc(2)-1]; % this line for MacBookPro
        occlines(:,21:22) = [horzOcc(3)+1 horzOcc(3)+1; horzOcc(4)+1 horzOcc(2)-1];   % this line for MacBookPro
end


white = WhiteIndex(window);
black = BlackIndex(window);
grayblack = floor(GrayIndex(window,0.15));
chargray = ceil(GrayIndex(window,0.30));
newgray = ceil(GrayIndex(window,0.35));
lgray = ceil(GrayIndex(window,0.75)); % light gray
dgray = floor(GrayIndex(window,0.5)); % dark gray
gray = dgray;	% medium gray
colors = {'blue','purple','orange','yellow','green','army','rose','cyan','red','gray','white'};
red = [255, 0, 0, 255];	% red
green = [0, 255, 0, 255];	% green
blue = [0, 0, 255, 255];	% blue
purple = [159, 31, 220, 255];	% purple
%purple = [200, 0, 200, 255];	% purple
orange = [255, 140, 0, 255];	% orange
yellow = [255, 255, 0, 255];	% yellow
army = [153, 204, 50, 255];	% army
rose = [255, 192, 192, 255];% rose
cyan = [0, 255, 255, 255];	% cyan


gray1_occ = chargray;   % create luminance values
gray2_occ = newgray;   % create luminance values
gray1_nocc = chargray;   % create luminance values
gray2_nocc = newgray;   % create luminance values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a little bit of texture generation for the displays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open an offscreen window for background
backwindow_fix=Screen('OpenOffscreenWindow', window, gray);
backwindow=Screen('OpenOffscreenWindow', window, gray);
backwindow_cap=Screen('OpenOffscreenWindow', window, gray);


%% Begin Main Stimulus Image Generation Loop
for ii =1:size(parsIn,1)
    if parsIn{ii,1}~="NA"
try
%% Parameters
%% Update pars below based on parsIn row for this pass..
orient=parsIn{ii,1};
objcon=parsIn{ii,2};
occluders=parsIn{ii,3};
change_screensize=parsIn{ii,4};
scale_f=parsIn{ii,5};
CR=parsIn{ii,6};
CR_penWidth=parsIn{ii,7};
CR_penLum=parsIn{ii,8};
CR_al=parsIn{ii,9};
CR_obj=parsIn{ii,10};
cr_long=parsIn{ii,11};
cr_short=parsIn{ii,12};
T=parsIn{ii,13};
T_al=parsIn{ii,14};
image_in=parsIn{ii,15};
imName=parsIn{ii,16};
rscale_f1=parsIn{ii,17};
rscale_f2=parsIn{ii,18};
L=parsIn{ii,19};
lum1=parsIn{ii,20};
lum2=parsIn{ii,21};
lumWt=parsIn{ii,22};
L_al=parsIn{ii,23};
window=parsIn{ii,24};
outdir=parsIn{ii,25};
outdir_base=parsIn{ii,26};
basename=parsIn{ii,27};
out_fmt=parsIn{ii,28};
cmdLineOut=parsIn{ii,29};
crop=parsIn{ii,30};
cropsz=parsIn{ii,31};
%% Output File Parameters
%basename = "revcorrBI";
%out_fmt = "png";

%% System Parameters
%--------------------------------------------------------------
%outdir = "out_test"; % Specify output directory (which will be created)
%outdir_base = '/Users/eduwell/OneDrive - mcw.edu/Documents/SNAP/projects/Rev_Corr/code/RevCorr/working_dir'; % specify base directory in which outdir will be created/saved..
start_dir = pwd; % save the initial starting directory..

% get unique timestamp
%today = datestr(now,30);	% get current time/date stamp for filename
db_mode = 0;
%% Condition General Parameters
%--------------------------------------------------------------

% orient = 'm'; % 'm' or 'z'
% objcon = 1; % 1 or 2

% Load above selections into trialprop
%% -------------------------------------
trialprop = cell(zeros(1,2)); % initialize
trialprop{1} = orient;
trialprop{2} = objcon;


%% -------------------------------------

%occluders = 0; % include occluders? 0 = yes, 1 = no

%change_screensize = 0;
%scale_f = 0.5;

%% Common Region Parameters
%--------------------------------------------------------------
% CR = 0; %Include Common Region Box?
% CR_penWidth = 1; % If so, specify pen width for common region box..
% CR_penLum = 75; % specify the luminance value for the common region pen strokes too..
% CR_al = 1; % Specify whether the common region should align with the objects (1) or be perpendicular (0)
% CR_obj = 2; % Specify which "object" the common region should surround (2 for upper object or 1 for lower object)
% cr_long = 32;
% cr_short = 14;

%% Texture Parameters
%--------------------------------------------------------------
% T = 1; % Run texture version? (1=yes, 0=no)
% T_al = 1; % Specify whether the texture should align with the objects (1) or be perpendicular (0)
% image_in = '/Users/eduwell/OneDrive - mcw.edu/Documents/SNAP/projects/Rev_Corr/code/RevCorr/working_dir/spheres_pic_rs2_nearest.jpg';
% imName = "spheres";
% rscale_f1 = 1.5; % texture scaling factor 1
% rscale_f2 = 0.5; % texture scaling factor 2

%generate strings for outfile 
if T==1
scalefactor1_str = num2str(rscale_f1);
%scalefactor1_str = strrep(scalefactor1_str,'.','_');
scalefactor2_str = num2str(rscale_f2);
%scalefactor2_str = strrep(scalefactor2_str,'.','_');
else
    scalefactor1_str="NA";
    scalefactor2_str="NA";
end
% Texture image normalization params:
normalize = norm;
norm_range = [0,255];

%% Luminance Parameters
%--------------------------------------------------------------
% L = 1; % Run luminance version? (1=yes, 0=no)
% lum1 = 120; % num between 0-255
% lum2 = 120; % num between 0-255
% lumWt = 0.7; % luminance weight (for texture/luminance combo.. note: texture weight will be 1-lumWt..) (num between 0-1)

%generate strings for outfile w/out decimal..
lumWt_str = num2str(lumWt);
lumWt_str = strrep(lumWt_str,'.','_');

% L_al = 0; % Specify whether the luminance should align with the objects (1) or be perpendicular (0)

%% Code
cd(outdir_base)

%Create output directory if it does not already exist
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

% enter output dir
cd(outdir)

% Clear the windows..
Screen('FillRect', window, gray);  % blank out the onscreen window
%Screen('FillRect', backwindow, gray);  % blank out the offscreen window
Screen('FillRect', backwindow_cap, gray);  % blank out the offscreen window

% Create Base Image
[backwindow] = PORClinkLTCR_v1(mainsrc, window, backwindow,vertOcc,horzOcc,boxUL,boxUR,boxLL,boxLR,gray,lum1,lum2,trialprop,occluders,CR,CR_penWidth,CR_penLum,CR_al,CR_obj,scrnWidHei,rscale_f1, rscale_f2, image_in, db_mode,normalize,norm_range,maindim,T,T_al,L,L_al,lumWt,cr_long,cr_short,cmdLineOut);

% Add the fixation point and lines
[backwindow_fix,backwindow] = PORCfix4a(backwindow_fix,backwindow,fixrect,occlines,penW,grayblack);

% Write to main window/rotate..
Screen('DrawTexture', backwindow_cap, backwindow, mainsrc, maindest, 45, 0, 1);  % draw stimulus to the on-screen window, rotated & zoomed
% Screen('Flip', window);% show linking display %

% Capture window as image..
img_snap =Screen('GetImage', backwindow_cap); %capture image of backwindow as it stands now...
% Convert to greyscale
img_snap = rgb2gray(img_snap);


% Save image to file
% contatenate parts extracted above along with the interp method 
% and format extension to form the descriptive output filename..
out_name = strcat(basename,"_occ_",num2str(occluders),"_ori_",orient,"_","_CR_",num2str(CR),"_CRpw_",num2str(CR_penWidth),"_CRpl_",num2str(CR_penLum),"_CRal_",num2str(CR_al),"_CRob_",num2str(CR_obj),"_CRdm_",num2str(cr_long),num2str(cr_short),"_T_",num2str(T),imName,"_Tr1_",scalefactor1_str,"_Tr2_",scalefactor2_str,"_Tal_",num2str(T_al),"_L_",num2str(L),"_lum1_",num2str(lum1),"_lum2_",num2str(lum2),"_LWT_",num2str(lumWt),"_Lal_",lumWt_str,"_Ocon_",num2str(objcon),".",out_fmt);

if crop == 1
    img_snap = imgSizeEdit(img_snap,cropsz,method);
end

addNoize=0;
ranPerPass=1;
specificNoiseNum=42;

if addNoize==1
% random_integer = randi([1, 1000]);
% noiseimage=strcat("/Users/eduwell/Library/CloudStorage/SynologyDrive-SNAP/projects/revCorrStim/noise/white_test1/noiseSample_",num2str(random_integer),".jpg");
% noiseimage=imread(noiseimage);

if ranPerPass==1
random_integer = randi([1, 100]);
%random_integer = randi([1, 1000]);
else
random_integer=specificNoiseNum;
end

% noiseimage=strcat("/Users/eduwell/Library/CloudStorage/SynologyDrive-SNAP/projects/revCorrStim/noise/whiteNoise512-512/noiseSample_",num2str(random_integer),".jpg");
% noiseimage=imread(noiseimage);

noiseimage=strcat("/Users/eduwell/matlabSandBox/texture/textureSynth/BaseImg512by512/noise/wtdAvgSmplz1-100_r",num2str(random_integer),".png");
noiseimage=imread(noiseimage);

[img_snap] = revCorrAddNoise_v1(img_snap,noiseimage,0.5);
img_snap = im2uint8(img_snap);

end

%Add this pass to imOutMat
imOutMat{1,ii} = img_snap;
imOutMat{2,ii} = {parsIn{ii,:}};

imwrite(img_snap,out_name,out_fmt);

% Clear the windows..
%Screen('FillRect', window, gray);  % blank out the onscreen window
Screen('FillRect', backwindow, gray);  % blank out the offscreen window
Screen('FillRect', backwindow_cap, gray);  % blank out the offscreen window

catch
rethrow(lasterror);
Screen('CloseAll');
% Clear up ..
cd(start_dir);
ListenChar(0);
ShowCursor;
%clear
end
    else
        imOutMat{1,ii} = "NA";
        imOutMat{2,ii} = "NA";
    end
end
%GENERATE THE FIXATION ONLY IMAGE
RevCorFixImGen(window, backwindow_fix,backwindow,backwindow_cap, fixrect,occlines,penW,grayblack,crop,cropsz,method,gray,mainsrc, maindest,baseOutDir,out_fmt)

% Clear up ..
cd(start_dir);
ListenChar(0);
ShowCursor;
Screen('CloseAll');
Screen('Close');
%clear
end