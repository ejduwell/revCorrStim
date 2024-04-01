function GenRevCorrBaseIms
% Wrapper function for calling RevCorrBaseImGen_v1 to generate stand-alone reverse correlation base images. 
% Specify parameters below in parameters section. These are fed to RevCorrBaseImGen_v1.m
% Other dependencies: PORClinkLTCR_v1.m, PORCfix4a.m
% RevCorrBaseImGen_v1 to generate a base image file.
%% Set up path
% Detect and set directory/path params:
% -------------------------------------
%Save starting directory location ..
strtDir = pwd;

% get matlab directory path by "which-ing" for this file..
% NOTE: key assumption: there is only one copy of this on the path.. (that
% should always be the case to avoid other conflic;ts..)
thisFile = "GenRevCorrBaseIms.m"; % string should correspond to the name of this file's name..
mlab_dir = fileparts(which(thisFile));

% go to the matlab dir..
cd(mlab_dir)
% -------------------------------------
%% Parameters
    %% Output File Parameters
    basename = "revcorrBI";
    out_fmt = "png";
    
    %% System Parameters
    %--------------------------------------------------------------
    outdir = "out_test"; % Specify output directory (which will be created)
    outdir_base = mlab_dir; % specify base directory in which outdir will be created/saved.. (currently set to the base directory auto-detected above under "Set up path" to avoid path annoyances across machines accessing shared folder..)
    
    %% Condition-General Parameters
    %--------------------------------------------------------------
    orient = 'm'; % 'm' or 'z'
    objcon = 2; % 1 or 2
    occluders = 1; % include occluders? 0 = yes, 1 = no
    
    change_screensize = 0; % if 1, will resize the screen dimensions by scale factor "scale_f" defined below.., if 0, does nothing
    scale_f = 0.5; % Screen rescale factor (number between 0 and 1) 
    cmdLineOut =0; % if 1, will print the parameters used for creating the base image out on the command line
    %% Common Region Parameters
    %--------------------------------------------------------------
    CR = 0; %Include Common Region Box?
    CR_penWidth = 1; % If so, specify pen width for common region box..
    CR_penLum = 0; % specify the luminance value for the common region pen strokes too..
    CR_al = 1; % Specify whether the common region should align with the objects (1) or be perpendicular (0)
    CR_obj = 2; % Specify which "object" the common region should surround (2 for upper object or 1 for lower object)
    cr_long = 32; % Long dimension of common region box expressed in % of screen height
    cr_short = 14; % Short dimension of common region box expressed in % of screen height
    
    %% Texture Parameters
    %--------------------------------------------------------------
    T = 0; % Run texture version? (1=yes, 0=no)
    T_al = 1; % Specify whether the texture should align with the objects (1) or be perpendicular (0)
    
    image_in = 'gravel_highres_rs10_nearest.jpg';

    imdirBase = strcat(mlab_dir,"/","texImgs");

    image_in = strcat(imdirBase,"/",image_in); % add the image directory base path to the filename
    
    imName = "spheres";
    rscale_f1 = 1.5; % texture scaling factor 1
    rscale_f2 = 0.5; % texture scaling factor 2
    
    %% Luminance Parameters
    %--------------------------------------------------------------
    L = 1; % Run luminance version? (1=yes, 0=no)
    lum1 = 90; % num between 0-255 90
    lum2 = 150; % num between 0-255 150
    %lumWt = 0.7; % luminance weight (for texture/luminance combo.. note: texture weight will be 1-lumWt..) (num between 0-1)
    lumWt = 1; % luminance weight (for texture/luminance combo.. note: texture weight will be 1-lumWt..) (num between 0-1)
    L_al = 1; % Specify whether the luminance should align with the objects (1) or be perpendicular (0)

%% Set Up Windows

if change_screensize == 1
    screenrect = screenrect * scale_f;
else
    screenrect = Screen('Rect',0);	% get the size of the display screen
end
scrn_wid = screenrect(3);
scrn_hei = screenrect(4);
scrnWidHei = [scrn_wid,scrn_hei];

%%% choose appropriate resolution
% critstimdisplay = 100;
% rwidth = screenrect(3);	% requested resolution width
% rheight = screenrect(4);	% requested resolution height
% zoom_factor = screenrect(4)/rheight;
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

Screen('Preference', 'SkipSyncTests', 1); %EJD: added to skip sync tests
%window = PsychImaging('OpenWindow', 0, GrayIndex(0,0.5), screenrect);	% Open generic on-screen window
window = PsychImaging('OpenWindow', 0, GrayIndex(0,0.5),screenrect);

% PsychImaging('PrepareConfiguration');
% PsychImaging('AddTask', 'General', 'UseFastOffscreenWindows');
Screen(window,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);    % allow for transparency!

%% Run the RevCorrBaseImGen_v1 with the selected parameters..
RevCorrBaseImGen_v1(orient,objcon,occluders,change_screensize,scale_f,CR,CR_penWidth,CR_penLum,CR_al,CR_obj,cr_long,cr_short,T,T_al,image_in,imName,rscale_f1,rscale_f2,L,lum1,lum2,lumWt,L_al,window,outdir,outdir_base,basename,out_fmt,cmdLineOut)

cd(strtDir); % make sure we end at the starting dir
end