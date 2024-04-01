function [trialtm,ts_out,qs_out,sd_out,nr_ideals, params_out, tdf_out] = RevCorrQuest5_int_v5(room,kboard,numtrials,link,noise,window,image_dir,db_mode,db_mode_screen,nr_params, imtag, pthrds, n_int, qst_startParams,t_prob,maxmin,maxmin_seg,main_dir,olCon,parTagz,obj1Par,obj2Par,qstParStr,noise_dir,nWt,nTag,sessionDir,giveInstrxn)

Screen('Preference', 'SkipSyncTests', 1); %EJD: added to skip sync tests
Screen('Preference','SuppressAllWarnings', 1);
rng("shuffle"); % shuffle random number generator..


% ***NEED TO WORK ON TIMING
%
% **EJD THOUGHT: THE TASK IN THIS CASE IN SOME WAYS INTERUPTS THE QUEST..
% WE'RE TRYING TO DIAL IN THEIR APE-HUMAN FRAME DETECTION THRESHOLD.. BUT
% THE DESCRIMINATION THEYRE MAKING IN THE TASK IS "APE ONLY VS.
% APE-HUMAN".. THIS MEANS THAT THERE WILL, BY DEFINITION BE A PROPONDERANCE
% OF APE-ONLY TRIALS.. AND THAT ON THOSE TRIALS, WE WILL NEED TO DISREGARD
% THE RECOMENDATION OF QUEST.. (AND SIMPLY DISPLAY THE APE ONLY IMAGE..)**
% WHAT DOES THIS MEAN FOR THE FUNCTION/VALIDITY OF THE QUEST?
%
%function[trialtm,t1,t2,t3,t4,t5,q1,q2,q3,q4,q5,sd1,sd2,sd3,sd4,sd5,nr_ideal_t1,nr_ideal_t2,nr_ideal_t3,nr_ideal_t4,nr_ideal_t5,param1,
%param2, param3, param4, param5] =
%FACEquest5_int(room,kboard,numtrials,link,noise,occ,window,image_dir,db_mode,nr_params)%
%ORIG..

%
% PORCquest4a.m
%
% script for Perceptual Organization Reverse Correlation experiments in PTB-3
% this script uses Quest to find the thresholds for linking cue discriminatipracticeon
% Makes calls to: POfix4a.m POlinkB4a.m
%
% this version works well w/ recording responses
%
%
% Written by Adam Greenberg, JHU/PBS
% Oct, 2010 modified from POquest1a
% asg 11/2/2010 modified from PORCquest1a to change pThreshold to 70% (from 75%)
% asg 4/20/2011 modified from PORCquest1b to include Boundary collinearity stimulus
% asg 6/8/2011 modified from PORCquest2a to combine Boundary collinearity & Luminance

%% Debug Mode Params
db_respV = 3; % Specify what sort of automatic response version you want in debug mode: 0 = all random, 1 = all correct, 2 = all incorrect, 3 = incorrect below db_respV3_thr & correct above db_respV3_thr
db_wait = 0; % Specify the amount of time (in sec) you want to wait between automatic responses in db_mode to view frames..
db_ans = 1; % if 1, program will print the answer on each frame during debugging..
%db_respV3_thr = 0.7; % threshold for correct/incorrect responses in db_respV=3
%db_respV3_thr_per = 70; % percent correct at threshold db_respV3_thr
db_press2advance = 0;
db_mode_correct_trials = zeros(1, numtrials);
dbmode_prtThrsh = 0;
dbskip = 1;

% if db_respV = 3, set up response probability function params
if db_respV == 3
    %Orig before Structure was implimented.. (cmntd)
    %     nr_rmax = 0.5;
    %     nr_c50 = 0.5;
    %     nr_b = 0.5;
    %     nr_n = 2;
    %
    % Grab the parameters out of the nr_params structure passed in as input..
    nr_rmax = nr_params.nr_rmax;
    nr_c50 = nr_params.nr_c50;
    nr_b = nr_params.nr_b;
    nr_n = nr_params.nr_n;
end

% Trial Params
presTime=0.6; % time after trial start that stimuli are presented ..(sec)

presTime2=presTime+0.7; % time after trial start that stimuli are removed ..(sec)

resp_cut = 2; % Max time after trial start that responses are recorded

ITI	=0.5; % inter-trial interval time (sec).. (time btw trials)
fback_On = 1; % switch for feedback
fb_correct = '$';
fb_incorrect = ' ';
%fback_txtSize = 40;
feedback = [fb_correct,fb_incorrect]; % feedback characters

if fback_On == 1
    feedbackTime =0.25; % amount of time (sec) feed back is present/displayed
else
    feedbackTime=0;
end

beep_wrong = 0;
try
    %% some standard things that we've got to set up

    Screen('Preference','VisualDebugLevel', 0);
    Screen('Preference','SuppressAllWarnings', 1);
    Screen('Preference','Tick0Secs',nan);
    if db_mode == 0
        HideCursor;
    end
    ListenChar(2);  % stop throwing characters to matlab windows
    GetSecs;						% Pre-load and pre-allocate
    CharAvail;						% Pre-load and pre-allocate
    %     KbName('UnifyKeyNames');
    %     KbCheck(-1);

   %%% choose appropriate resolution
    critstimdisplay = 100;
    rwidth = 1024;	% requested resolution width
    rheight = 768;	% requested resolution height
    
    screenrect = Screen('Rect',0);	% get the size of the display screen
    zoom_factor = screenrect(4)/rheight;
    comp = Screen('Computer');

    if db_mode_screen == 1

        if strcmp(comp.machineName,'smog')
        screenrect(3)=1680;
        end
        %scale_f = 1;

        if strcmp(comp.machineName,'tron')
        scale_f = 0.5;
        % ejd attempt to rescale whole double screen into a single screen..
        %screenrect = screenrect * scale_f; %ORIG
        screenrect(3) = screenrect(3) * scale_f;
        end
    end

    % save some screen values for easy access..
    scrn_top = screenrect(2);
    scrn_bot = screenrect(4);
    scrn_left = screenrect(1);
    scrn_right = screenrect(3);
    scrn_1percY = (abs(scrn_bot - scrn_top))/100;
    scrn_1percX = (abs(scrn_left - scrn_right))/100;

%     PsychImaging('PrepareConfiguration');
%     PsychImaging('AddTask', 'General', 'UseFastOffscreenWindows');
%     window = PsychImaging('OpenWindow', 0, GrayIndex(0,0.5), screenrect);	% Open generic on-screen window
%     Screen(window,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);    % allow for transparency!
%     ifi = Screen('GetFlipInterval', window, 60);
    
    % Start with the "stimulus loading page"..
    white = WhiteIndex(window);
    black = BlackIndex(window);
    gray=white/2;
    % A) Create Page
    Screen('FillRect', window, gray); % gray out window
    % add formatted text..
    txtDistFrmTop1=50*scrn_1percY; % specify distance of text from top as percent of screen dimension..
    Screen('TextSize',window, 50);
    speech1 = ['Task is loading. Please wait... \n'];
    DrawFormattedText(window, speech1,'center', (scrn_top+txtDistFrmTop1), black);
    Screen('Flip', window);% show linking display

    %% main parameters (could prompt for/calculate)
    intro_txtHeaderSize = round(60*(zoom_factor/2));
    intro_txtSize = round(40*(zoom_factor/2));
    fback_txtSize = round(40*(zoom_factor/2));
    %screenrect = Screen('Rect',0);	% get the size of the display screen
    maindim = 768;	% dimension of mainrect size;max=750 (must jive with surfpatch=>5*surfpatch=maindim)

    surfpatch = round(maindim);	% dimension of stimulus surface patch squares
    %surfpatch = round(maindim/5);	% dimension of stimulus surface patch squares

    %fixdim = round(maindim/1000);	% dimension of fixation box
    fixdim = round(maindim/25);	% dimension of fixation box

    occlength = round(surfpatch*0.7); % length of occluder outlines
    %magzoom = zoom_factor;   % zooming/magnification multiplier for display of stimulus % ***!!! EJD CHANGED FROM 1 to 2
    magzoom=1; %EJD TEST
    penW = 1;
    penH = 1;

    trial_con = ['z','m'];  %!! : ETHAN GENERALIZED..  this had been "orientation"
    % is now "trial_con" (for trial condition).. characters specify both the
    % button press and which type of trial all in one..
    % Now (for face version): trial_con(1)=APE ONLY, trial_con(2)=HUMAN PRESENT..

    key1 = KbName(trial_con(1)); % NOW: APE ONLY WAS:Vertical (LL-UR) orientation
    key2 = KbName(trial_con(2)); % NOW: APE+HUMAN WAS:Horizontal (UL-LR) orientation

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check screen resolution
    %[wid,hei]=Screen('WindowSize',0);

    % if db_mode == 0
    %     %HideCursor; %EJD COMMENTED ..
    % end

    % determine rect sizes for stimuli
    mainrect = CenterRect([0 -rheight rwidth 0],screenrect); % rect defining outer corners of occluded objects
    %mainrect = CenterRect([0 -maindim maindim 0],screenrect); % rect defining outer corners of occluded objects

    mainsrc = CenterRect([0 -rheight rwidth 0],screenrect); % rect defining outer corners of non-zoomed display
    %mainsrc = CenterRect([0 -(maindim+4) maindim+4 0],screenrect); % rect defining outer corners of non-zoomed display
    
    %****!!!
    %maindest = CenterRect([0 -rheight rwidth 0],screenrect); % rect defining outer corners of zoomed display % EJD commented.. uncommented line below.. 
    %maindest = CenterRect([0 -magzoom*(maindim+4) magzoom*(maindim+4) 0],screenrect); % rect defining outer corners of zoomed display
    maindest = CenterRect([0 -magzoom*(rheight) magzoom*(rwidth) 0],screenrect); % rect defining outer corners of zoomed display
    %****!!!

    [x,y]=RectCenter(screenrect);
    %bigbox = CenterRectOnPoint([0 0 3*surfpatch 3*surfpatch],x,y); % box that just includes the four surface patches
    bigbox = CenterRectOnPoint([0 0 rheight rwidth],x,y); % box that just includes the four surface patches


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
    % fixlines(:,1:2) = [fixrect(1)-round(penW/2) fixrect(3)+round(penW/2);fixrect(2) fixrect(2)];   % top of fixation square
    % fixlines(:,3:4) = [fixrect(1)-round(penW/2) fixrect(3)+round(penW/2);fixrect(4) fixrect(4)];   % bottom of fixation square
    % fixlines(:,5:6) = [fixrect(1) fixrect(1);fixrect(2) fixrect(4)];   % left of fixation square
    % fixlines(:,7:8) = [fixrect(3) fixrect(3);fixrect(2) fixrect(4)];   % right of fixation square
    % fixlines(:,9:10) = [fixrect(1) fixrect(3);fixrect(2) fixrect(4)];   % upper left to lower right of fixation square
    % fixlines(:,11:12) = [fixrect(1) fixrect(3);fixrect(4) fixrect(2)];   % lower left to upper right of fixation square
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

    % Get the fixation point only image...
    %hard-code-path
    %fix_im_path = "/Users/testadmin/Documents/tDCS_stim_MATLAB/images/square.png"; %Image Path..
    fix_im_path = strcat(image_dir,"/","fixation.png");
    %fix_im_path = strcat(main_dir,"/images/square.png"); %Image Path.. %generalized
    fix_im = imread(fix_im_path);    % read in image

    % Using the fixation only image as a starting point, 
    % generate an image with 0s everywhere but where the fixation point and 
    % inducer lines are (ie everything with a value darker than 127.. 
    % everything that isn't gray..) This will be combined logically with 
    % each BI+Noise image so the appearance of the fixation point/inducers
    % stays constant regardless of the noise..
    fixIndOlayIm=fix_im<127;
    fixIndOlayIm=uint8(double(fixIndOlayIm).*double(fix_im));
    nonFxnIndMask=uint8(fixIndOlayIm==0); %mask with 1s where fixIndOlayIm is 0 (no fixation marker/inducers) and 0s where fixIndOlayIm is non-zero.

    %% create some trial properties we need later

    trial_tmp = [1, 2];
    trial_prob = t_prob;

    trialprop = cell(numtrials,5);
    rng('shuffle')
    for x = 1:numtrials
        %rng('shuffle')
        tcon_idx = randsample(trial_tmp,1,true,[(1-trial_prob) trial_prob]);
        trialprop(x,1) = {trial_con(tcon_idx)};   % orientation EJD Added ,1 to select v or h...

        %rng('shuffle')
        tcon_idx = randsample(trial_tmp,1,true,[(1-trial_prob) trial_prob]);
        trialprop(x,2) = {trial_con(tcon_idx)};   % orientation EJD Added ,1 to select v or h...

        %rng('shuffle')
        tcon_idx = randsample(trial_tmp,1,true,[(1-trial_prob) trial_prob]);
        trialprop(x,3) = {trial_con(tcon_idx)};   % orientation EJD Added ,1 to select v or h...

        %rng('shuffle')
        tcon_idx = randsample(trial_tmp,1,true,[(1-trial_prob) trial_prob]);
        trialprop(x,4) = {trial_con(tcon_idx)};   % orientation EJD Added ,1 to select v or h...

        %rng('shuffle')
        tcon_idx = randsample(trial_tmp,1,true,[(1-trial_prob) trial_prob]);
        trialprop(x,5) = {trial_con(tcon_idx)};   % orientation EJD Added ,1 to select v or h...
    end

    % for x = 1:numtrials
    %     trialprop(x,1) = {randsample(trial_con,1)};   % orientation EJD Added ,1 to select v or h...
    %     trialprop(x,2) = {randsample(trial_con,1)};   % orientation EJD Added ,1 to select v or h...
    %     trialprop(x,3) = {randsample(trial_con,1)};   % orientation EJD Added ,1 to select v or h...
    %     trialprop(x,4) = {randsample(trial_con,1)};   % orientation EJD Added ,1 to select v or h...
    %     trialprop(x,5) = {randsample(trial_con,1)};   % orientation EJD Added ,1 to select v or h...
    % end

    %ifi = Screen('GetFlipInterval', window, 100);
    % CLUT stuff - make sure it matches what's in TMSfly2.m
    white = WhiteIndex(window);
    black = BlackIndex(window);
    grayblack = floor(GrayIndex(window,0.15));
    chargray = ceil(GrayIndex(window,0.30));
    newgray = ceil(GrayIndex(window,0.35));
    lgray = ceil(GrayIndex(window,0.75)); % light gray
    dgray = floor(GrayIndex(window,0.5)); % dark grayFlip
    
    %gray = [146, 145, 145, 255];	% medium gray
    gray = white / 2;
    
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % a little bit of texture generation for the displaysFlip
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % open an offscreen window for background
    backwindow_fix=Screen('OpenOffscreenWindow', window, gray);
    backwindow=Screen('OpenOffscreenWindow', window, gray);
    % another one for noise
    if noise
        backwindow_n=Screen('OpenOffscreenWindow', window, gray);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ok, now we're ready for some threshold estimation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % some constants
    % stimtime = .125; % length of time (in seconds) that target is on the screen
    % response = zeros(numtrials,1);  % pre-allocate
    beephi = MakeBeep(2000,.1);     % correct beep, for feedback
    beeplo = MakeBeep(500,.1);      % incorrect beep, for feedback

    %% Setup Quest first (depends on linking cue)
    % common to all
    beta=3.5;delta=0.01;gamma=0.5;
    switch link
        case 'F' %
            if n_int >= 1 
            % Staircase1
            tGuess1 = qst_startParams.q1.tGuess;
            tGuessSd1 = qst_startParams.q1.tGuessSd;  % std. dev. on initial guess of threshold
            grain1 = qst_startParams.q1.grain;
            range1 = qst_startParams.q1.range;
            pThreshold1=pthrds.thr1;

            % Delete q1 if it exists..
            if exist("q1","var")
                clear(q1);
            end

            % Set up the quest
            q1=QuestCreate(tGuess1,tGuessSd1,pThreshold1,beta,delta,gamma,grain1,range1);

            q1.normalizePdf=1; % This adds a few ms per call to QuestUpdate, but otherwise the pdf will underflow after about 1000 trials.

            end

            if n_int >= 2
            %Staircase2
            tGuess2 = qst_startParams.q2.tGuess;
            tGuessSd2 = qst_startParams.q2.tGuessSd;  % std. dev. on initial guess of threshold
            grain2 = qst_startParams.q2.grain;
            range2 = qst_startParams.q2.range;
            pThreshold2=pthrds.thr2;

            % Delete q2 if it exists..
            if exist("q2","var")
                clear(q2);
            end

            % Set up the quest
            q2=QuestCreate(tGuess2,tGuessSd2,pThreshold2,beta,delta,gamma,grain2,range2);

            q2.normalizePdf=1; % This adds a few ms per call to QuestUpdate, but otherwise the pdf will underflow after about 1000 trials.

            end

            if n_int >= 3
            %Staircase3
            tGuess3 = qst_startParams.q3.tGuess;
            tGuessSd3 = qst_startParams.q3.tGuessSd;  % std. dev. on initial guess of threshold
            grain3 = qst_startParams.q3.grain;
            range3 = qst_startParams.q3.range;
            pThreshold3 = pthrds.thr3;

            % Delete q3 if it exists..
            if exist("q3","var")
                clear(q3);
            end

            % Set up the quest
            q3=QuestCreate(tGuess3,tGuessSd3,pThreshold3,beta,delta,gamma,grain3,range3);

            q3.normalizePdf=1; % This adds a few ms per call to QuestUpdate, but otherwise the pdf will underflow after about 1000 trials.

            end

            %Staircase4
            if n_int >= 4
            tGuess4 = qst_startParams.q4.tGuess;
            tGuessSd4 = qst_startParams.q4.tGuessSd;  % std. dev. on initial guess of threshold
            grain4 = qst_startParams.q4.grain;
            range4 = qst_startParams.q4.range;
            pThreshold4=pthrds.thr4;

            % Delete q4 if it exists..
            if exist("q4","var")
                clear(q4);
            end

            % Set up the quest
            q4=QuestCreate(tGuess4,tGuessSd4,pThreshold4,beta,delta,gamma,grain4,range4);

            q4.normalizePdf=1; % This adds a few ms per call to QuestUpdate, but otherwise the pdf will underflow after about 1000 trials.

            end

            %Staircase5
            if n_int >= 5
            tGuess5 = qst_startParams.q5.tGuess;
            tGuessSd5 = qst_startParams.q5.tGuessSd;  % std. dev. on initial guess of threshold
            grain5 = qst_startParams.q5.grain;
            range5 = qst_startParams.q5.range;
            pThreshold5=pthrds.thr5;

            % Delete q5 if it exists..
            if exist("q5","var")
                clear(q5);
            end

            % Set up the quest
            q5=QuestCreate(tGuess5,tGuessSd5,pThreshold5,beta,delta,gamma,grain5,range5);

            q5.normalizePdf=1; % This adds a few ms per call to QuestUpdate, but otherwise the pdf will underflow after about 1000 trials.
            
            end

            % Get the expected nr_sigmoid x value to give you the desired
            % pThresholds..... ie what param1 do you need to feed in to get
            % 0.7 accuracy as an output?

            syms c
            %r(c) = nr_rmax*((c^n)/(nr_c50^n + c^n)) + nr_b; % naka rushton model function
            expected_param1_value = double(max(solve(nr_rmax*((c^nr_n)/(nr_c50^nr_n + c^nr_n)) + nr_b == pThreshold1,c)));
            if n_int >= 1
            nr_ideal_t1 = double(max(solve(nr_rmax*((c^nr_n)/(nr_c50^nr_n + c^nr_n)) + nr_b == pThreshold1,c)));
            else
            nr_ideal_t1=[];
            end
            
            if n_int >= 2
            nr_ideal_t2 = double(max(solve(nr_rmax*((c^nr_n)/(nr_c50^nr_n + c^nr_n)) + nr_b == pThreshold2,c)));
            else
            nr_ideal_t2=[];
            end
            
            if n_int >= 3
            nr_ideal_t3 = double(max(solve(nr_rmax*((c^nr_n)/(nr_c50^nr_n + c^nr_n)) + nr_b == pThreshold3,c)));
            else
            nr_ideal_t3=[];
            end
            
            if n_int >= 4
            nr_ideal_t4 = double(max(solve(nr_rmax*((c^nr_n)/(nr_c50^nr_n + c^nr_n)) + nr_b == pThreshold4,c)));
            else
            nr_ideal_t4=[];
            end
            
            if n_int >= 5
            nr_ideal_t5 = double(max(solve(nr_rmax*((c^nr_n)/(nr_c50^nr_n + c^nr_n)) + nr_b == pThreshold5,c)));
            else
            nr_ideal_t5=[];
            end
            
            % Get the dimensions of the offscreen window as an image
            %Screen('FillRect', backwindow, grey);  % black out the offscreen window
            im_4_WinDims = Screen('GetImage', window); %Save the window as an image
            [s1, s2, s3] = size(im_4_WinDims); %save it's dimensions..

            % LOAD IN LIST OF APE-HUMAN FACE IMAGES AHEAD OF TRIALS..
            % Go to the directory where the images live...
            cd(image_dir)
            %Check whether occluded or unoccluded versions were requested..
            % select accordingly..

            % get a list of the images whose names contain the tag/string "imtag"
            if olCon == 0
                % just occ
                cd("occ") %NOTE..MAY WANT TO EVENTUALLY GENERALIZE..
                [imFiles, numFilesProcessed] = getAllFileNamesRec(pwd,imtag);
                cd(image_dir);
            elseif olCon == 1
                % just nocc
                cd("nocc") %NOTE..MAY WANT TO EVENTUALLY GENERALIZE..
                [imFiles, numFilesProcessed] = getAllFileNamesRec(pwd,imtag);
                cd(image_dir);
            elseif olCon == 2
                % both..
                cd("occ") %NOTE..MAY WANT TO EVENTUALLY GENERALIZE..
                [imFiles1, numFilesProcessed1] = getAllFileNamesRec(pwd,imtag);
                cd(image_dir);
                cd("nocc") %NOTE..MAY WANT TO EVENTUALLY GENERALIZE..
                [imFiles2, numFilesProcessed2] = getAllFileNamesRec(pwd,imtag);
                cd(image_dir);
                % Combine occ and nocc outputs into same arrays..
                numFilesProcessed=numFilesProcessed1+numFilesProcessed2;
                imFiles=vertcat(imFiles1,imFiles2);
            end

            fprintf('We found %d files ...\n', ...
                numFilesProcessed);

            % Now read in the images contained in imFiles
            %img_array = cell(size(imFiles,1),3);
            img_array = cell(size(imFiles,1),2);
            for ii=1:size(imFiles,1)
                fname = imFiles{ii,1}; % get filename
                img = imread(fname);    % read in image
                img_array(ii,1) = {img};  % save image in img_array cell array column1..
                img_array(ii,2) = {fname}; % save full image path/name in img_array cell array column2..

                % Extract the parameter values for tags in in 'parTagz' with
                % xtractParFrmFilename.m
                [parValz, parTagz] = xtractParsFrmFilename(fname,parTagz,"num");
                ori_con = xtractParsFrmFilename(fname,["ori"],"str");%orientation condition
                [occ_con, ~] = xtractParsFrmFilename(fname,"occ","num");%occlusion condition
                countr=1; %initialize countr
                for j = 3:(3+(length(parValz)-1))
                    img_array(ii,j) = {parValz(1,countr)};  % save param in
                    % img_array cell array column3.
                    countr=countr+1;
                end
                clear countr
                img_array(ii,(3+(length(parValz)-1))+1)={ori_con};
                img_array(ii,(3+(length(parValz)-1))+2)={occ_con};
            end
            clear ii

            if noise==1
                % READ IN THE NOISE IMAGES
                %nTag=".png";
                [noiseFiles, nNoizeFiles] = getAllFileNamesRec(noise_dir,nTag);

                % Check subject directory to see whether the usedNoise.mat
                % file exists. If so, read in the list of frames that have been
                % used in previous sessions and take them off the list for this
                % session to prevent any given frame from being used more than
                % once per subject.
                %--------------------------------------------------------------
                %--------------------------------------------------------------
                strtDirTmp=pwd;% save start position
                cd(sessionDir);% go to the output dir for this session
                % now go up 4 levels..
                % this should be the main subject directory
                cd ..
                cd ..
                cd ..
                cd ..

                % Check whether the usedNoise.mat file exists..
                usedNoiseFileName="usedNoise.mat";
                if isfile(usedNoiseFileName)
                    isUsedNoiseFile=1;
                    disp(" ")
                    disp(strcat(usedNoiseFileName," found. Reading in list of noise images used in previous sessions and removing them from the list for this session.."));
                    disp(" ")
                    usedNoise=load(usedNoiseFileName);
                    usedFrmList=usedNoise.usedFrmList;
                    

                else
                    isUsedNoiseFile=0;
                    disp(" ")
                    disp(strcat(usedNoiseFileName," not found. Assuming this is the first run with this subject/setup..."));
                    disp(" ")
                end
                cd(strtDirTmp); % return to start point..
                %--------------------------------------------------------------
                %--------------------------------------------------------------

                % Now load in the image file paths contained in imFiles
                % into a cell array.
                %img_array = cell(size(imFiles,1),3);
                noise_array = cell(size(noiseFiles,1),2);
                noise_arrayIdx=zeros(size(noiseFiles,1),1);
                for ii=1:size(noiseFiles,1)
                    fname = noiseFiles{ii,1}; % get filename
                    %img = imread(fname);    % read in image
                    %noise_array(ii,1) = {};  % leave idx 1 (where we will later save image} blank
                    noise_array(ii,2) = {fname}; % save full image path/name in img_array cell array column2..
                    noise_arrayIdx(ii,1)=ii;
                end
                clear ii

                % if Used Noise File was present cross-check against the
                % used noise file list and eliminate noise images that were
                % used previously..
                if isUsedNoiseFile==1
                    
                    % Generalize paths for comparison 
                    % (cut off the "main_dir" (everything before "revCorrStim")
                    % string from the beginning of each in case some trials 
                    % were collected on other machines...)
                    nzArryGenzd = generalizePathz('revCorrStim',noise_array(:,2));
                    usedFrmzGenzd = generalizePathz("revCorrStim",usedFrmList);

                    [~, ~, ~, NonMatchIndexs] = compareCellStrLists(nzArryGenzd,usedFrmzGenzd);
                    %[~, ~, ~, NonMatchIndexs] = compareCellStrLists(noise_array(:,2),usedFrmList);
                    
                    noise_arrayNonMatch=cell(size(NonMatchIndexs,1),2);
                    noise_arrayIdx=zeros(size(noise_arrayNonMatch,1),1);
                    for pp=1:size(NonMatchIndexs,1)
                        TmpIdx=NonMatchIndexs{pp,1};
                        noise_arrayNonMatch{pp,1}=noise_array{TmpIdx,1};
                        noise_arrayNonMatch{pp,2}=noise_array{TmpIdx,2};
                        noise_arrayIdx(pp,1)=pp;    
                    end

                    % set noise_array equal to noise_arrayNonMatch
                    noise_array=noise_arrayNonMatch;
                    % remove noise_arrayNonMatch/other variables..
                    clear noise_arrayNonMatch;
                    clear nzArryGenzd
                    clear usedFrmzGenzd
                    
                end

                % randomize noise_arrayIdx row order
                % and grab the top ntrials indices
                % to use for the trials in this session..
                noise_arrayIdx = noise_arrayIdx(randperm(size(noise_arrayIdx, 1)), :);
                totalTrials=n_int*numtrials;

                if size(noise_arrayIdx,1)>=totalTrials

                    noise_arrayIdx=noise_arrayIdx(1:totalTrials,:);

                    % finally, reduce noise_array to only the indices contained
                    % in the random index selection above. While doing this read
                    % the noise files whose paths are in the selected
                    % row indices above in as images and save in their
                    % respective cell positions.
                    noise_arrayTmp=cell(size(noise_arrayIdx,1),size(noise_array,2)); % preallocate
                    for ii=1:size(noise_arrayTmp,1)
                        noise_arrayTmp(ii,:)=noise_array(noise_arrayIdx(ii,1),:); % grab the row
                        noise_arrayTmp(ii,1)={imread(noise_arrayTmp{ii,2})}; %read in the image filepath in col 2 and save in col1
                    end
                    % set noise_array variable equal to limited/randomized
                    % noise_arrayTmp set above and also
                    % remove the temp variable noise_arrayTmp
                    noise_array=noise_arrayTmp;
                    clear noise_arrayTmp;

                else
                    disp("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                    disp("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                    disp("!!! THIS SUBJECT HAS USED ALL AVAILABLE NOISE !!!");
                    disp("!!!    FRAMES CREATED FOR THIS EXPERIMENT     !!!");
                    disp("!!!        TIME TO COOK UP SOME MORE          !!!");
                    disp("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                    disp("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                    
                    assert(size(noise_arrayIdx,1)>=totalTrials,'Out of unique noise frames!');
                end
            end

            %EJD CURRENTLY HERE ***4/7/23
            % FOR LUM1/LUM2 ... MAY NEED TO UPDATE TO DEAL WITH OTHER PAR
            % CONFIGURATIONS..
            %#############################################################
            %par_inxs = cell2mat(img_array(1:end,3:4)); % save the parameters column for index referencing later (2:end bc we don't want to include the "Ape only" image at 1)


            % (Difference Parameter QUESTS..)
            % For Luminance
            if qstParStr == "lumDiff"
                whichVersion='lum';
                % compute absolute luminance difference..
                difArray = (abs(cell2mat(img_array(1:end,3))-cell2mat(img_array(1:end,4))));
                img_array(:,end+1) = num2cell(difArray); % add the difference values to img_array
            end

            % For Texture
            if qstParStr == "texDiff"
                whichVersion='tex';
                % For textures the "difference" is a scaling difference..
                % Compute relative scaling difference: (scaleFactor1/scaleFactor2)
                % difArray outputs will therefore represent the multiplicative
                % scaling factor between the textures windowed withing the
                % objects
                difArray = ((cell2mat(img_array(1:end,3)))./(cell2mat(img_array(1:end,4))));
                img_array(:,end+1) = num2cell(difArray); % add the difference values to img_array
            end


            % For Common Region (Individual Object Parameter QUESTS..)
            if qstParStr == "CRpl"
                whichVersion='cr';
                difArray = (abs(cell2mat(img_array(1:end,3))-gray)); % get difference between background and cr box luminance...
                img_array(:,end+1) = num2cell(difArray);
            end


            % Create a logical index that is true for rows with each orientation
            idx1 = strcmp([img_array{:,2+length(parValz)+1}], [trial_con(1)]); %for orientation 1
            idx2 = strcmp([img_array{:,2+length(parValz)+1}], trial_con(2)); %for orientation 2
            
            % Use the logical index to select only the rows with each orientation
            img_arrayC1 = img_array(idx1, :);
            img_arrayC2 = img_array(idx2, :);
            
            if qstParStr == "lumDiff" || qstParStr == "texDiff"
            %par_inxs = horzcat(par_inxs,difArray);
            par_inxsC1 = cell2mat(img_arrayC1(1:end,3:4));
            par_inxsC1 = horzcat(par_inxsC1,cell2mat(img_arrayC1(1:end,size(img_arrayC1,2)))); % add on the difference column/QUEST parameter column

            par_inxsC2 = cell2mat(img_arrayC2(1:end,3:4));
            par_inxsC2 = horzcat(par_inxsC2,cell2mat(img_arrayC2(1:end,size(img_arrayC1,2)))); % add on the difference column/QUEST parameter column
            end

            if qstParStr == "CRpl"
            %par_inxs = horzcat(par_inxs,difArray);
            par_inxsC1 = cell2mat(img_arrayC1(1:end,3));
            par_inxsC1 = horzcat(par_inxsC1,cell2mat(img_arrayC1(1:end,size(img_arrayC1,2)))); % add on the difference column/QUEST parameter column

            par_inxsC2 = cell2mat(img_arrayC2(1:end,3));
            par_inxsC2 = horzcat(par_inxsC2,cell2mat(img_arrayC2(1:end,size(img_arrayC1,2)))); % add on the difference column/QUEST parameter column
            end
            %#############################################################

            % set min and max params..
            maxparam = maxmin{1,1};
            minparam = maxmin{1,2};
    end

    if giveInstrxn == 1
    %revCorrInstructions_v1(window,screenrect,intro_txtHeaderSize,intro_txtSize,scrn_top,scrn_1percY,black,gray,whichVersion);
    revCorrInstructions_v1(window,screenrect,intro_txtHeaderSize,intro_txtSize,scrn_top,scrn_bot,scrn_1percY,scrn_1percX,scrn_left,scrn_right,black,gray,whichVersion);

    end

    %% Start trials and add results to Quest databases
    %Screen('FillRect', window, [grayblack grayblack grayblack 0]);  % blank out the screen
    Screen('FillRect', window, gray); % gray out window
    
    trialtm = zeros(5,numtrials);   % pre-allocate
    Screen('TextSize',window, intro_txtSize);
    speech1 = ' Press ENTER/RETURN to begin presentation of stimulus.';

%     if db_mode == 0
%         %reply = Ask(window, speech1, red, gray, ['KbWait(' num2str(kboard) ')'], RectLeft, RectTop, 15);
%         reply = Ask(window,speech1,red,gray, ['KbWait(' num2str(kboard) ')'],RectLeft,RectTop,intro_txtSize);
%     end

    reply = Ask(window,speech1,red,gray,'GetString',RectLeft,RectTop,intro_txtSize);
    
    %Screen('FillRect', window, [grayblack grayblack grayblack 0]);  % blank out the screen
    Screen('FillRect', window, gray); % gray out window

    % put fixation display up
    %[backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack);
    %response = zeros(numtrials,5);

    response = zeros(numtrials,5);

    %Start KbQueue
    KbQueueCreate();

    % preallocate space for responses
    %corr_resp = string(1,numtrials);
    im_mat = zeros(numtrials,5);
    noizMat= zeros(numtrials,5);
    im_mat = string(im_mat);
    noizMat= string(noizMat);

    con_mat = zeros(numtrials,5);
    con_mat = string(con_mat);
    resp_mat = zeros(numtrials,5);
    resp_mat = string(resp_mat);
    corr_resp_mat = zeros(numtrials,5);
    corr_resp_mat = string(resp_mat);
    tim_mat = zeros(numtrials,3,5);
    tim_mat2 = zeros(numtrials,4,5);
    ImParMat = zeros(numtrials,5);
    ObjParMat = zeros(numtrials,length(parTagz),5);
    %occCon_mat = zeros(numtrials,5);

    nzFrmCountr=1; % initialize variable for counting noise frames/tracking index..

    for trials = 1:numtrials
        %% Staircase 1

        if link=='L'
            param1(trials)=round(10^(QuestQuantile(q1)));   % get next level of parameter to test
        elseif link=='B'
            param1(trials)=quant(10^(QuestQuantile(q1)),0.1);  % get next level of parameter to test
        elseif link=='F'
            param1(trials)=quant(10^(QuestQuantile(q1)),0.1);  % get next level of parameter to test
        end
        if isnan(param1(trials))
            param1(trials)= param1(trials-1);
        end

        %EJD COMMENTED/EDITED BELOW
        %param1(trials) = max(minparam,min(maxparam,param1(trials)));% properly restrict value to useable range

        %if maxmin_seg, set up the max-min segments
        if maxmin_seg == 1

            if n_int == 2
                seg1min=minparam;
                seg1max=maxparam/2;
                seg2min=maxparam/2;
                seg2max=maxparam;
            end

            if n_int == 3
                seg1min=minparam;
                seg1max=maxparam/3;
                seg2min=maxparam/3;
                seg2max=2*(maxparam/3);
                seg3min=2*(maxparam/3);
                seg3max=maxparam;
            end

            if n_int == 4
                seg1min=minparam;
                seg1max=maxparam/4;
                seg2min=maxparam/4;
                seg2max=2*(maxparam/4);
                seg3min=2*(maxparam/4);
                seg3max=3*(maxparam/4);
                seg4min=3*(maxparam/4);
                seg4max=maxparam;
            end

            if n_int == 5
                seg1min=minparam;
                seg1max=maxparam/5;
                seg2min=maxparam/5;
                seg2max=2*(maxparam/5);
                seg3min=2*(maxparam/5);
                seg3max=3*(maxparam/5);
                seg4min=3*(maxparam/5);
                seg4max=4*(maxparam/5);
                seg5min=4*(maxparam/5);
                seg5max=maxparam;
            end

        end

        if maxmin_seg == 1
            % limit to range btw maxparam and minparam
            if param1(trials) > seg1max
                param1(trials) = seg1max;
            end
            if param1(trials) < seg1min
                param1(trials) = seg1min;
            end
        end

        if maxmin_seg == 0
            % limit to range btw maxparam and minparam
            if param1(trials) > maxparam
                param1(trials) = maxparam;
            end
            if param1(trials) < minparam
                param1(trials) = minparam;
            end
        end

        switch link
            case 'F'
                tcon = trialprop{trials,1}; % EJD NOTE: was "orient" now.. "face"...
                corr_resp_mat(trials,1) = tcon;
                [backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack,fix_im,window);
                Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed
                Screen('Flip',window);
                strt = GetSecs;
                tim_mat2(trials,1,1) = strt;
                
                if tcon==trial_con(2)
                    % grab the image from  img_array with parameter closest to one
                    % suggested by from quest (param1)
                    a=par_inxsC2(:,end)';
                    n=param1(end); % assign the most recent param in param1 to "n"
                    [val,idx]=min(abs(a-n));
                    closest_param=a(idx);
                    % Check if there are multiple images which give this
                    % value.. if so pick one at random..
                    [row, col] = find(a == closest_param);
                    idxTmp = randi([1, length(col)]);
                    idx=col(idxTmp);

                    param1(trials)=closest_param;  % get next level of parameter to test
                    con_mat(trials,1) = "Quest1";
                    im_mat(trials,1) = img_arrayC2{idx,2}; % Save image name/path
                    ImParMat(trials,1) = closest_param;
                    %Get the object parameters at index idx too..
                    ObjParMat(trials,1:end,1)=par_inxsC2(idx,1:length(parTagz));
                    %occCon_mat(trials,1)=img_arrayC2{idx,6}; % save the occlusion condition
                    %!!??
                    image_in = img_arrayC2{idx,1};
                end

               
                if tcon==trial_con(1)
                    % grab the image from  img_array with parameter closest to one
                    % suggested by from quest (param1)
                    a=par_inxsC1(:,end)';
                    n=param1(end); % assign the most recent param in param1 to "n"
                    [val,idx]=min(abs(a-n));
                    closest_param=a(idx);
                    % Check if there are multiple images which give this
                    % value.. if so pick one at random..
                    [row, col] = find(a == closest_param);
                    idxTmp = randi([1, length(col)]);
                    idx=col(idxTmp);

                    param1(trials)=closest_param;  % get next level of parameter to test
                    con_mat(trials,1) = "Quest1";
                    im_mat(trials,1) = img_arrayC1{idx,2}; % Save image name/path
                    ImParMat(trials,1) = closest_param;
                    %Get the object parameters at index idx too..
                    ObjParMat(trials,1:end,1)=par_inxsC1(idx,1:length(parTagz));
                    %occCon_mat(trials,1)=img_arrayC1{idx,6}; % save the occlusion condition
                    %!!??
                    image_in = img_arrayC1{idx,1};
                end
                

                if noise==1
%                    % create a random permutation of the noise row indices
%                     NoiseRanIdxOrder= randperm(size(noise_arrayIdx, 1));
%                     % reorder the rows of the noise_arrayIdx matrix using the permutation
%                     noise_arrayIdx = noise_arrayIdx(NoiseRanIdxOrder, :);
%                     noiseImgIn = noise_array{noise_arrayIdx(1,1),1}; %get the noise image for the index on top..
                    
                    
                    noiseImgIn = noise_array{nzFrmCountr,1};
                    [image_in] = revCorrAddNoise_v1(image_in,noiseImgIn,nWt,fixIndOlayIm,nonFxnIndMask); % add the noise to the image..
                    noizMat(trials,1) = noise_array{nzFrmCountr,2};
                    %noise_arrayIdx = noise_arrayIdx(2:end,1); % now eliminate that index from the list so the same noise image won't be selected again..
                    nzFrmCountr=nzFrmCountr+1; % update countr..
                end

                % draw linking display to offscreen window
                [backwindow] = PORClinkRevCor_v1(mainsrc, window, backwindow,trialprop(trials,1), trial_con, image_in,fixrect,grayblack);
                %Screen('DrawTexture', window, backwindow, mainsrc, maindest, 0, 0, 1);  % draw stimulus to the on-screen window, rotated & zoomed
                Screen('DrawTexture', window, backwindow, mainsrc, maindest, 0, 0, 1);  % draw stimulus to the on-screen window, rotated & zoomed
                %[backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack,fix_im);
                %Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed

                keyIsDown = 0;
                swtch = 0;
                swtch2 = 0;
                respns = "";
                %t_trial0 = GetSecs;
                % While stimulus is on the screen and no response
                while ((GetSecs-strt) <= presTime2) && keyIsDown == 0 %trialResp == 0

                    % Draw stimulus
                    if ((GetSecs-strt)>=presTime) && swtch == 0
                        Screen('DrawTexture', window, backwindow, mainsrc, maindest, 0, 0, 1);  % draw stimulus to the on-screen window, rotated & zoomed
                        Screen('Flip', window);% show linking display
                        t_disp = GetSecs; % record time display came up..
                        %start listening for responses
                        KbQueueStart();
                        swtch = 1;
                    end

                    [keyIsDown, firstPress, ~, ~, ~] = KbQueueCheck();
                end
                % Check to see if we exited the While loop above due to a
                % response.. if not, raise swtch2 to indicate we still don't have a response..
                if keyIsDown == 0
                    swtch2 = 1;
                end

                % remove stimulus image/swap for just fixation..
                [backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack, fix_im,window);
                Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed
                Screen('Flip',window);
                %[blank_time StimulusOnsetTime FlipTimestamp Missed Beampos] = Screen('Flip', window, link_time+0.4-(0.5*ifi));%blank comes 400ms after linking display


                % NOTE!!!!!: ASK ABOUT WHETHER WE NEED TO WAIT UNTIL THE MAX RESPONSE WINDOW TO PROCEED. IF SO: GET RID OF THIS IF STATEMENT. IF NOT, KEEP AND GET RID OF IF IN WHILE LOOP BELOW AS ITS REDUNDANT..

                % wait until the response time cutoff is reached if there has not been a response yet..
                while ((GetSecs-strt) < resp_cut) && keyIsDown == 0
                    [keyIsDown, firstPress, ~, ~, ~] = KbQueueCheck();
                end
                if (keyIsDown == 1) && (swtch2 == 1)
                    t_resp = max(firstPress)-t_disp;
                    swtch2 = 0; % turn off switch bc a response occured
                end

                %stop listening for responses
                KbQueueStop();
                tim_mat(trials,2,1) = GetSecs-strt; % save total trial time..
                tim_mat2(trials,2,1) = GetSecs;
                % Record Response
                if swtch2 == 0
                    tim_mat(trials,1,1) = max(firstPress)-t_disp; % record time of response
                    resp_mat(trials,1) = KbName(firstPress); % record button press
                    respns = KbName(firstPress);
                else
                    tim_mat(trials,1,1) = -999; % record time of response
                    resp_mat(trials,1) = "!"; % record time of response
                end

        end


        % If feedback is on, provide the appropriate feedback
        switch fback_On
            case 1
                if respns == tcon   % vert key pressed
                    response(trials,1) = 1;
                    [backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack,fix_im,window);
                    answer = feedback(1);
                    %answer = strcat(answer," ","button pressed:", tdesc_file{trials,9});
                    tfback0 = GetSecs;
                    while (GetSecs-tfback0)<=feedbackTime
                        Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed
                        Screen('TextSize',window, fback_txtSize);
                        DrawFormattedText(window, answer,'center', scrn_top+30*scrn_1percY, black);
                        Screen('Flip',window);
                    end
                else
                    response(trials,1) = 0;
                    [backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack,fix_im,window);
                    answer = feedback(2);
                    %answer = strcat(answer,' ','button pressed:', tdesc_file{trials,9});
                    tfback0 = GetSecs;
                    while (GetSecs-tfback0)<=feedbackTime
                        Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed
                        Screen('TextSize',window, fback_txtSize);
                        DrawFormattedText(window, answer,'center', scrn_top+30*scrn_1percY, black);
                        Screen('Flip',window);
                    end
                end
            case 0
                if respns == tcon
                    response(trials,1) = 1;
                else
                    response(trials,1) = 0;
                end
        end

        % Inter-trial Interval..
        [backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack,fix_im,window);
        iti_0 = GetSecs;
        while (GetSecs-iti_0)<ITI
            Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed
            Screen('Flip',window);
        end

        tim_mat2(trials,4,1) = GetSecs;
        tim_mat2(trials,3,1) = iti_0;


        tim_mat(trials,3,1) = tim_mat2(trials,4,1)-strt; % save total trial time after ITI too..

        % Flush the Response Queue..
        KbQueueFlush();

        % Add the new datum
        q1 = QuestUpdate(q1,log10(param1(trials)),response(trials,1));

        %% Staircase 2
        if n_int>1
            if link=='L'
                param2(trials)=round(10^(QuestQuantile(q2)));   % get next level of parameter to test
            elseif link=='B'
                param2(trials)=quant(10^(QuestQuantile(q2)),0.1);  % get next level of parameter to test
            elseif link=='F'
                param2(trials)=quant(10^(QuestQuantile(q2)),0.1);  % get next level of parameter to test
            end
            if isnan(param2(trials))
                param2(trials)= param2(trials-1);
            end

            %EJD COMMENTED/EDITED BELOW
            %param1(trials) = max(minparam,min(maxparam,param1(trials)));% properly restrict value to useable range

            if maxmin_seg == 1
                % limit to range btw maxparam and minparam
                if param2(trials) > seg2max
                    param2(trials) = seg2max;
                end
                if param2(trials) < seg2min
                    param2(trials) = seg2min;
                end
            end

            if maxmin_seg == 0
                % limit to range btw maxparam and minparam
                if param2(trials) > maxparam
                    param2(trials) = maxparam;
                end
                if param2(trials) < minparam
                    param2(trials) = minparam;
                end
            end

            switch link
                case 'F'
                    tcon = trialprop{trials,2}; % EJD NOTE: was "orient" now.. "face"...
                    corr_resp_mat(trials,2) = tcon;
                    [backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack,fix_im,window);
                    Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed
                    Screen('Flip',window);
                    strt = GetSecs;
                    tim_mat2(trials,1,2) = strt;

                if tcon==trial_con(2)
                    % grab the image from  img_array with parameter closest to one
                    % suggested by from quest (param1)
                    a=par_inxsC2(:,end)';
                    n=param2(end); % assign the most recent param in param1 to "n"
                    [val,idx]=min(abs(a-n));
                    closest_param=a(idx);
                    % Check if there are multiple images which give this
                    % value.. if so pick one at random..
                    [row, col] = find(a == closest_param);
                    idxTmp = randi([1, length(col)]);
                    idx=col(idxTmp);

                    param2(trials)=closest_param;  % get next level of parameter to test
                    con_mat(trials,2) = "Quest2";
                    im_mat(trials,2) = img_arrayC2{idx,2}; % Save image name/path
                    ImParMat(trials,2) = closest_param;
                    %Get the object parameters at index idx too..
                    ObjParMat(trials,1:end,2)=par_inxsC2(idx,1:length(parTagz));
                    %occCon_mat(trials,2)=img_arrayC2{idx,6}; % save the occlusion condition
                    %!!??
                    image_in = img_arrayC2{idx,1};
                end

               
                if tcon==trial_con(1)
                    % grab the image from  img_array with parameter closest to one
                    % suggested by from quest (param1)
                    a=par_inxsC1(:,end)';
                    n=param2(end); % assign the most recent param in param1 to "n"
                    [val,idx]=min(abs(a-n));
                    closest_param=a(idx);
                    % Check if there are multiple images which give this
                    % value.. if so pick one at random..
                    [row, col] = find(a == closest_param);
                    idxTmp = randi([1, length(col)]);
                    idx=col(idxTmp);

                    param2(trials)=closest_param;  % get next level of parameter to test
                    con_mat(trials,2) = "Quest2";
                    im_mat(trials,2) = img_arrayC1{idx,2}; % Save image name/path
                    ImParMat(trials,2) = closest_param;
                    %Get the object parameters at index idx too..
                    ObjParMat(trials,1:end,2)=par_inxsC1(idx,1:length(parTagz));
                    %occCon_mat(trials,2)=img_arrayC1{idx,6}; % save the occlusion condition
                    %!!??
                    image_in = img_arrayC1{idx,1};
                end

                if noise==1
%                     % create a random permutation of the noise row indices
%                     NoiseRanIdxOrder= randperm(size(noise_arrayIdx, 1));
%                     % reorder the rows of the noise_arrayIdx matrix using the permutation
%                     noise_arrayIdx = noise_arrayIdx(NoiseRanIdxOrder, :);
%                     noiseImgIn = noise_array{noise_arrayIdx(1,1),1}; %get the noise image for the index on top..
                    
                    noiseImgIn = noise_array{nzFrmCountr,1};

                    [image_in] = revCorrAddNoise_v1(image_in,noiseImgIn,nWt,fixIndOlayIm,nonFxnIndMask); % add the noise to the image..
                    noizMat(trials,2) = noise_array{nzFrmCountr,2};
                    %noise_arrayIdx = noise_arrayIdx(2:end,1); % now eliminate that index from the list so the same noise image won't be selected again..
                    nzFrmCountr=nzFrmCountr+1; % update countr..
                end

                    % draw linking display to offscreen window
                    [backwindow] = PORClinkRevCor_v1(mainsrc, window, backwindow,trialprop(trials,1), trial_con, image_in,fixrect,grayblack);
                    %Screen('DrawTexture', window, backwindow, mainsrc, maindest, 0, 0, 1);  % draw stimulus to the on-screen window, rotated & zoomed
                    Screen('DrawTexture', window, backwindow, mainsrc, maindest, 0, 0, 1);  % draw stimulus to the on-screen window, rotated & zoomed
                    %[backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack,fix_im);
                    %Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed

                    keyIsDown = 0;
                    swtch = 0;
                    swtch2 = 0;
                    respns = "";
                    %t_trial0 = GetSecs;
                    % While stimulus is on the screen and no response
                    while ((GetSecs-strt) <= presTime2) && keyIsDown == 0 %trialResp == 0

                        % Draw stimulus
                        if ((GetSecs-strt)>=presTime) && swtch == 0
                            Screen('DrawTexture', window, backwindow, mainsrc, maindest, 0, 0, 1);  % draw stimulus to the on-screen window, rotated & zoomed
                            Screen('Flip', window);% show linking display
                            t_disp = GetSecs; % record time display came up..
                            %start listening for responses
                            KbQueueStart();
                            swtch = 1;
                        end

                        [keyIsDown, firstPress, ~, ~, ~] = KbQueueCheck();
                    end
                    % Check to see if we exited the While loop above due to a
                    % response.. if not, raise swtch2 to indicate we still don't have a response..
                    if keyIsDown == 0
                        swtch2 = 1;
                    end

                    % remove stimulus image/swap for just fixation..
                    [backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack, fix_im,window);
                    Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed
                    Screen('Flip',window);
                    %[blank_time StimulusOnsetTime FlipTimestamp Missed Beampos] = Screen('Flip', window, link_time+0.4-(0.5*ifi));%blank comes 400ms after linking display

                    % NOTE!!!!!: ASK ABOUT WHETHER WE NEED TO WAIT UNTIL THE MAX RESPONSE WINDOW TO PROCEED. IF SO: GET RID OF THIS IF STATEMENT. IF NOT, KEEP AND GET RID OF IF IN WHILE LOOP BELOW AS ITS REDUNDANT..

                    % wait until the response time cutoff is reached if there has not been a response yet..
                    while ((GetSecs-strt) < resp_cut) && keyIsDown == 0
                        [keyIsDown, firstPress, ~, ~, ~] = KbQueueCheck();
                    end
                    if (keyIsDown == 1) && (swtch2 == 1)
                        t_resp = max(firstPress)-t_disp;
                        swtch2 = 0; % turn off switch bc a response occured
                    end

                    %stop listening for responses
                    KbQueueStop();
                    tim_mat(trials,2,2) = GetSecs-strt; % save total trial time..
                    tim_mat2(trials,2,2) = GetSecs;
                    % Record Response
                    if swtch2 == 0
                        tim_mat(trials,1,2) = max(firstPress)-t_disp; % record time of response
                        resp_mat(trials,2) = KbName(firstPress); % record button press
                        respns = KbName(firstPress);
                    else
                        tim_mat(trials,1,2) = -999; % record time of response
                        resp_mat(trials,2) = "!"; % record time of response
                    end

            end

            % If feedback is on, provide the appropriate feedback
            switch fback_On
                case 1
                    if respns == tcon   % vert key pressed
                        response(trials,2) = 1;
                        [backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack,fix_im,window);
                        answer = feedback(1);
                        %answer = strcat(answer," ","button pressed:", tdesc_file{trials,9});
                        tfback0 = GetSecs;
                        while (GetSecs-tfback0)<=feedbackTime
                            Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed
                            Screen('TextSize',window, fback_txtSize);
                            DrawFormattedText(window, answer,'center', scrn_top+30*scrn_1percY, black);
                            Screen('Flip',window);
                        end
                    else
                        response(trials,2) = 0;
                        [backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack,fix_im,window);
                        answer = feedback(2);
                        %answer = strcat(answer,' ','button pressed:', tdesc_file{trials,9});
                        tfback0 = GetSecs;
                        while (GetSecs-tfback0)<=feedbackTime
                            Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed
                            Screen('TextSize',window, fback_txtSize);
                            DrawFormattedText(window, answer,'center', scrn_top+30*scrn_1percY, black);
                            Screen('Flip',window);
                        end
                    end
                case 0
                    if respns == tcon
                        response(trials,2) = 1;
                    else
                        response(trials,2) = 0;
                    end
            end

            % Inter-trial Interval..
            [backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack,fix_im,window);
            iti_0 = GetSecs;
            while (GetSecs-iti_0)<ITI
                Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed
                Screen('Flip',window);
            end

            tim_mat2(trials,4,2) = GetSecs;
            tim_mat2(trials,3,2) = iti_0;

            tim_mat(trials,3,2) = tim_mat2(trials,4,2)-strt; % save total trial time after ITI too..

            % Flush the Response Queue..
            KbQueueFlush();

            % Add the new datum
            q2 = QuestUpdate(q2,log10(param2(trials)),response(trials,2));

        end
        %% Staircase 3
        iQuestNum=3;
        if n_int>2
            if link=='L'
                param3(trials)=round(10^(QuestQuantile(q3)));   % get next level of parameter to test
            elseif link=='B'
                param3(trials)=quant(10^(QuestQuantile(q3)),0.1);  % get next level of parameter to test
            elseif link=='F'
                param3(trials)=quant(10^(QuestQuantile(q3)),0.1);  % get next level of parameter to test
            end
            if isnan(param3(trials))
                param3(trials)= param3(trials-1);
            end

            %EJD COMMENTED/EDITED BELOW
            %param1(trials) = max(minparam,min(maxparam,param1(trials)));% properly restrict value to useable range

            if maxmin_seg == 1
                % limit to range btw maxparam and minparam
                if param3(trials) > seg3max
                    param3(trials) = seg3max;
                end
                if param3(trials) < seg3min
                    param3(trials) = seg3min;
                end
            end

            if maxmin_seg == 0
                % limit to range btw maxparam and minparam
                if param3(trials) > maxparam
                    param3(trials) = maxparam;
                end
                if param3(trials) < minparam
                    param3(trials) = minparam;
                end
            end

            switch link
                case 'F'
                    tcon = trialprop{trials,3}; % EJD NOTE: was "orient" now.. "face"...
                    corr_resp_mat(trials,3) = tcon;
                    [backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack,fix_im,window);
                    Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed
                    Screen('Flip',window);
                    strt = GetSecs;
                    tim_mat2(trials,1,3) = strt;

                if tcon==trial_con(2)
                    % grab the image from  img_array with parameter closest to one
                    % suggested by from quest (param1)
                    a=par_inxsC2(:,end)';
                    n=param3(end); % assign the most recent param in param1 to "n"
                    [val,idx]=min(abs(a-n));
                    closest_param=a(idx);
                    % Check if there are multiple images which give this
                    % value.. if so pick one at random..
                    [row, col] = find(a == closest_param);
                    idxTmp = randi([1, length(col)]);
                    idx=col(idxTmp);

                    param3(trials)=closest_param;  % get next level of parameter to test
                    con_mat(trials,3) = "Quest3";
                    im_mat(trials,3) = img_arrayC2{idx,2}; % Save image name/path
                    ImParMat(trials,3) = closest_param;
                    %Get the object parameters at index idx too..
                    ObjParMat(trials,1:end,3)=par_inxsC2(idx,1:length(parTagz));
                    %occCon_mat(trials,3)=img_arrayC2{idx,6}; % save the occlusion condition
                    %!!??
                    image_in = img_arrayC2{idx,1};
                end

               
                if tcon==trial_con(1)
                    % grab the image from  img_array with parameter closest to one
                    % suggested by from quest (param1)
                    a=par_inxsC1(:,end)';
                    n=param3(end); % assign the most recent param in param1 to "n"
                    [val,idx]=min(abs(a-n));
                    closest_param=a(idx);
                    % Check if there are multiple images which give this
                    % value.. if so pick one at random..
                    [row, col] = find(a == closest_param);
                    idxTmp = randi([1, length(col)]);
                    idx=col(idxTmp);

                    param3(trials)=closest_param;  % get next level of parameter to test
                    con_mat(trials,3) = "Quest3";
                    im_mat(trials,3) = img_arrayC1{idx,2}; % Save image name/path
                    ImParMat(trials,3) = closest_param;
                    %Get the object parameters at index idx too..
                    ObjParMat(trials,1:end,3)=par_inxsC1(idx,1:length(parTagz));
                    %occCon_mat(trials,3)=img_arrayC1{idx,6}; % save the occlusion condition
                    %!!??
                    image_in = img_arrayC1{idx,1};
                end

                if noise==1
%                    % create a random permutation of the noise row indices
%                     NoiseRanIdxOrder= randperm(size(noise_arrayIdx, 1));
%                     % reorder the rows of the noise_arrayIdx matrix using the permutation
%                     noise_arrayIdx = noise_arrayIdx(NoiseRanIdxOrder, :);
%                     noiseImgIn = noise_array{noise_arrayIdx(1,1),1}; %get the noise image for the index on top..

                    noiseImgIn = noise_array{nzFrmCountr,1};
                    
                    [image_in] = revCorrAddNoise_v1(image_in,noiseImgIn,nWt,fixIndOlayIm,nonFxnIndMask); % add the noise to the image..
                    noizMat(trials,3) = noise_array{nzFrmCountr,2};
                    %noise_arrayIdx = noise_arrayIdx(2:end,1); % now eliminate that index from the list so the same noise image won't be selected again..
                    nzFrmCountr=nzFrmCountr+1; % update countr..
                end

                    % draw linking display to offscreen window
                    [backwindow] = PORClinkRevCor_v1(mainsrc, window, backwindow,trialprop(trials,1), trial_con, image_in,fixrect,grayblack);
                    %Screen('DrawTexture', window, backwindow, mainsrc, maindest, 0, 0, 1);  % draw stimulus to the on-screen window, rotated & zoomed
                    Screen('DrawTexture', window, backwindow, mainsrc, maindest, 0, 0, 1);  % draw stimulus to the on-screen window, rotated & zoomed
                    %[backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack,fix_im);
                    %Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed

                    keyIsDown = 0;
                    swtch = 0;
                    swtch2 = 0;
                    respns = "";
                    %t_trial0 = GetSecs;
                    % While stimulus is on the screen and no response
                    while ((GetSecs-strt) <= presTime2) && keyIsDown == 0 %trialResp == 0

                        % Draw stimulus
                        if ((GetSecs-strt)>=presTime) && swtch == 0
                            Screen('DrawTexture', window, backwindow, mainsrc, maindest, 0, 0, 1);  % draw stimulus to the on-screen window, rotated & zoomed
                            Screen('Flip', window);% show linking display
                            t_disp = GetSecs; % record time display came up..
                            %start listening for responses
                            KbQueueStart();
                            swtch = 1;
                        end

                        [keyIsDown, firstPress, ~, ~, ~] = KbQueueCheck();
                    end
                    % Check to see if we exited the While loop above due to a
                    % response.. if not, raise swtch2 to indicate we still don't have a response..
                    if keyIsDown == 0
                        swtch2 = 1;
                    end

                    % remove stimulus image/swap for just fixation..
                    [backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack, fix_im,window);
                    Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed
                    Screen('Flip',window);
                    %[blank_time StimulusOnsetTime FlipTimestamp Missed Beampos] = Screen('Flip', window, link_time+0.4-(0.5*ifi));%blank comes 400ms after linking display

                    % NOTE!!!!!: ASK ABOUT WHETHER WE NEED TO WAIT UNTIL THE MAX RESPONSE WINDOW TO PROCEED. IF SO: GET RID OF THIS IF STATEMENT. IF NOT, KEEP AND GET RID OF IF IN WHILE LOOP BELOW AS ITS REDUNDANT..

                    % wait until the response time cutoff is reached if there has not been a response yet..
                    while ((GetSecs-strt) < resp_cut) && keyIsDown == 0
                        [keyIsDown, firstPress, ~, ~, ~] = KbQueueCheck();
                    end
                    if (keyIsDown == 1) && (swtch2 == 1)
                        t_resp = max(firstPress)-t_disp;
                        swtch2 = 0; % turn off switch bc a response occured
                    end

                    %stop listening for responses
                    KbQueueStop();
                    tim_mat(trials,2,3) = GetSecs-strt; % save total trial time..
                    tim_mat2(trials,2,3) = GetSecs;
                    % Record Response
                    if swtch2 == 0
                        tim_mat(trials,1,3) = max(firstPress)-t_disp; % record time of response
                        resp_mat(trials,3) = KbName(firstPress); % record button press
                        respns = KbName(firstPress);
                    else
                        tim_mat(trials,1,3) = -999; % record time of response
                        resp_mat(trials,3) = "!"; % record time of response
                    end

            end

            % If feedback is on, provide the appropriate feedback
            switch fback_On
                case 1
                    if respns == tcon   % vert key pressed
                        response(trials,3) = 1;
                        [backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack,fix_im,window);
                        answer = feedback(1);
                        %answer = strcat(answer," ","button pressed:", tdesc_file{trials,9});
                        tfback0 = GetSecs;
                        while (GetSecs-tfback0)<=feedbackTime
                            Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed
                            Screen('TextSize',window, fback_txtSize);
                            DrawFormattedText(window, answer,'center', scrn_top+30*scrn_1percY, black);
                            Screen('Flip',window);
                        end
                    else
                        response(trials,3) = 0;
                        [backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack,fix_im,window);
                        answer = feedback(2);
                        %answer = strcat(answer,' ','button pressed:', tdesc_file{trials,9});
                        tfback0 = GetSecs;
                        while (GetSecs-tfback0)<=feedbackTime
                            Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed
                            Screen('TextSize',window, fback_txtSize);
                            DrawFormattedText(window, answer,'center', scrn_top+30*scrn_1percY, black);
                            Screen('Flip',window);
                        end
                    end
                case 0
                    if respns == tcon
                        response(trials,3) = 1;
                    else
                        response(trials,3) = 0;
                    end
            end

            % Inter-trial Interval..
            [backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack,fix_im,window);
            iti_0 = GetSecs;
            while (GetSecs-iti_0)<ITI
                Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed
                Screen('Flip',window);
            end

            tim_mat2(trials,4,3) = GetSecs;
            tim_mat2(trials,3,3) = iti_0;

            tim_mat(trials,3,3) = tim_mat2(trials,4,3)-strt; % save total trial time after ITI too..

            % Flush the Response Queue..
            KbQueueFlush();

            % Add the new datum
            q3 = QuestUpdate(q3,log10(param3(trials)),response(trials,3));
        end
        %% Staircase 4
        if n_int>3
            if link=='L'
                param4(trials)=round(10^(QuestQuantile(q4)));   % get next level of parameter to test
            elseif link=='B'
                param4(trials)=quant(10^(QuestQuantile(q4)),0.1);  % get next level of parameter to test
            elseif link=='F'
                param4(trials)=quant(10^(QuestQuantile(q4)),0.1);  % get next level of parameter to test
            end
            if isnan(param4(trials))
                param4(trials)= param4(trials-1);
            end

            %EJD COMMENTED/EDITED BELOW
            %param1(trials) = max(minparam,min(maxparam,param1(trials)));% properly restrict value to useable range

            if maxmin_seg == 1
                % limit to range btw maxparam and minparam
                if param4(trials) > seg4max
                    param4(trials) = seg4max;
                end
                if param4(trials) < seg4min
                    param4(trials) = seg4min;
                end
            end


            if maxmin_seg == 0
                % limit to range btw maxparam and minparam
                if param4(trials) > maxparam
                    param4(trials) = maxparam;
                end
                if param4(trials) < minparam
                    param4(trials) = minparam;
                end
            end

            switch link
                case 'F'
                    tcon = trialprop{trials,4}; % EJD NOTE: was "orient" now.. "face"...
                    corr_resp_mat(trials,4) = tcon;
                    [backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack,fix_im,window);
                    Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed
                    Screen('Flip',window);
                    strt = GetSecs;
                    tim_mat2(trials,1,4) = strt;


                 if tcon==trial_con(2)
                    % grab the image from  img_array with parameter closest to one
                    % suggested by from quest (param1)
                    a=par_inxsC2(:,end)';
                    n=param4(end); % assign the most recent param in param1 to "n"
                    [val,idx]=min(abs(a-n));
                    closest_param=a(idx);
                    % Check if there are multiple images which give this
                    % value.. if so pick one at random..
                    [row, col] = find(a == closest_param);
                    idxTmp = randi([1, length(col)]);
                    idx=col(idxTmp);

                    param4(trials)=closest_param;  % get next level of parameter to test
                    con_mat(trials,4) = "Quest4";
                    im_mat(trials,4) = img_arrayC2{idx,2}; % Save image name/path
                    ImParMat(trials,4) = closest_param;
                    %Get the object parameters at index idx too..
                    ObjParMat(trials,1:end,4)=par_inxsC2(idx,1:length(parTagz));
                    %occCon_mat(trials,4)=img_arrayC2{idx,6}; % save the occlusion condition
                    %!!??
                    image_in = img_arrayC2{idx,1};
                end

               
                if tcon==trial_con(1)
                    % grab the image from  img_array with parameter closest to one
                    % suggested by from quest (param1)
                    a=par_inxsC1(:,end)';
                    n=param4(end); % assign the most recent param in param1 to "n"
                    [val,idx]=min(abs(a-n));
                    closest_param=a(idx);
                    % Check if there are multiple images which give this
                    % value.. if so pick one at random..
                    [row, col] = find(a == closest_param);
                    idxTmp = randi([1, length(col)]);
                    idx=col(idxTmp);

                    param4(trials)=closest_param;  % get next level of parameter to test
                    con_mat(trials,4) = "Quest4";
                    im_mat(trials,4) = img_arrayC1{idx,2}; % Save image name/path
                    ImParMat(trials,4) = closest_param;
                    %Get the object parameters at index idx too..
                    ObjParMat(trials,1:end,4)=par_inxsC1(idx,1:length(parTagz));
                    %occCon_mat(trials,4)=img_arrayC1{idx,6}; % save the occlusion condition
                    %!!??
                    image_in = img_arrayC1{idx,1};
                end

                if noise==1
%                     % create a random permutation of the noise row indices
%                     NoiseRanIdxOrder= randperm(size(noise_arrayIdx, 1));
%                     % reorder the rows of the noise_arrayIdx matrix using the permutation
%                     noise_arrayIdx = noise_arrayIdx(NoiseRanIdxOrder, :);
%                     noiseImgIn = noise_array{noise_arrayIdx(1,1),1}; %get the noise image for the index on top..

                    noiseImgIn = noise_array{nzFrmCountr,1};
                    
                    [image_in] = revCorrAddNoise_v1(image_in,noiseImgIn,nWt,fixIndOlayIm,nonFxnIndMask); % add the noise to the image..
                    noizMat(trials,4) = noise_array{nzFrmCountr,2};
                    %noise_arrayIdx = noise_arrayIdx(2:end,1); % now eliminate that index from the list so the same noise image won't be selected again..
                    nzFrmCountr=nzFrmCountr+1; % update countr..
                end

                    % draw linking display to offscreen window
                    [backwindow] = PORClinkRevCor_v1(mainsrc, window, backwindow,trialprop(trials,1), trial_con, image_in,fixrect,grayblack);
                    %Screen('DrawTexture', window, backwindow, mainsrc, maindest, 0, 0, 1);  % draw stimulus to the on-screen window, rotated & zoomed
                    Screen('DrawTexture', window, backwindow, mainsrc, maindest, 0, 0, 1);  % draw stimulus to the on-screen window, rotated & zoomed
                    %[backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack,fix_im);
                    %Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed

                    keyIsDown = 0;
                    swtch = 0;
                    swtch2 = 0;
                    respns = "";
                    %t_trial0 = GetSecs;
                    % While stimulus is on the screen and no response
                    while ((GetSecs-strt) <= presTime2) && keyIsDown == 0 %trialResp == 0

                        % Draw stimulus
                        if ((GetSecs-strt)>=presTime) && swtch == 0
                            Screen('DrawTexture', window, backwindow, mainsrc, maindest, 0, 0, 1);  % draw stimulus to the on-screen window, rotated & zoomed
                            Screen('Flip', window);% show linking display
                            t_disp = GetSecs; % record time display came up..
                            %start listening for responses
                            KbQueueStart();
                            swtch = 1;
                        end

                        [keyIsDown, firstPress, ~, ~, ~] = KbQueueCheck();
                    end
                    % Check to see if we exited the While loop above due to a
                    % response.. if not, raise swtch2 to indicate we still don't have a response..
                    if keyIsDown == 0
                        swtch2 = 1;
                    end

                    % remove stimulus image/swap for just fixation..
                    [backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack, fix_im,window);
                    Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed
                    Screen('Flip',window);
                    %[blank_time StimulusOnsetTime FlipTimestamp Missed Beampos] = Screen('Flip', window, link_time+0.4-(0.5*ifi));%blank comes 400ms after linking display


                    % NOTE!!!!!: ASK ABOUT WHETHER WE NEED TO WAIT UNTIL THE MAX RESPONSE WINDOW TO PROCEED. IF SO: GET RID OF THIS IF STATEMENT. IF NOT, KEEP AND GET RID OF IF IN WHILE LOOP BELOW AS ITS REDUNDANT..

                    % wait until the response time cutoff is reached if there has not been a response yet..
                    while ((GetSecs-strt) < resp_cut) && keyIsDown == 0
                        [keyIsDown, firstPress, ~, ~, ~] = KbQueueCheck();
                    end
                    if (keyIsDown == 1) && (swtch2 == 1)
                        t_resp = max(firstPress)-t_disp;
                        swtch2 = 0; % turn off switch bc a response occured
                    end

                    %stop listening for responses
                    KbQueueStop();
                    tim_mat(trials,2,4) = GetSecs-strt; % save total trial time..
                    tim_mat2(trials,2,4) = GetSecs;
                    % Record Response
                    if swtch2 == 0
                        tim_mat(trials,1,4) = max(firstPress)-t_disp; % record time of response
                        resp_mat(trials,4) = KbName(firstPress); % record button press
                        respns = KbName(firstPress);
                    else
                        tim_mat(trials,1,4) = -999; % record time of response
                        resp_mat(trials,4) = "!"; % record time of response
                    end

            end

            % If feedback is on, provide the appropriate feedback
            switch fback_On
                case 1
                    if respns == tcon   % vert key pressed
                        response(trials,4) = 1;
                        [backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack,fix_im,window);
                        answer = feedback(1);
                        %answer = strcat(answer," ","button pressed:", tdesc_file{trials,9});
                        tfback0 = GetSecs;
                        while (GetSecs-tfback0)<=feedbackTime
                            Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed
                            Screen('TextSize',window, fback_txtSize);
                            DrawFormattedText(window, answer,'center', scrn_top+30*scrn_1percY, black);
                            Screen('Flip',window);
                        end
                    else
                        response(trials,4) = 0;
                        [backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack,fix_im,window);
                        answer = feedback(2);
                        %answer = strcat(answer,' ','button pressed:', tdesc_file{trials,9});
                        tfback0 = GetSecs;
                        while (GetSecs-tfback0)<=feedbackTime
                            Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed
                            Screen('TextSize',window, fback_txtSize);
                            DrawFormattedText(window, answer,'center', scrn_top+30*scrn_1percY, black);
                            Screen('Flip',window);
                        end
                    end
                case 0
                    if respns == tcon
                        response(trials,4) = 1;
                    else
                        response(trials,4) = 0;
                    end
            end

            % Inter-trial Interval..
            [backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack,fix_im,window);
            iti_0 = GetSecs;
            while (GetSecs-iti_0)<ITI
                Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed
                Screen('Flip',window);
            end

            tim_mat2(trials,4,4) = GetSecs;
            tim_mat2(trials,3,4) = iti_0;

            tim_mat(trials,3,4) = tim_mat2(trials,4,4)-strt; % save total trial time after ITI too..

            % Flush the Response Queue..
            KbQueueFlush();

            % Add the new datum
            q4 = QuestUpdate(q4,log10(param4(trials)),response(trials,4));
        end
        %% Staircase 5
        if n_int>4
            if link=='L'
                param5(trials)=round(10^(QuestQuantile(q5)));   % get next level of parameter to test
            elseif link=='B'
                param5(trials)=quant(10^(QuestQuantile(q5)),0.1);  % get next level of parameter to test
            elseif link=='F'
                param5(trials)=quant(10^(QuestQuantile(q5)),0.1);  % get next level of parameter to test
            end
            if isnan(param5(trials))
                param5(trials)= param5(trials-1);
            end

            %EJD COMMENTED/EDITED BELOW
            %param1(trials) = max(minparam,min(maxparam,param1(trials)));% properly restrict value to useable range

            if maxmin_seg == 1
                % limit to range btw maxparam and minparam
                if param5(trials) > seg5max
                    param5(trials) = seg5max;
                end
                if param5(trials) < seg5min
                    param5(trials) = seg5min;
                end
            end

            if maxmin_seg == 0
                % limit to range btw maxparam and minparam
                if param5(trials) > maxparam
                    param5(trials) = maxparam;
                end
                if param5(trials) < minparam
                    param5(trials) = minparam;
                end
            end

            switch link
                case 'F'
                    tcon = trialprop{trials,5}; % EJD NOTE: was "orient" now.. "face"...
                    corr_resp_mat(trials,5) = tcon;
                    [backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack,fix_im,window);
                    Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed
                    Screen('Flip',window);
                    strt = GetSecs;
                    tim_mat2(trials,1,5) = strt;

                if tcon==trial_con(2)
                    % grab the image from  img_array with parameter closest to one
                    % suggested by from quest (param1)
                    a=par_inxsC2(:,end)';
                    n=param5(end); % assign the most recent param in param1 to "n"
                    [val,idx]=min(abs(a-n));
                    closest_param=a(idx);
                    % Check if there are multiple images which give this
                    % value.. if so pick one at random..
                    [row, col] = find(a == closest_param);
                    idxTmp = randi([1, length(col)]);
                    idx=col(idxTmp);

                    param5(trials)=closest_param;  % get next level of parameter to test
                    con_mat(trials,5) = "Quest5";
                    im_mat(trials,5) = img_arrayC2{idx,2}; % Save image name/path
                    ImParMat(trials,5) = closest_param;
                    %Get the object parameters at index idx too..
                    ObjParMat(trials,1:end,5)=par_inxsC2(idx,1:length(parTagz));
                    %occCon_mat(trials,5)=img_arrayC2{idx,6}; % save the occlusion condition
                    %!!??
                    image_in = img_arrayC2{idx,1};
                end

               
                if tcon==trial_con(1)
                    % grab the image from  img_array with parameter closest to one
                    % suggested by from quest (param1)
                    a=par_inxsC1(:,end)';
                    n=param5(end); % assign the most recent param in param1 to "n"
                    [val,idx]=min(abs(a-n));
                    closest_param=a(idx);
                    % Check if there are multiple images which give this
                    % value.. if so pick one at random..
                    [row, col] = find(a == closest_param);
                    idxTmp = randi([1, length(col)]);
                    idx=col(idxTmp);

                    param5(trials)=closest_param;  % get next level of parameter to test
                    con_mat(trials,5) = "Quest5";
                    im_mat(trials,5) = img_arrayC1{idx,2}; % Save image name/path
                    ImParMat(trials,5) = closest_param;
                    %Get the object parameters at index idx too..
                    ObjParMat(trials,1:end,5)=par_inxsC1(idx,1:length(parTagz));
                    %occCon_mat(trials,5)=img_arrayC1{idx,6}; % save the occlusion condition
                    %!!??
                    image_in = img_arrayC1{idx,1};
                end

                if noise==1
%                     % create a random permutation of the noise row indices
%                     NoiseRanIdxOrder= randperm(size(noise_arrayIdx, 1));
%                     % reorder the rows of the noise_arrayIdx matrix using the permutation
%                     noise_arrayIdx = noise_arrayIdx(NoiseRanIdxOrder, :);
%                     noiseImgIn = noise_array{noise_arrayIdx(1,1),1}; %get the noise image for the index on top..

                    noiseImgIn = noise_array{nzFrmCountr,1};
                    
                    [image_in] = revCorrAddNoise_v1(image_in,noiseImgIn,nWt,fixIndOlayIm,nonFxnIndMask); % add the noise to the image..
                    noizMat(trials,5) = noise_array{nzFrmCountr,2};
                    %noise_arrayIdx = noise_arrayIdx(2:end,1); % now eliminate that index from the list so the same noise image won't be selected again..
                    nzFrmCountr=nzFrmCountr+1; % update countr..
                end

                    % draw linking display to offscreen window
                    [backwindow] = PORClinkRevCor_v1(mainsrc, window, backwindow,trialprop(trials,1), trial_con, image_in,fixrect,grayblack);
                    %Screen('DrawTexture', window, backwindow, mainsrc, maindest, 0, 0, 1);  % draw stimulus to the on-screen window, rotated & zoomed
                    Screen('DrawTexture', window, backwindow, mainsrc, maindest, 0, 0, 1);  % draw stimulus to the on-screen window, rotated & zoomed
                    %[backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack,fix_im);
                    %Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed

                    keyIsDown = 0;
                    swtch = 0;
                    swtch2 = 0;
                    respns = "";
                    %t_trial0 = GetSecs;
                    % While stimulus is on the screen and no response
                    while ((GetSecs-strt) <= presTime2) && keyIsDown == 0 %trialResp == 0

                        % Draw stimulus
                        if ((GetSecs-strt)>=presTime) && swtch == 0
                            Screen('DrawTexture', window, backwindow, mainsrc, maindest, 0, 0, 1);  % draw stimulus to the on-screen window, rotated & zoomed
                            Screen('Flip', window);% show linking display
                            t_disp = GetSecs; % record time display came up..
                            %start listening for responses
                            KbQueueStart();
                            swtch = 1;
                        end

                        [keyIsDown, firstPress, ~, ~, ~] = KbQueueCheck();
                    end
                    % Check to see if we exited the While loop above due to a
                    % response.. if not, raise swtch2 to indicate we still don't have a response..
                    if keyIsDown == 0
                        swtch2 = 1;
                    end

                    % remove stimulus image/swap for just fixation..
                    [backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack, fix_im,window);
                    Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed
                    Screen('Flip',window);
                    %[blank_time StimulusOnsetTime FlipTimestamp Missed Beampos] = Screen('Flip', window, link_time+0.4-(0.5*ifi));%blank comes 400ms after linking display

                    % NOTE!!!!!: ASK ABOUT WHETHER WE NEED TO WAIT UNTIL THE MAX RESPONSE WINDOW TO PROCEED. IF SO: GET RID OF THIS IF STATEMENT. IF NOT, KEEP AND GET RID OF IF IN WHILE LOOP BELOW AS ITS REDUNDANT..

                    % wait until the response time cutoff is reached if there has not been a response yet..
                    while ((GetSecs-strt) < resp_cut) && keyIsDown == 0
                        [keyIsDown, firstPress, ~, ~, ~] = KbQueueCheck();
                    end
                    if (keyIsDown == 1) && (swtch2 == 1)
                        t_resp = max(firstPress)-t_disp;
                        swtch2 = 0; % turn off switch bc a response occured
                    end

                    %stop listening for responses
                    KbQueueStop();
                    tim_mat(trials,2,5) = GetSecs-strt; % save total trial time..
                    tim_mat2(trials,2,5) = GetSecs;
                    % Record Response
                    if swtch2 == 0
                        tim_mat(trials,1,5) = max(firstPress)-t_disp; % record time of response
                        resp_mat(trials,5) = KbName(firstPress); % record button press
                        respns = KbName(firstPress);
                    else
                        tim_mat(trials,1,5) = -999; % record time of response
                        resp_mat(trials,5) = "!"; % record time of response
                    end

            end

            % If feedback is on, provide the appropriate feedback
            switch fback_On
                case 1
                    if respns == tcon   % vert key pressed
                        response(trials,5) = 1;
                        [backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack,fix_im,window);
                        answer = feedback(1);
                        %answer = strcat(answer," ","button pressed:", tdesc_file{trials,9});
                        tfback0 = GetSecs;
                        while (GetSecs-tfback0)<=feedbackTime
                            Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed
                            Screen('TextSize',window, fback_txtSize);
                            DrawFormattedText(window, answer,'center', scrn_top+30*scrn_1percY, black);
                            Screen('Flip',window);
                        end
                    else
                        response(trials,5) = 0;
                        [backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack,fix_im,window);
                        answer = feedback(2);
                        %answer = strcat(answer,' ','button pressed:', tdesc_file{trials,9});
                        tfback0 = GetSecs;
                        while (GetSecs-tfback0)<=feedbackTime
                            Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed
                            Screen('TextSize',window, fback_txtSize);
                            DrawFormattedText(window, answer,'center', scrn_top+30*scrn_1percY, black);
                            Screen('Flip',window);
                        end
                    end
                case 0
                    if respns == tcon
                        response(trials,5) = 1;
                    else
                        response(trials,5) = 0;
                    end
            end

            % Inter-trial Interval..
            [backwindow_fix,backwindow] = PLOTfix(backwindow_fix,backwindow,fixrect,grayblack,fix_im,window);
            iti_0 = GetSecs;
            while (GetSecs-iti_0)<ITI
                Screen('DrawTexture', window, backwindow_fix, mainsrc, maindest, 0, 0);  % draw fixation screen to the on-screen window, rotated & zoomed
                Screen('Flip',window);
            end

            tim_mat2(trials,4,5) = GetSecs;
            tim_mat2(trials,3,5) = iti_0;

            tim_mat(trials,3,5) = tim_mat2(trials,4,5)-strt; % save total trial time after ITI too..

            % Flush the Response Queue..
            KbQueueFlush();

            % Add the new datum
            q5 = QuestUpdate(q5,log10(param5(trials)),response(trials,5));
        end
    end
    % clear screen after last trial
    Screen('FillRect', window, gray);
    Screen('Flip',window);
    %strt1 = GetSecs;
    %[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed Beampos] = Screen('Flip', window, strt1+(1.5*ifi));

    %% Build tdf_out output array
    
    nParz=length(parTagz);
    tdf_out = cell(numtrials*n_int,(19+nParz));
    
    countr = 1;
    for ii = 1:numtrials
        for zz = 1:n_int
            tdf_out(countr,1) = cellstr(im_mat(ii,zz));
            tdf_out(countr,2) = cellstr(con_mat(ii,zz));
            % add the occlusion conditions
            %tdf_out(countr,3) = num2cell(occCon_mat(ii,zz));
            % add the image matrices in too
            tdf_out{countr,5} = imread(tdf_out{countr,1});
            tdf_out(countr,7) = cellstr(corr_resp_mat(ii,zz));
            %add scoring column
            tdf_out(countr,12) = num2cell(response(ii,zz));
            tdf_out(countr,8:10) = num2cell(tim_mat(ii,:,zz));
            tdf_out(countr,11) = cellstr(resp_mat(ii,zz));
            % add the raw time data..
            tdf_out(countr,13:16) = num2cell(tim_mat2(ii,:,zz));
            
            %add the noise image data..
            tdf_out(countr,17) = cellstr(noizMat(ii,zz));
            tdf_out{countr,18} = imread(tdf_out{countr,17});
            
            % Add the parameters columns..
            tdf_out(countr,19) = num2cell(ImParMat(ii,zz)); %Difference Parameter/Questpar
            % Raw Object Parameters
            tdf_out(countr,20:end)=num2cell(ObjParMat(ii,1:end,zz));
            
            countr = countr+1;
        end 
    end
    clear countr
    
    %convert coded values for no response to "NA"
    for ii = 1:size(tdf_out,1)
        if tdf_out{ii,8} == -999
            tdf_out{ii,8} = "NA";
        end

        if tdf_out{ii,11} == '!'
            tdf_out{ii,11}= "NA";
        end
    end

    if qstParStr == "lumDiff" || qstParStr == "texDiff"
                
                headers = {"ImgFile_Path","QUEST_Con","occCon","questPar","Image","Randomization_Col","Correct_RespKey", "T_Resp","T_Trial","T_Trial_wITI","Subj_RespKey","Correctness","tStart_Raw","tEnd_Raw","tITIstart_Raw","tITIend_Raw","noiseImgFile","noiseImg",qstParStr,obj1Par,obj2Par};
    end
    
    if qstParStr == "CRpl"
                headers = {"ImgFile_Path","QUEST_Con","occCon","questPar","Image","Randomization_Col","Correct_RespKey", "T_Resp","T_Trial","T_Trial_wITI","Subj_RespKey","Correctness","tStart_Raw","tEnd_Raw","tITIstart_Raw","tITIend_Raw","noiseImgFile","noiseImg",strcat(qstParStr,"_difFrmBkgd"),obj1Par};
    end
    tdf_out(:,4)={qstParStr}; % add the quest parameter string to each row in the "questPar" column..

    % Add the occlusion condition data/pull from filenames..
    for ww=1:size(tdf_out,1)
        fname=tdf_out{ww,1};
        [occConTmp, ~] = xtractParsFrmFilename(fname, "occ","num");
        tdf_out(ww,3)={occConTmp};
    end

    % Check if any additional parameters are being varied other than the
    % quest-modulated-parameter. If so, get the values out of the filename
    % and save the values in a set of additional columns. Then contatenate
    % them horizontally onto the end.
    %% ----------------------------------------------------------------------
    parOnOffTagz={["CR","CRpl"],["T","texDiff"],["L","lumDiff"]};
    subParz={{"CRpw","CRpl","CRal","CRob","CRdm"},{"Tr1","Tr2","Tal"},{"lum1","lum2","LWT","Lal"}};
    nonQstParz={};
    for ww=1:length(parOnOffTagz)
    parOnOffTagz_pair=parOnOffTagz{1,ww};
    % If this is not the QUEST parameter, check to see if this is being
    % used/varied randomly outside of QUEST.
    if parOnOffTagz_pair(2)~=qstParStr
        fname=tdf_out{1,1}; % use first file name at at top column. **(Assuming all images either use or do not use a parameter set..)
        [parOnOff, ~] = xtractParsFrmFilename(fname, parOnOffTagz_pair(1),"num");
        nonQstParz{1,end+1}=parOnOffTagz_pair;
        nonQstParz{2,size(nonQstParz,2)}=parOnOff;
        nonQstParz{3,size(nonQstParz,2)}=subParz{ww};
    end
    end

    nrwz=size(tdf_out,1);
    xtraHeadrz1={};% for on/of cols
    xtraHeadrz2={};% for pars
    onOffColz={};
    xtraColCount=0;
    for ww=1:size(nonQstParz,2)
        parIsOn=nonQstParz{2,ww}(1);
        xtraHeadrz1=horzcat(xtraHeadrz1,{nonQstParz{1,ww}(1)}); % add the on/of par
        xtraColCount=xtraColCount+1; % add this on off col to to total extra col tally
        onOffColz(1:nrwz,end+1)={nonQstParz{2,ww}};
        % If this parameter is on add its
        if parIsOn==1
            % count extra subpar cols and add to tally too..
            valList=nonQstParz{3,ww};
            ncolzTmp=length(valList);
            xtraColCount=xtraColCount+ncolzTmp;
            % add them to xtraHeadrz2
            xtraHeadrz2=horzcat(xtraHeadrz2,valList);
        end
    end

    if length(xtraHeadrz2)>0

    xtraCols=cell(size(tdf_out,1),length(xtraHeadrz2));
    xtraParz=string(zeros(1,length(xtraHeadrz2)));
    for uu=1:size(xtraParz,2)
        xtraParz(1,uu)=xtraHeadrz2{uu};
    end
    
    for vv=1:size(xtraCols,1)
        fname=tdf_out{vv,1}; % get the next file name
        [parzPass, ~] = xtractParsFrmFilename(fname, xtraParz,"num");
        for ss=1:size(parzPass,2)
            xtraCols(vv,ss)={parzPass(ss)};
        end
    
    end

    xtraHeadrz=horzcat(xtraHeadrz1,xtraHeadrz2);
    xtraCols=horzcat(onOffColz,xtraCols);
    else
    xtraHeadrz=xtraHeadrz1;
    xtraCols=onOffColz;
    end
    headers=horzcat(headers,xtraHeadrz);
    tdf_out=horzcat(tdf_out,xtraCols);
    pause="";

    %[parValz, parTagz] = xtractParsFrmFilename(fname,parOnOffTagz_pair(1),"num");
    %ori_con = xtractParsFrmFilename(fname,["ori"],"str");%orientation condition
    %----------------------------------------------------------------------
    name_cells = {headers, tdf_out};
    tdf_out = cat(1, name_cells{:});

    %% Now compute thresholds
    if link=='L'
        t1=round(10^(QuestMean(q1)));		% Recommended by Pelli (1989) and King-Smith et al.(1994).
        sd1=round(10^(QuestSd(q1)));
    elseif link=='B'
        t1=quant(10^(QuestMean(q1)),0.1);		% Recommended by Pelli (1989) and King-Smith et al.(1994).
        sd1=quant(10^(QuestSd(q1)),0.1);
    elseif link=='F' %??? DO WE NEED TO TWEAK THIS? WHAT PARTS OF THE SCRIPT USE "t1" and "sd1"
        %t1=quant(10^(QuestMean(q1)),0.1);		% Recommended by Pelli (1989) and King-Smith et al.(1994).
        %sd1=quant(10^(QuestSd(q1)),0.1);

        if n_int>=1
        t1=quant(10^(QuestMean(q1)),0.0001);		% Recommended by Pelli (1989) and King-Smith et al.(1994).
        sd1=quant(10^(QuestSd(q1)),0.0001);
        else
            t1=[];
            sd1=[];
            q1=[];
        end

        if n_int>=2
        t2=quant(10^(QuestMean(q2)),0.0001);		% Recommended by Pelli (1989) and King-Smith et al.(1994).
        sd2=quant(10^(QuestSd(q2)),0.0001);
        else
            t2=[];
            sd2=[];
            q2=[];
        end

        if n_int>=3
        t3=quant(10^(QuestMean(q3)),0.0001);		% Recommended by Pelli (1989) and King-Smith et al.(1994).
        sd3=quant(10^(QuestSd(q3)),0.0001);
        else
            t3=[];
            sd3=[];
            q3=[];
        end

        if n_int>=4
        t4=quant(10^(QuestMean(q4)),0.0001);		% Recommended by Pelli (1989) and King-Smith et al.(1994).
        sd4=quant(10^(QuestSd(q4)),0.0001);
        else
            t4=[];
            sd4=[];
            q4=[];
        end

        if n_int>=5
        t5=quant(10^(QuestMean(q5)),0.0001);		% Recommended by Pelli (1989) and King-Smith et al.(1994).
        sd5=quant(10^(QuestSd(q5)),0.0001);
        else
            t5=[];
            sd5=[];
            q5=[];
        end
    end

    %% Finish trials and print threshold estimate to screen

    WaitSecs(0.1);
    % WILL NEED TO UPDATE THIS FOR FACE...
    % speechend = sprintf('Your Day 1 session is now complete. Thanks for participating. Please hit ENTER.');
    %if db_mode == 0
    %    speechend = sprintf(' This block has been completed. Please hit ENTER to continue.');
    %elseif db_mode == 1 && dbmode_prtThrsh == 1
        t1_text = num2str(t1);
        pThreshold1_text = num2str(pThreshold1);

        % set unused "param#"s to [] based on # of interleaved quests selected..
        %thrCondz=1; % initialize threshold condition counter
        if n_int<5
            param5=[];
            speechend6 = "No threshold 5 selected..";
        else
            t5_text = num2str(t5);
            pThreshold5_text = num2str(pThreshold5);
            %thrCondz=5;
            speechend6 = strcat("The QUEST estimate of the ",pThreshold5_text," threshold parameter value is: ", t5_text);
        end

        if n_int<4
            param4=[];
            speechend5 = "No threshold 4 selected..";
            
        else
            t4_text = num2str(t4);
            pThreshold4_text = num2str(pThreshold4);
            %thrCondz=4;
            speechend5 = strcat("The QUEST estimate of the ",pThreshold4_text," threshold parameter value is: ", t4_text);

        end

        if n_int<3
            param3=[];
            speechend4 = "No threshold 3 selected..";
            
        else
            t3_text = num2str(t3);
            pThreshold3_text = num2str(pThreshold3);
            %thrCondz=3;
            speechend4 = strcat("The QUEST estimate of the ",pThreshold3_text," threshold parameter value is: ", t3_text);
        end

        if n_int<2
            param2=[];
            speechend3 = "No threshold 2 selected..";
        else
            t2_text = num2str(t2);
            pThreshold2_text = num2str(pThreshold2);
            %thrCondz=2;
            speechend3 = strcat("The QUEST estimate of the ",pThreshold2_text," threshold parameter value is: ", t2_text);
        end

        Quest_conv_value_in_nr = nr_sigmoid(t1, nr_rmax, nr_c50, nr_b, nr_n);


        Quest_conv_value_in_nr = Quest_conv_value_in_nr(2,1);
        Quest_conv_value_in_nr_txt = num2str(Quest_conv_value_in_nr);
        expected_param1_value_txt = num2str(expected_param1_value);

        perc_correct = mean(db_mode_correct_trials);
        perc_correct_txt = num2str(perc_correct);

        speechend = sprintf('This block has been completed. Please hit ENTER to continue.');
        speechend1 = strcat('This block has been completed.');
        speechend2 = strcat("The QUEST estimate of the ",pThreshold1_text," threshold parameter value is: ", t1_text);

        %speechend2 = strcat('The value QUEST converged on was: '," " , t1_text, ' ... Plugged into the NR_function this gives % accuracy of:', " ", Quest_conv_value_in_nr_txt, " ", 'It was supposed to converge on:'," ", expected_param1_value_txt,"...", " ", "to give % accuracy of:", " ", pThreshold1_text, ' ... Please hit ENTER to continue.');
        %speechend3 = strcat("The emperical average percent correct across trails in this run (# correct/# trials) was:", " ", perc_correct_txt);
        
        starz = '*************************************************************';

        disp(starz)
        disp(starz)
        disp(speechend1)
        disp(speechend2)
        disp(speechend3)
        disp(speechend4)
        disp(speechend5)
        disp(speechend6)
        disp(starz)
        disp(starz)

    %end

    %if db_mode == 0 || (db_mode == 1 && dbmode_prtThrsh == 1)
        Screen('TextSize',window, intro_txtSize);
        %Ask(window,speechend,red,gray, ['KbWait(' num2str(kboard) ')'],RectLeft,RectTop,15);
        
        %Ask(window,speechend,red,gray, ['KbWait(' num2str(kboard)')'],RectLeft,RectTop,intro_txtSize); % ORIG
        replyJnk=Ask(window,speechend,red, gray,'GetString',RectLeft,RectTop,intro_txtSize);
        
        % Screen('FillRect', window, [grayblack grayblack grayblack 0]);
        Screen('FillRect', window, BlackIndex(window));
        %-------------------------------------------------------------------------
    %end

catch
    PTB_bp_regain_control()

    %% close stuff up

    ListenChar(0);
    ShowCursor;
    Screen('CloseAll');

end

%% wrap up the output params.
ts_out = [t1,t2,t3,t4,t5];
qs_out = [q1,q2,q3,q4,q5];
sd_out = [sd1,sd2,sd3,sd4,sd5];
% set unused "param#"s to [] based on # of interleaved quests selected..
if n_int<5
    param5=[];
    if n_int<4
        param4=[];
        if n_int<3
            param3=[];
            if n_int<2
                param2=[];
            end
        end
    end
end
params_out = {param1, param2, param3, param4, param5};
nr_ideals = [nr_ideal_t1,nr_ideal_t2,nr_ideal_t3,nr_ideal_t4,nr_ideal_t5];

end
