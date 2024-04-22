%function FACEquest_main5_EJD()
%
% PORCmain4a.m
%
% main script for Perceptual Organization Reverse Correlation experiments in PTB-3
% Makes calls to: PORCquest4a.m PORCpractice4a.m PORClinkB4a.m
%
% this version works well w/ recording responses
%
% Structure of Data Array:
%						   Cue	     |Target      |    Linking          |	 Trial	   | Target |   Key   |
%			 |  Config.	|  Loc.	     |	Loc.      |   	Type            |   Validity   |  Type  |Parameter|
% trialprop: ----------------------------------------------------------------------------------------------
%			 |Vertical/	|UL(1)/UR(2)/|UL(1)/UR(2)/|Motion/Luminance/    |Valid/Same/   |conveX/ |  1/2/3  |
%			 |Horizontal|LR(3)/LL(4) |LR(3)/LL(4) |Region/Texture/Bound.|              |        |   4/5   |
%            |          |            |            |coll./colr+bound=D   |Different     |concaVe |         |
%			 ----------------------------------------------------------------------------------------------
%
%
% Written by Adam Greenberg, CMU/Psych.
% Sept, 2010
% asg 11/2/10 modified from PORCmain1a.m to (a)track 70% in staircase, (b) reduce No_Occluder sessions by 1/2, & (c) use 2 base images (instead of 1) for each orientation (horz. & vert.)
% asg 4/19/11 modified from PORCmain1b.m to include boundary collinearity stimuli
% asg 6/9/11  modified from PORCmain2a.m to include combined luminance/collinearity stimuli

% EJD began updating/developing in 02/2022

%% Parameters
%clc
%clear mex
%clear all
clear;

Screen('Preference', 'SkipSyncTests', 1); %EJD: added to skip sync tests
Screen('Preference','SuppressAllWarnings', 1);
%matdir = pwd;

% Detect and set directory/path params:
% -------------------------------------
% get matlab directory path by "which-ing" for this file..
% NOTE: key assumption: there is only one copy of this on the path.. (that
% should always be the case to avoid other conflicts..)
mlab_dir = fileparts(which('RevCorr_QSTmain6.m'));

diaryfile=strcat(mlab_dir,"/commandWindowDiary");
diary(diaryfile);

cd(mlab_dir);
% cd 2 directory to enter the main dir..
cd ..
cd ..
% save the path..
main_dir = pwd;

% go back to the matlab dir..
cd(mlab_dir)
% -------------------------------------

% get the matlab release runnung on this machine
[mlb_v, mlb_d] = version;
expression = '(R[0-9][0-9][0-9][0-9].)';
mlb_version = regexp(mlb_v,expression,'match');
mlb_year=str2double(mlb_version{1,1}(2:end-1));

% Serial Testing/Development Params..
serial_test = 0;

%quest_testreps = 3;
quest_testreps_init = 1; % number of reps for the initial/first quest block
quest_testreps_reg = 3; % number of reps for the rest of the quest blocks (non-first-blocks)


test_initials = "XXX"; %phoney initials for testing...
print_figs = 1;
print_fig_gif = 1;
gif_delay = 0.1;
mean_sd_fitCrv_fig = 1;
%figs_dir = "/Users/eduwell/OneDrive - mcw.edu/Documents/SNAP/projects/Priyanka_Faces/quest_dev/quest_pics/test2_bin10_thrs70-75-84/";
if serial_test == 0
    %quest_testreps = 1;
    %print_figs = 0;
end

red = [255, 0, 0, 255];	% red
link = 'F'; % 'B'=boundary collinearity; 'L'=luminance; 'Z'=luminance/collinearity combined; 'T' = Texture (EJD ADDED); 'F' = Faces

if serial_test == 0
subj = upper(input('Please type subject initials and press return. ','s'));
end
if serial_test == 1
subj = test_initials;
end

% Ask which version they want to run (occ vs nocc)..
whichVer = input('Do you want to run the occluded or unoccluded version? Type "occ" for occluded, "nocc" for unoccluded, or "both" for both. Then press return.  ','s');
% Don't let them pass until they enter either a 'occ' or an 'nocc'
while (whichVer ~= "occ") && (whichVer ~= "nocc") && (whichVer ~= "both")
whichVer = input('Input not recognized. Please type either "occ" for occluded or "nocc" for unoccluded. Then press return.  ','s');
end

% GET USER INPUT FOR STIMULUS PARAMETER VERSION
%==========================================================================
% Initialize variables
L = 0;
T = 0;
C = 0;

% Loop until input is valid
validInput = false;
while ~validInput
    disp(" ");
    disp('Specify which parameters you want to be included/varied (Luminance, Texture, Common Region, or some combo)');
    userInput = input('To do so, enter a string containing only L, T, C, or some combo of these characters: ', 's');
    
    % Check for invalid characters
    if all(ismember(userInput, ['L', 'T', 'C']))
        validInput = true;
        
        % Check for each character and set variables
        if contains(userInput, 'L')
            L = 1;
        else
            L = 0;
        end
        
        if contains(userInput, 'T')
            T = 1;
        else
            T = 0;
        end
        
        if contains(userInput, 'C')
            C = 1;
        else
            C = 0;
        end
    else
        fprintf('Invalid input. Please enter a string containing only L, T, or C.\n');
    end
end
% Display variables for confirmation
fprintf('L = %d\n', L);
fprintf('T = %d\n', T);
fprintf('C = %d\n', C);
disp(" ");
%==========================================================================

% Ask which noise they want..
whichNoise = input('Which noise do you want to use? Type "white" for white noise or "kernel" for kernel noise. Then press return.  ','s');
% Don't let them pass until they enter either a 'white' or an 'kernel'
while (whichNoise ~= "white") && (whichNoise ~= "kernel")
whichNoise = input('Input not recognized. Type "white" for white noise or "kernel" for kernel noise. Then press return.  ','s');
end
disp(" ");
if whichNoise=="krnlNz"
noizType="krnlNz";
end
if whichNoise=="white"
noizType="white";
end


% Ask if they want to run the practice trials..
run_pracTrls = input('Do you want to run practice trials before the experiment? Press "y" for yes or "n" for no. Then press return.  ','s');
% Don't let them pass until they enter either a 'y' or an 'n'
while (run_pracTrls ~= "y") && (run_pracTrls ~= "n")
run_pracTrls = input('Input not recognized. Please press either "y" for yes or "n" for no. Then press return.  ','s');
end
disp(" ");

db_mode = 1;
db_mode_screen = 1;
dbm_skip = 1;

quest_db = 1;
start_dir = pwd;

%hard-code-path
%pract_FM_im ="/Users/testadmin/Documents/tDCS_stim_MATLAB/images/morph_levels.png"; % full path to facemorph example image for practice mode..
pract_FM_im =strcat(main_dir,"/images/morph_levels.png"); % full path to facemorph example image for practice mode.. % generalized..

%Ethan 100 im deck..
%image_dir ="/Users/eduwell/OneDrive - mcw.edu/Documents/SNAP/projects/Priyanka_Faces/quest_dev/faces/out/test3_frames_28-07-2022-01-40-29";

%EJD NOTE: WILL NEED TO GENERALIZE CODE TO WORK FOR ALL REVCOR VERSIONS
%W/OUT NEEDING TO ALTER PARS BELOW MANUALLY.. (4/4/23)
%##########################################################################
%##########################################################################

%noise_dir = strcat(main_dir,"/noise/white_test1");
%noise_dir = strcat(main_dir,"/noise/whiteBin5_rs0-15");

% Tag for finding/selecting base images within the base image directory...
imtag = ".png"; % unique tag in filename iding the images you want...

% Path to Noise Directory..
%noizType="krnlNz";
%noizType="krnlNz";

% if noizType=="white"
% noise_dir = strcat(main_dir,"/noise/512by512_whiteNoise_60000frms");
% nTag=".png";
% end

if noizType=="white"
    if L==1 && T ==0 && C==0
        noise_dir = strcat(main_dir,"/noise/512by512_whiteNoise_20000frms_smpl1");
        nTag=".png";
    end
    if T==1 && L==0 && L==0
        noise_dir = strcat(main_dir,"/noise/512by512_whiteNoise_20000frms_smpl2");
        nTag=".png";
    end
    if C==1 && L==0 && T==0
        noise_dir = strcat(main_dir,"/noise/512by512_whiteNoise_20000frms_smpl3");
        nTag=".png";
    end
end

if noizType=="krnlNz"
    if L==1 && T ==0 && C==0
        noise_dir = strcat(main_dir,"/noise/lumOnlyBI_20000frms_krnlNz_imzPerKrnl_111111100");
        nTag=".png";
    end
    if T==1 && L==0 && L==0
        noise_dir = strcat(main_dir,"/noise/texOnlyBI_20000frms_krnlNz_imzPerKrnl_111111100");
        nTag=".png";
    end
    if C==1 && L==0 && T==0
        noise_dir = strcat(main_dir,"/noise/crOnlyBI_20000frms_krnlNz_imzPerKrnl_111111100");
        nTag=".png";
    end
end

% noise_dir = strcat(main_dir,"/noise/lumBI_krnlNz_imzPerKrnl_111111000");
% % Tag for finding/selecting noise images within the directory above..
% nTag=".png";

%                      QUEST OPTIMIZED PARAMETER:
%##########################################################################

%  ---------  Difference Parameter Between-Objects Scenarios:  ---------  
% -------------------------------------------------------------------------

if L==1 && T ==0 && C==0
stimVer="lum";
% LUMINANCE ONLY
image_dir = strcat(main_dir,"/images/test4_LumOnly-05-Mar-2024-10-12-48");
qstParStr="lumDiff"; % diff between lum1/lum2
parTagz = ["lum1","lum2"];
obj1Par="lum1";
obj2Par="lum2";
nWt=0.75;
end

if T==1 && L==0 && C==0
stimVer="tex";
% TEXTURE ONLY
image_dir = strcat(main_dir,"/images/test5_TexOnly-12-Mar-2024-18-48-17"); % generalized
qstParStr="texDiff"; % diff between texture resolution1/texture resolution2
parTagz = ["Tr1","Tr2"];
obj1Par="Tr1";
obj2Par="Tr2";
nWt=0.75;
end
% -------------------------------------------------------------------------

%   ---------  Individual Parameter Single-Object Scenarios:  ---------  
if C==1 && L==0 && T==0
stimVer="cr";
% COMMON REGION ONLY
image_dir = strcat(main_dir,"/images/test5_CROnly-11-Mar-2024-14-16-09"); % generalized
qstParStr="CRpl"; % common region boundary luminance 
parTagz = ["CRpl"];
obj1Par="CRpl";
obj2Par="CRpl";
nWt=0.75;
end
% -------------------------------------------------------------------------

%     --------------  Noise Optimization Scenarios:  --------------  
% image_dir = strcat(main_dir,"/images/CRLumConst4-22-Apr-2023-19-24-24"); % generalized
% noise_dir = strcat(main_dir,"/noise/Blobs_test1");
% qstParStr="nWt"; % common region boundary luminance 
% parTagz = ["nWt"];
% obj1Par="lum1";
% obj2Par="lum2";
% nWt=0.50;

% -------------------------------------------------------------------------

%##########################################################################



%##########################################################################
%##########################################################################

% Quest Parameters
run_quests = 1; % specifies whether quests will be run (1) or not (0)
n_int = 3; % specify the number of interleved quests you want.. if 1, 
% no interleaving will occur.. only threshold 1 will be used. Can range
% from 0-5.. 5 is the max # of quests you can currently interleave .. if
% you need more, we'll need to alter Facequest5_int.m ...

%quests_ntrials = 52; % specify number of trials in each quest test rep
quests_ntrials_init = 50; % specify number of trials in each quest test rep on the first quest block

quests_ntrials_reg = 100; % specify number of trials in each quest test rep  on the quest block after the first block

quests_ntrials_practice = 6; % specify number of trials in each quest test rep

%t_prob = 0.6667; % probablility/proportion of ape-human trials (vs. ape only) (FOR THE ACTUAL EXPERIMENT)
%t_prob_practice = 0.6667; % probablility/proportion of ape-human trials (vs. ape only) (FOR THE PRACTICE SESSION)
t_prob = 0.5; % probablility/proportion of right vs. left trials (FOR THE ACTUAL EXPERIMENT)
t_prob_practice = 0.5; % probablility/proportion of right vs. left trials (FOR THE PRACTICE SESSION)

%Set different max/min parameter limits based on experiment type...
if qstParStr=="lumDiff"
maxmin = cell(1,2); % initialize maxmin
maxmin{1,1} = 255; % specifies the maximum parameter value allowed in quest..
maxmin{1,2} = 0; % specifies the minimum parameter value allowed in quest..
maxmin_seg = 0; % if 1, will segment the curve into n_int pieces and set
%separate maxs/mins for each segment to try to ensure even coverage..
% NOTE: IF YOU SELECT MAXMIN_SEG, MAKE SURE YOUR PTHRESHOLDS ARE IN
% ASCENDING ORDER!!!
end
if qstParStr=="texDiff"
maxmin = cell(1,2); % initialize maxmin
maxmin{1,1} = 6.25; % specifies the maximum parameter value allowed in quest..
maxmin{1,2} = 0; % specifies the minimum parameter value allowed in quest..
maxmin_seg = 0; % if 1, will segment the curve into n_int pieces and set
% separate maxs/mins for each segment to try to ensure even coverage..
% NOTE: IF YOU SELECT MAXMIN_SEG, MAKE SURE YOUR PTHRESHOLDS ARE IN
% ASCENDING ORDER!!!
end
if qstParStr=="CRpl"
maxmin = cell(1,2); % initialize maxmin
maxmin{1,1} = 127.5; % specifies the maximum parameter value allowed in quest..
maxmin{1,2} = 0; % specifies the minimum parameter value allowed in quest..
maxmin_seg = 0; % if 1, will segment the curve into n_int pieces and set
%separate maxs/mins for each segment to try to ensure even coverage..
% NOTE: IF YOU SELECT MAXMIN_SEG, MAKE SURE YOUR PTHRESHOLDS ARE IN
% ASCENDING ORDER!!!
end

%Condition Pars
%Specify occlusion condition
if whichVer=="occ"
olCon=0; %(0=occluded, 1=non-occuluded, 2=both)
elseif whichVer=="nocc"
olCon=1; %(0=occluded, 1=non-occuluded, 2=both)
elseif whichVer=="both"
olCon=2; %(0=occluded, 1=non-occuluded, 2=both)
end

if maxmin_seg == 1
maxmin_seg_str = "_mmSeg";
else
maxmin_seg_str = "";
end

% Interleaved staircase threshold parameters..
% NOTE: DEPENDING ON THE # OF INTERLEAVED QUESTS SELECETED (WITH N_INT
% ABOVE) SOME OF THESE THRESHOLDS MAY OR MAY NOT BE USED.. IT USES THE
% PTHRESHOLDS 1-N_INT...
%Staircase1
pThreshold1=0.65; % NOTE: make this your threshold of interest.. this is the one compared on the graphs..
%Staircase2
pThreshold2=0.75;    
%Staircase3
pThreshold3=0.85;
%Staircase4
pThreshold4=0.85;
%Staircase5
pThreshold5=0.85;
% Load the above into a structure...
pthrds = struct();
pthrds.thr1 = pThreshold1;
pthrds.thr2 = pThreshold2;
pthrds.thr3 = pThreshold3;
pthrds.thr4 = pThreshold4;
pthrds.thr5 = pThreshold5;

DataDir=strcat(main_dir,"/data_master");
verStrVect=[noizType,stimVer,whichVer];

qparsOut = getPrevQuestPars(DataDir,subj,n_int,pthrds,verStrVect);

qst_startParams = struct(); % preallocate for Quest initial starting params..
if qparsOut.isPrevQuest==0
    disp(" ");
    disp("Using default QUEST starting parameter values.");
    disp(" ");

    quest_testreps=quest_testreps_init;
    quests_ntrials=quests_ntrials_init;

% if the qparsOut.isPrevQuest output of  getPrevQuestPars() is 0, there is no previous quest
% data.. in this case set the default/"first run" parameters..
for hh=1:n_int
    qStringTmp=strcat("q",num2str(hh));
%     qst_startParams.(qStringTmp).tGuess = 0.01; 
%     qst_startParams.(qStringTmp).tGuessSd = 20;  % std. dev. on initial guess of threshold
%     qst_startParams.(qStringTmp).grain = 0.01;
%     qst_startParams.(qStringTmp).range = 1000;
    qst_startParams.(qStringTmp).tGuess = maxmin{1,1}*0.01; % start at 1% of the max
    qst_startParams.(qStringTmp).tGuessSd = 0.3*(abs(maxmin{1,1}-maxmin{1,2}));  % set std. dev. to n% of the total possible parameter value range
    qst_startParams.(qStringTmp).grain = 0.01;
    qst_startParams.(qStringTmp).range = 1.5*(abs(maxmin{1,1}-maxmin{1,2})); % set quest range to 2 times total possible parameter range.
end

else
% if the qparsOut output of  getPrevQuestPars() is not 0, there is a previous
% quest session available.. if so, use this to set the starting parameters
% for each of the interleaved quests..

quest_testreps=quest_testreps_reg;
quests_ntrials=quests_ntrials_reg;

for hh=1:n_int
    qStringTmp=strcat("q",num2str(hh));
    qst_startParams.(qStringTmp).tGuess = qparsOut.pthrds(2,hh);   % pull initial thr estimate/guess from previous quest thr out
    %qst_startParams.(qStringTmp).tGuessSd = qparsOut.sdevs(1,hh);  % pull initial std. dev guess from previous quest sd out
    %qst_startParams.(qStringTmp).tGuessSd = (qparsOut.sdevs(1,hh))*5;  % pull initial std. dev guess from previous quest sd out.. multiply by 5 to make sure its able to move around enough at first..
    qst_startParams.(qStringTmp).tGuessSd = 0.3*(abs(maxmin{1,1}-maxmin{1,2}));  % set std. dev. to n% of the total possible parameter value range
    qst_startParams.(qStringTmp).grain = 0.01; % may want to make these tailored in future.. for now keep as is
    qst_startParams.(qStringTmp).range = 1.5*(abs(maxmin{1,1}-maxmin{1,2})); % set quest range to 2 times total possible parameter range.
end
end

repcheck = input(strcat("The total number of blocks/QUEST session repetitions is curently set to ",num2str(quest_testreps),".\n If that looks OK, type 'y' and press return. If not, type 'n' and then press return.  "),'s'); 
if repcheck == 'n' 
    disp("NOTE: IF you run more than one QUEST session rep, subsequent stimulus programs (experiment version) will use the average psychometric function from across all reps..")
    quest_testreps = str2double(input("Type the number of reps you wish to run. Then press return.  ",'s')); 
    disp(" ");
end

% Naka Rushton Automated Response Psychometric Curve Params:
nr_rmax = 0.5;
nr_c50 = 115;
nr_b = 0.5;
nr_n = 2;

% Load the above into a structure..
nr_params = struct();
nr_params.nr_rmax = nr_rmax;
nr_params.nr_c50 = nr_c50;
nr_params.nr_b = nr_b;
nr_params.nr_n = nr_n;

%NR model bounds
upper = maxmin{1,1}; 
lower = 0;
%incrmt = 1; % controls the x axis resolution/smallest incriment..
incrmt=0.1;
%c50fstart = 100; % specifies starting value in c50 optimization for NR function fit
c50fstart = (upper+lower)/2; % specifies starting value in c50 optimization for NR function fit
c50_upper = upper; % specifies upper limit value of c50 for NR function fit
c50_lower = lower; % specifies lower limit value of c50 for NR function fit
c50_initGrid = 1;
c50fstart_grid = linspace(c50_lower,c50_upper,10);

% Optimizable NR parameters.. (control which params are optimized/which are fixed..)
nr_opt = 1; % if 0, only c50 is allowed to vary freely in the optimization/fitting process.., if 1, others specified below are also optimized..
nr_opt_shown2 = 0; % if 1, will fit the version of the nr curve without n optimized too for comparison..
if nr_opt == 1
nr_n_upper = 10; % specifies the upper limit of the "n" parameter (exponent) in the naka rushton model
nr_n_lower = 2; % specifies lower limit of nr_n..
end

% Bin threshold parameters for binning and including datapoints on model
% fitting...
%n_bins = 100; % specifies the number of bins..
n_bins = 50; % specifies the number of bins..

excl_bin1 = 0; % if 1, excludes the first bin from the model fitting, if 0, includes all bins..
% above option was added when Ethan realized that the first bin (ape face
% only) would likely/potentially not conform with the psychometric curve..

%n_thresh = round((quests_ntrials*n_int)*0.04); % specifies the # of datapoints needed to be included as a "bin"
n_thresh = 4;

bin_weight = 1; % if 1, bin means will be "weighted" based on their N... ie if N=10, that bin's mean will be included 10x..
bweight_exp =1; % if>1, bin means will be weighted exponentially by N^bweight_exp...

% Session Parameters
run_session = 0; % Specifies whether sessions will be run in addition to 
% quests.. (CURRENTLY A PLACEHOLDER.. ONLY DEVELOPED QUESTS AS OF NOW..)

%hard-code-path
%figs_dir_base = "/Users/testadmin/Documents/tDCS_stim_MATLAB/data_master/";
figs_dir_base = strcat(main_dir,"/data_master/"); % generalized


%figs_dir=stimVer;
figs_dir=join(verStrVect,"/"); % combine version par strings with "/" to form subdirectory to contain data..

%figs_dir = strcat("bins",num2str(n_bins),"_binTh",num2str(n_thresh),"_nqint",num2str(n_int),"_ntrls",num2str(quests_ntrials),"_thrs",thr_str,"_grn",grain_str,"_tgst",tguess_str,"_rng",range_str,"_sdG",tGuessSd_str,"_tprob",num2str(t_prob*100),"_c50u",num2str(c50_upper),"_c50l",num2str(c50_lower),"_c50st",num2str(c50fstart),bweight_str,"_bwExp",num2str(bweight_exp),"_maxmin",num2str(maxmin{1,1}),"_",num2str(maxmin{1,2}),maxmin_seg_str);

%% Screen Setup
try
    %%% some standard things that we've got to set up
    Screen('Preference','Tick0Secs',nan);
    %ListenChar(2);  % stop throwing characters to matlab windows     %EJD commented
    %HideCursor;                                                      %EJD commented
    GetSecs;						% Pre-load and pre-allocate
    CharAvail;						% Pre-load and pre-allocate
    KbName('UnifyKeyNames');
    KbCheck(-1);
    
    %%% choose appropriate resolution
    critstimdisplay = 100;
    rwidth = 1024;	% requested resolution width
    rheight = 768;	% requested resolution height
    
    screenrect = Screen('Rect',0);	% get the size of the display screen
    zoom_factor = screenrect(4)/rheight;
    
%     res1=Screen('Resolutions', 0); %list of all resolutions possible
%     % grab the screen parameters closest to ones requested
%     % suggested by from quest (param1)
%     a=[res1(:).width];
%     b=[res1(:).height];
%     n=rheight; % assign the most recent param in param1 to "n"
%     [val,idx]=min(abs(b-n));
%     closest_param=b(idx);
%     rheight = closest_param;
%     rwidth = a(idx);
% 
%     res2=find([res1(:).width]==rwidth & [res1(:).height]==rheight);
%     res3=unique([res1(res2).hz]);
%     for r=1:length(res3)
%         rmod(r)=mod(critstimdisplay,100*(1/res3(r)));  % match frame rate to value stored in critstimdisplay
%     end
%     
%     zoom_factor = screenrect(4)/rheight;
%     [x,rmin]=min(rmod); % best frequency match
%     res4=find([res1(:).width]==rwidth & [res1(:).height]==rheight & [res1(:).hz]==res3(rmin)); % best freq. & size
%     [x,rmax]=max([res1(res4).pixelSize]);   % largest pixelSize value of good resolutions
%     newres = res1(res4(rmax));
    %oldResolution=Screen('Resolution', 0, newres.width, newres.height,newres.hz, newres.pixelSize); %set it % EJD COMMENTED
    
    %screenrect = Screen('Rect',0);	% get the size of the display screen
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
%         directory = '/Users/agreenb/Documents/Files/CMU/PORC1/Data/';
        directory = '/Users/agreenb/Dropbox/OBA/PORC4/Data/';
        kboard = GetKeyboardIndices([],[],73400320);    % internal laptop keyboard
        room = 'office';
    elseif strcmp(comp.machineName,'AdamsAppleMB')	% this is the older (2007) Intel MacBook Pro
        directory = '/Users/agreenb/Dropbox/OBA/PORC4/Data/';
        kboard = GetKeyboardIndices([],[],487587840);    % external laptop keyboard (right USB location)
%       kboard = GetKeyboardIndices([],[],????);    % external laptop keyboard (left USB location)
        room = 'office';
    elseif strcmp(comp.machineName,'kern-mbp2-2021')	% Ethan's Macbook from Kern...
        %directory = '/Users/eduwell/OneDrive - mcw.edu/Documents/SNAP/projects/Rev_Corr/code/RevCorr/working_dir/data/';
        directory = figs_dir_base;
        kboard = GetKeyboardIndices();    % external laptop keyboard (right USB location)
%       kboard = GetKeyboardIndices([],[],????);    % external laptop keyboard (left USB location)
        room = 'office';
        Screen('Preference', 'SkipSyncTests', 1); %EJD: added to skip sync tests
        %oldEnableFlag=Screen('Preference', 'EmulateOldPTB'); %EJD added to try to deal with compatability issues with old psychtoolbox...
    elseif strcmp(comp.machineName,'NRLG-TBRC-12455')	% Priyanka's Macbook for tdCS...
        %directory = '/Users/eduwell/OneDrive - mcw.edu/Documents/SNAP/projects/Rev_Corr/code/RevCorr/working_dir/data/';
        directory = figs_dir_base;
        kboard = GetKeyboardIndices();   
%       kboard = GetKeyboardIndices([],[],????);    
        room = 'office';
        Screen('Preference', 'SkipSyncTests', 1); %EJD: added to skip sync tests
        %oldEnableFlag=Screen('Preference', 'EmulateOldPTB'); %EJD added to try to deal with compatability issues with old psychtoolbox...   
    elseif strcmp(comp.machineName,'tron')	% Priyanka's Macbook for tdCS...
        %directory = '/Users/eduwell/OneDrive - mcw.edu/Documents/SNAP/projects/Rev_Corr/code/RevCorr/working_dir/data/';
        directory = figs_dir_base;
        kboard = GetKeyboardIndices();   
%       kboard = GetKeyboardIndices([],[],????);    
        room = 'office';
        Screen('Preference', 'SkipSyncTests', 1); %EJD: added to skip sync tests
        Screen('Preference','VisualDebugLevel', 0);
        Screen('Preference','SuppressAllWarnings', 1);
    elseif strcmp(comp.machineName,'smog')	% Priyanka's Macbook for tdCS...
        %directory = '/Users/eduwell/OneDrive - mcw.edu/Documents/SNAP/projects/Rev_Corr/code/RevCorr/working_dir/data/';
        directory = figs_dir_base;
        kboard = GetKeyboardIndices();   
%       kboard = GetKeyboardIndices([],[],????);    
        room = 'office';
        Screen('Preference', 'SkipSyncTests', 1); %EJD: added to skip sync tests
        Screen('Preference','VisualDebugLevel', 0);
        Screen('Preference','SuppressAllWarnings', 1);
        %oldEnableFlag=Screen('Preference', 'EmulateOldPTB'); %EJD added to try to deal with compatability issues with old psychtoolbox...   
    elseif strcmp(comp.machineName,'psy-agreenb-r7')
        directory = figs_dir_base;
        kboard = GetKeyboardIndices();    
        room = 'TBRC_C2026_1'; % subject testing room 2026, machine closest to door..
        Screen('Preference', 'SkipSyncTests', 1); %EJD: added to skip sync tests
        Screen('Preference','VisualDebugLevel', 0);
        Screen('Preference','SuppressAllWarnings', 1);  
    else
        
       fprintf('******************************\nYou''re not on a known machine!\n******************************\n');
       
       junk = input('However... If you would like to proceed, type Y, if not press N, and the program will terminate.. ','s');
        
        if junk == 'Y'
            fprintf('******************************\nOKIE DOKIE HERE WE GO....KIE!\n******************************\n');
        else
        return;
        end
        
    end
    
    %% Start Stimulus Code
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create onscreen window
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    screenrect = Screen('Rect',0);	% get the size of the display screen
    
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

    
    window = PsychImaging('OpenWindow', 0, GrayIndex(0,0.5), screenrect);	% Open generic on-screen window
    Screen(window,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);    % allow for transparency!
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'UseFastOffscreenWindows');

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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %today = datestr(now,15);	% get current time/date stamp for filename
    today = datestr(now,30);
    %%% must begin by figuring out whether this is day1 or not
    olddir = cd;
    cd(directory); % change to data directory
    if not(isfolder(subj))
        mkdir(subj)
    end
    cd(subj)
    if not(isfolder(figs_dir))
        mkdir(figs_dir)
    end
    cd(figs_dir)
    todays_dir = strcat(subj,"_",today,"_QST");
    mkdir(todays_dir)
    cd(todays_dir)
    out_dir = pwd;
    fileid = strcat(out_dir,"/",subj);
    save(fileid,'today','link',"-v7.3");
    
    day=1;  % index for determining which day of testing this is for this subj
    
%     if exist(strcat(subj,'1.mat'),'file')
%         %there already exists a day1 file with that Ss ID
%         load(strcat(subj,'1.mat'),'day','quest*','stair'); % load the order of conditions, day number, threshold, & stairs status for this subject
%         day = day + 1;  % advance the day counter
%         fileid1 = strcat(directory,subj,num2str(day));
%         save(fileid1,'day','-append');  % save the new day number
%         fileid = strcat(directory,subj,num2str(day));
%         save(fileid,'today');
%     else
%         %this must be day1 for this subject
%         day = 1; %EJD UNCOMMENTED
%         %stair = 'N';  %EJD UNCOMMENTED
%         fileid = strcat(directory,subj,'1');
%         save(fileid,'today','link');
%     end
    

    % save room location for each day
    eval(['room' num2str(day) ' = room;']);
    eval(['save(fileid,''room' num2str(day) ''',''-append'');']);
    

    if run_quests==1
        cd(olddir); % change back to code dir
            
    % initialize output arrays..
    trialtm_out = cell(1,quest_testreps);
    ts_out_out = cell(1,quest_testreps);
    qs_out_out = cell(1,quest_testreps);
    sd_out_out = cell(1,quest_testreps);
    nr_ideals_out = cell(1,quest_testreps);
    params_out_out = cell(1,quest_testreps);
    functions_out = cell(1,quest_testreps);
    fit_t1_out = cell(1,quest_testreps);
    figs_out = cell(1,quest_testreps);
    fit_curve_out = cell(1,quest_testreps);
    ver_curve_out = cell(1,quest_testreps);
    tdfs = cell(1,quest_testreps);
    
%% Run practice trials if requested..
    if run_pracTrls == "y"
    % run the practice session:
    %[tdf_out_practice] = RevCorrQuest5_int_v3_practice(room,kboard,quests_ntrials_practice,link,1,window,image_dir,db_mode,db_mode_screen,nr_params, imtag, pthrds,n_int,qst_startParams,t_prob_practice,maxmin,maxmin_seg,pract_FM_im,main_dir,olCon,parTagz,obj1Par,obj2Par,qstParStr,noise_dir,nWt,nTag);

    [~,~,~,~,~, ~,tdf_out_practice] = RevCorrQuest5_int_v6_practice(room,kboard,quests_ntrials_practice,link,1,window,image_dir,db_mode,db_mode_screen,nr_params, imtag, pthrds,n_int,qst_startParams,t_prob,maxmin,maxmin_seg,main_dir,olCon,parTagz,obj1Par,obj2Par,qstParStr,noise_dir,nWt,nTag,out_dir);
    else
        tdf_out_practice = [];
    end

%% Now run the actual quests...
    for vv = 1:quest_testreps
                
                if run_pracTrls == "y"
                    giveInstrxn=0;                
                elseif ((run_pracTrls == "n") && (vv==1))
                    giveInstrxn=1;
                elseif ((run_pracTrls == "n") && (vv>1))
                    giveInstrxn=0;
                end
                
                [trialtm,ts_out,qs_out,sd_out,nr_ideals, params_out,tdf_out] = RevCorrQuest5_int_v7(room,kboard,quests_ntrials,link,1,window,image_dir,db_mode,db_mode_screen,nr_params, imtag, pthrds,n_int,qst_startParams,t_prob,maxmin,maxmin_seg,main_dir,olCon,parTagz,obj1Par,obj2Par,qstParStr,noise_dir,nWt,nTag,out_dir,giveInstrxn);
                
                % Create "waiting" Page
                Screen('FillRect', window, gray); % gray out window
                % add formatted text..
                txtDistFrmTop1=50*scrn_1percY; % specify distance of text from top as percent of screen dimension..
                Screen('TextSize',window, 50);
                speech1 = ['Saving Data. Please wait... \n'];
                DrawFormattedText(window, speech1,'center', (scrn_top+txtDistFrmTop1), black);
                Screen('Flip', window);% show linking display

                %[trialtm,ts_out,qs_out,sd_out,nr_ideals, params_out] = FACEquest5_int(room,kboard,quests_ntrials,link,0,window,image_dir,db_mode,db_mode_screen,nr_params, imtag, pthrds,n_int,qst_startParams,t_prob,maxmin,maxmin_seg);
                save(strcat(fileid,".mat"),"-append"); % save entire workspace after each block just to be to be safe..
%% Process outputs from Quest                
                % Unpack the outputs from the quests
                % quest estimated thresholds
                st1_1_qst = ts_out(1);
                if n_int>1
                    t2_1_qst = ts_out(2);
                    if n_int>2
                        t3_1_qst = ts_out(3);
                        if n_int>3
                            t4_1_qst = ts_out(4);
                            if n_int>4
                                t5_1_qst = ts_out(5);
                            end
                        end
                    end
                end

                % quest structures
                q1_1_qst = qs_out(1);               
                if n_int>1
                    q2_1_qst = qs_out(2);
                    if n_int>2
                        q3_1_qst = qs_out(3);
                        if n_int>3
                            q4_1_qst = qs_out(4);
                            if n_int>4
                                q5_1_qst = qs_out(5);
                            end
                        end
                    end
                end 
                
                % standard deviations
                sd1_1_qst = sd_out(1);
                if n_int>1
                    sd2_1_qst = sd_out(2);
                    if n_int>2
                        sd3_1_qst = sd_out(3);
                        if n_int>3
                            sd4_1_qst = sd_out(4);
                            if n_int>4
                                sd5_1_qst = sd_out(5);
                            end
                        end
                    end
                end 
 
                % "ideal" naka rushton parameters from automated
                % responses..
                nr_ideal_t1_qst = nr_ideals(1);
                if n_int>1
                    nr_ideal_t2_qst = nr_ideals(2);
                    if n_int>2
                        nr_ideal_t3_qst = nr_ideals(3);
                        if n_int>3
                            nr_ideal_t4_qst = nr_ideals(4);
                            if n_int>4
                                nr_ideal_t5_qst = nr_ideals(5);
                            end
                        end
                    end
                end

                
                paramq1 = params_out{1,1};
                paramqz=paramq1;
                if n_int>1
                    paramq2 = params_out{1,2};
                    paramqz=vertcat(paramqz,paramq2);
                    if n_int>2
                        paramq3 = params_out{1,3};
                        paramqz=vertcat(paramqz,paramq3);
                        if n_int>3
                            paramq4 = params_out{1,4};
                            paramqz=vertcat(paramqz,paramq4);
                            if n_int>4
                                paramq5 = params_out{1,5};
                                paramqz=vertcat(paramqz,paramq5);
                            end
                        end
                    end
                end 
                
                % Combine parampz in proper order by interleaving by
                % trial..
                paramqzInt=[];%initialize
                for yy=1:size(paramqz,2)
                    trialSetTmp=paramqz(1:end,yy)';
                    paramqzInt=horzcat(paramqzInt,trialSetTmp);
                end
                
                % Extract the raw response data from the interleaved quests
                % and fit naka rushton function to the raw data.. 
                % With this, then give the estimated psychometric thresholds 
                % from the fit NR model 
                q_list = who("q*_*_qst");
                for ii = 1:length(q_list)
                    
%                     cmd_str_qparam = strcat("q_param = ", "paramq",num2str(ii),"(1:(",q_list{ii},".trialCount));");
%                     cmd_str_qresp = strcat("q_resp = ", q_list{ii},".response(1:(",q_list{ii},".trialCount));");
                    cmd_str_qparam = strcat("q_param = ", "paramqzInt","(1:(",q_list{ii},".trialCount));");
                    cmd_str_qresp = strcat("q_resp = ", q_list{ii},".response(1:(",q_list{ii},".trialCount));");
                    eval(cmd_str_qparam);
                    eval(cmd_str_qresp);
                    %q1_resp = ii.response(1:(ii.trialCount));
                    q_temp = vertcat(q_param,q_resp);
                    if ii == 1
                        q_respmat = q_temp;
                    else
                        q_respmat = horzcat(q_respmat,q_temp);
                    end
                end
                clear ii      
                

                % bucket the data ..
                q_respmat = q_respmat';
                %q_respmat = abs(q_respmat);
                
                
                % Find the bin placement based on the X values
                [N,edges,bins] = histcounts(q_respmat(:,1),n_bins);

                bin_means=zeros(2,n_bins);
                for n = 1:n_bins
                    %bin_means(:,n) = mean(q_respmat(bins==n,:))';
                    bin_means(1,n) = mean(q_respmat(bins==n,1))';
                    bin_means(2,n) = mean(q_respmat(bins==n,2))';
                end
                   
                % Deal with NaNs.. (eliminate bins without any values)..
                bin_means_cl = []; % initialize by starting with first col
                bin_means_cl = vertcat(bin_means_cl,bin_means_cl);
                for ii=1:(length(bin_means))
                    if isnan(bin_means(1,ii)) == 1
                        % Do Nothing...
                    else
                        if N(ii) >= n_thresh % if this bucket has above n_thresh samples, include it..
                        tmp = bin_means(:,ii);
                        bin_means_cl = horzcat(bin_means_cl,tmp);
                        end
                    end
                end
                clear ii
%                 if (isnan(bin_means_cl(1,1)) == 1) || (N(1)<=n_thresh)
%                     bin_means_cl = bin_means_cl(:,2:end);
%                 end
                
                bin_means_cl2 = bin_means_cl;

                %for ii=1:length(bin_means_cl2)
                %    bin_means_cl2(2,ii) = (1-bin_means_cl2(2,ii));
                %end

                %do same for N.. (eliminate positions that are 0 or whose # of samples is below the threshold)
                N_noZero = [];
                for ii=1:(length(N))
                    if (N(1,ii) >= n_thresh) && (isnan(bin_means(1,ii)) ~= 1)
                        tmp = N(1,ii);
                        N_noZero = horzcat(N_noZero,tmp);
                    end
                end
                clear ii
                clear x
                
                % if bin_weight == 1, loop through and duplicate each
                % datapoint in bin_means_cl N number of times (the number
                % of datapoints in the the bin)
                if bin_weight == 1
                    bin_means_cl_wtd1 = [];
                    bin_means_cl_wtd2 = [];
                    bin_means_cl_wtd3 = [];
                    bin_means_cl_wtd = vertcat(bin_means_cl_wtd1,bin_means_cl_wtd2);
                    bin_means_cl_wtd = vertcat(bin_means_cl_wtd,bin_means_cl_wtd3);
                    
                    for ee=1:size(bin_means_cl,2)
                        n_tmp = N_noZero(ee);
                        bintmp1 = bin_means_cl(1,ee);
                        bintmp2 = bin_means_cl(2,ee);
                        bintmpN = N_noZero(1,ee);
                        bintmp = vertcat(bintmp1,bintmp2);
                        bintmp = vertcat(bintmp,bintmpN);
                            for tt = 1:((n_tmp)^(bweight_exp))
                                bin_means_cl_wtd = horzcat(bin_means_cl_wtd,bintmp);
                            end
                    end
                    bin_means_cl = bin_means_cl_wtd; % reassign bin_means_cl to be bin_means_cl_wtd
                end

                % Now set up the naka rushton model function and fit it to
                % the binned data
                
                if excl_bin1 == 1
                x = bin_means_cl(1,2:end); % EJD Changed to 2:end (from :) .. to eliminate the first bin (monkey only.. you get near 100% on that..)
                y = bin_means_cl(2,2:end); % EJD Changed to 2:end (from :) .. to eliminate the first bin (monkey only.. you get near 100% on that..)
                else
                x = bin_means_cl(1,1:end); % EJD Changed to 2:end (from :) .. to eliminate the first bin (monkey only.. you get near 100% on that..)
                y = bin_means_cl(2,1:end); % EJD Changed to 2:end (from :) .. to eliminate the first bin (monkey only.. you get near 100% on that..)
                end

                 % define "veridical" NR psychometric curve for automated
                 % responses based on nr parameters specified in params
                 % section..
                 syms c
                 nr_fncn(c) = nr_rmax*((c^nr_n)/(nr_c50^nr_n + c^nr_n)) + nr_b;

                %upper = max(x);
                %lower = min(x);

                interval = lower:incrmt:upper;
                
                nr_verid = double(nr_fncn(interval));
                
                % fit the Naka Rushton model to the data..
                
                if nr_opt == 0 || nr_opt_shown2 ==1 % If nr_opt == 0, only optimize c50/fix the exponent param
                % "n" at 2...
                
                if c50_initGrid == 1
                    c50_initGrid_gofMat = cell(3,size(c50fstart_grid,2));
                    for ii = 1:size(c50fstart_grid,2)
                        c50_initGrid_gofMat{1,ii} = ii;
                        c50fstart = c50fstart_grid(ii);
                        [f,gof] = fit(x',y',"nr_rmax*((x^nr_n)/(nr_c50^nr_n + x^nr_n)) + nr_b","StartPoint",[nr_b,c50fstart, nr_n, nr_rmax],'Lower',[nr_b, c50_lower, nr_n, nr_rmax],'Upper',[nr_b, c50_upper, nr_n, nr_rmax],'Robust', 'Bisquare');
                        c50_initGrid_gofMat{2,ii} = gof.sse;
                        c50_initGrid_gofMat{3,ii} = f;
                    end
                    gof_mat = cell2mat(c50_initGrid_gofMat(2,1:end));
                    [val,idx]=min(gof_mat);
                    f = c50_initGrid_gofMat{3,idx};
                    pause = "";
                    clear ii
                    
                    
                else
                f = fit(x',y',"nr_rmax*((x^nr_n)/(nr_c50^nr_n + x^nr_n)) + nr_b","StartPoint",[nr_b,c50fstart, nr_n, nr_rmax],'Lower',[nr_b, c50_lower, nr_n, nr_rmax],'Upper',[nr_b, c50_upper, nr_n, nr_rmax],'Robust', 'Bisquare');
                end

                end

                %----------------------------------------------------------
                if nr_opt == 1 % If nr_opt == 1

                
                if c50_initGrid == 1
                    c50_initGrid_gofMat = cell(3,size(c50fstart_grid,2));
                    for ii = 1:size(c50fstart_grid,2)
                        c50_initGrid_gofMat{1,ii} = ii;
                        c50fstart = c50fstart_grid(ii);
                        [f2,gof2] = fit(x',y',"nr_rmax*((x^nr_n)/(nr_c50^nr_n + x^nr_n)) + nr_b","StartPoint",[nr_b,c50fstart, (nr_n_upper + nr_n_lower)/2, nr_rmax],'Lower',[nr_b, c50_lower, nr_n_lower, nr_rmax],'Upper',[nr_b, c50_upper, nr_n_upper, nr_rmax],'Robust', 'Bisquare');
                        c50_initGrid_gofMat{2,ii} = gof2.sse;
                        c50_initGrid_gofMat{3,ii} = f2;
                    end
                    gof_mat2 = cell2mat(c50_initGrid_gofMat(2,1:end));
                    [val,idx]=min(gof_mat2);
                    f2 = c50_initGrid_gofMat{3,idx};
                    pause = "";
                    clear ii
                    
                    
                else
                f2 = fit(x',y',"nr_rmax*((x^nr_n)/(nr_c50^nr_n + x^nr_n)) + nr_b","StartPoint",[nr_b,c50fstart, (nr_n_upper + nr_n_lower)/2, nr_rmax],'Lower',[nr_b, c50_lower, nr_n_lower, nr_rmax],'Upper',[nr_b, c50_upper, nr_n_upper, nr_rmax],'Robust', 'Bisquare');
                end

                end
                %----------------------------------------------------------

%                 if nr_opt_shown2 == 1 % If nr_opt == 0, only optimize c50/fix the exponent param
%                 % "n" at 2...
%                 f2 = fit(x',y',"nr_rmax*((x^nr_n)/(nr_c50^nr_n + x^nr_n)) + nr_b","StartPoint",[nr_b,c50fstart, nr_n, nr_rmax],'Lower',[nr_b, c50_lower, nr_n, nr_rmax],'Upper',[nr_b, c50_upper, nr_n, nr_rmax],'Robust', 'Bisquare');
%                 end
% 
%                 
%                 if nr_opt == 1 % If nr_opt == 1, optimize both c50 and nr_n..
%                 f = fit(x',y',"nr_rmax*((x^nr_n)/(nr_c50^nr_n + x^nr_n)) + nr_b","StartPoint",[nr_b,c50fstart, (nr_n_upper + nr_n_lower)/2, nr_rmax],'Lower',[nr_b, c50_lower, nr_n_lower, nr_rmax],'Upper',[nr_b, c50_upper, nr_n_upper, nr_rmax],'Robust', 'Bisquare');
%                 end
                
%% Make Figs
                figure()
                hold
                counts=bin_means_cl(3,:);
                if nr_opt == 0
                fitted_curve = double(f(interval));
                end
                if nr_opt == 1
                fitted_curve = double(f2(interval));
                end

                s = scatter(x',y',25,counts,'filled');
                colorbar;
                %plot(x',y','o','Color','blue')
                plot(interval,fitted_curve,'Color','red');
                %plot(f,x',y')
                
                if nr_opt == 0
                fit_t1 = double(max(solve(f.nr_rmax*((c^f.nr_n)/(f.nr_c50^f.nr_n + c^f.nr_n)) + f.nr_b == pThreshold1,c)));
                end
                
                if nr_opt == 1
                    % define NR psychometric function for based on the parameters
                    % optimized in function f after the fitting above..
                    if mlb_year > 2021
                        syms cf [1 1] real
                    else
                        syms cf real
                    end
                    %nr_fncn_f(cf) = nr_rmax*((cf^nr_nf)/(nr_c50f^nr_nf + cf^nr_nf)) + nr_b;
                    
                    fit_t1 = double(max(solve(f2.nr_rmax*((cf^f2.nr_n)/(f2.nr_c50^f2.nr_n + cf^f2.nr_n)) + f2.nr_b == pThreshold1,cf)));
                    
                end

                if db_mode == 1
                nr_ideal_t1 = double(max(solve(nr_rmax*((c^nr_n)/(nr_c50^nr_n + c^nr_n)) + nr_b == pThreshold1,c)));
                perct1_fit_err = ((abs(nr_ideal_t1-fit_t1))/nr_ideal_t1)*100; 
                end

                th_str = num2str(pThreshold1);
                plot(fit_t1,pThreshold1,'o','LineWidth',3,'Color','red')
                
                if db_mode == 1
                %plot(interval,nr_verid,'blue')
                %plot(nr_ideal_t1,pThreshold1,'o','LineWidth',3, 'Color','blue')
                end
                
                if nr_opt == 1 && nr_opt_shown2 == 1
                fitted_curve2 = double(f(interval));
                plot(interval,fitted_curve2,'Color','green')
                end

                if db_mode == 1

                    if nr_opt == 0
                    %legend("data","fitted curve", strcat("fitted curve"," ",th_str," threshold"), "veridical curve", strcat("veridical"," ",th_str," threshold"))
                    legend("data","fitted curve", strcat("fitted curve"," ",th_str," threshold"))
                    end
                    if nr_opt == 1
                    %legend("data","fitted curve", strcat("fitted curve"," ",th_str," threshold"), "veridical curve", strcat("veridical"," ",th_str," threshold"), "fitted curve c50 opt. only..")
                    legend("data","fitted curve", strcat("fitted curve"," ",th_str," threshold"))
                    end
                elseif db_mode == 0
                legend("data","fitted curve", strcat("fitted curve"," ",th_str," threshold"))
                end
                
                xlim([lower upper])
                ylim([0.5 1])
                %idealNR_comparison_fig = gcf;
                
                
                %pack up outputs..
                trialtm_out{1,vv} = trialtm;
                ts_out_out{1,vv} = ts_out;
                qs_out_out{1,vv} = qs_out;
                sd_out_out{1,vv} = sd_out;
                nr_ideals_out{1,vv} = nr_ideals;
                params_out_out{1,vv} = params_out;
                if nr_opt==0
                functions_out{1,vv} = f;
                end
                if nr_opt==1
                functions_out{1,vv} = f2;
                end
                
                fit_t1_out{1,vv} = fit_t1;
                fit_curve_out{1,vv} = fitted_curve;
                ver_curve_out{1,vv} = nr_verid;
                tdfs{1,vv} = tdf_out;

                % Print out figs if requested..
                if print_figs == 1
                    %fig_fname = strcat(figs_dir,"_rep_",num2str(vv),".tif");
                    fig_fname = strcat(join(verStrVect,"_"),"_rep_",num2str(vv),".tif");
                    
                    cd(directory); % change to data directory
                    cd(subj)
                    cd(figs_dir)
                    cd(todays_dir)
                    frameTmp=getframe(gcf);
                    imwrite(frameTmp.cdata, fig_fname);
                    fig_readback = imread(strcat(fig_fname));
                    cd(start_dir)
                end
                
                figs_out{1,vv} = fig_readback;
                % Save the output parameters/data to file..
                if db_mode == 1
                    save(fileid,'trialtm_out','ts_out_out','qs_out_out','sd_out_out','nr_ideals_out','nr_params','params_out_out','functions_out','fit_t1_out','ver_curve_out','figs_out','-append');
                    save(fileid,'trialtm_out','ts_out_out','qs_out_out','sd_out_out','params_out_out','functions_out','fit_t1_out','fit_curve_out','tdfs','tdf_out_practice','figs_out','-append');
                    %save(fileid,'trialtm','ts_out','qs_out','sd_out','nr_ideals','nr_params','params_out','f','fit_t1','idealNR_comparison_fig','-append');
                end
                if db_mode == 0
                    save(fileid,'trialtm_out','ts_out_out','qs_out_out','sd_out_out','params_out_out','functions_out','fit_t1_out','fit_curve_out','tdfs','tdf_out_practice','figs_out','-append');                  
                end

                % Save/update info regarding which noise frames have been
                % used for this subject.
                %==========================================================
                %==========================================================
                strtDirTmp=pwd; % save start point as temp string var

                % change to data directory for this subject/noise/setup
                cd(directory); 
                cd(subj);

                usedNoiseFileName="usedNoise.mat";
                if isfile(usedNoiseFileName)
                    disp(" ")
                    disp(strcat(usedNoiseFileName," found. Loading data and appending list of noise frames used this run..."));
                    disp(" ")
                    usedNoise=load(usedNoiseFileName);

                    % save noise frame list, but make sure the column is
                    % correct/hasn't changed first.. check that the header
                    % of col 17 is "noiseImgFile"
                    if tdf_out{1,17}=="noiseImgFile"
                    usedFrmPrev=usedNoise.usedFrmList;
                    usedFrmList=vertcat(tdf_out(2:end,17),usedFrmPrev); % vertically concatenate this list onto the existing one
                    save(usedNoiseFileName,'usedFrmList','-append');            
                    else
                        disp(" ");
                        disp("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                        disp("!!! Noise Frame List Not Found in the Expected Column !!!");
                        disp("!!!     Could Not Save List of Used Noise Frames      !!!");
                        disp("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                        disp(" ");
                    end

                else
                    disp(" ")
                    disp(strcat(usedNoiseFileName," not found. Assuming this is the first run with this subject/setup..."));
                    disp(" ")
                    
                    % save noise frame list, but make sure the column is
                    % correct/hasn't changed first.. check that the header
                    % of col 17 is "noiseImgFile"
                    if tdf_out{1,17}=="noiseImgFile"
                    usedFrmList=tdf_out(2:end,17);

                    save(usedNoiseFileName,'usedFrmList');

                    else
                        disp(" ");
                        disp("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                        disp("!!! Noise Frame List Not Found in the Expected Column !!!");
                        disp("!!!     Could Not Save List of Used Noise Frames      !!!");
                        disp("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                        disp(" ");
                    end
                    
                end

                cd(strtDirTmp);
                %==========================================================
                %==========================================================
    end
    % Create "waiting" Page
    Screen('FillRect', window, gray); % gray out window
    % add formatted text..
    txtDistFrmTop1=50*scrn_1percY; % specify distance of text from top as percent of screen dimension..
    Screen('TextSize',window, 50);
    speech1 = ['You''ve reached the end of this set of blocks! \n'... 
        'Please wait while we wrap things up...'];
    DrawFormattedText(window, speech1,'center', (scrn_top+txtDistFrmTop1), black);
    Screen('Flip', window);% show linking display
    end

    %% Save output figs/gifs
    % compile the images in figs_out into a .gif file for rapid viewing if requested...
    if print_fig_gif == 1
        cd(directory); % change to data directory
        cd(subj)
        cd(figs_dir)
        cd(todays_dir)
        nImages = size(figs_out,2);
        %gif_filename = strcat(figs_dir,'_reps1-',num2str(nImages),'.gif'); % Specify the output file name
        gif_filename = strcat(join(verStrVect,"_"),'_reps1-',num2str(nImages),'.gif'); % Specify the output file name
        for idx = 1:nImages
            [A,map] = rgb2ind(figs_out{1,idx},256);
            if idx == 1
                imwrite(A,map,gif_filename,'gif','LoopCount',Inf,'DelayTime',gif_delay);
            else
                imwrite(A,map,gif_filename,'gif','WriteMode','append','DelayTime',gif_delay);
            end
        end
        cd(start_dir)
    end
    
    if mean_sd_fitCrv_fig == 1
        cd(directory); % change to data directory
        cd(subj)
        cd(figs_dir)
        cd(todays_dir)
        % preallocate space 
        all_fitCrvs = zeros(size(fitted_curve,1),quest_testreps);
        for ii = 1:quest_testreps
            all_fitCrvs(1:end,ii) = fit_curve_out{1,ii};
        end
        clear ii
        
        mean_fitCrv = zeros(size(fit_t1,1),1);
        stdev_fitCrv = zeros(size(fit_t1,1),1);
        min_fitCrv = zeros(size(fit_t1,1),1);
        max_fitCrv = zeros(size(fit_t1,1),1);

        for ii = 1:size(fitted_curve,1)
            mean_fitCrv(ii,1) = mean(all_fitCrvs(ii,1:end));
            stdev_fitCrv(ii,1) = std(all_fitCrvs(ii,1:end));
            min_fitCrv(ii,1) = min(all_fitCrvs(ii,1:end));
            max_fitCrv(ii,1) = max(all_fitCrvs(ii,1:end));
        end
        clear ii
        
        figure
        [lineOut, fillOut] = stdshade(all_fitCrvs',0.5,[0.5,0,0],interval,0);
        hold
%         plot(interval,nr_verid,'LineWidth',2,'Color','blue')
%         alpha 0.5
%         plot(nr_ideal_t1,pThreshold1,'o','LineWidth',2, 'Color','blue')
        alpha 0.5
        plot(interval,min_fitCrv,"--",'LineWidth',1,'Color','black')
        plot(interval,max_fitCrv,"--",'LineWidth',1,'Color','black')

        xlim([lower upper])
        ylim([0.5 1])
        %fig_fname = strcat(figs_dir,"_MeanStDev",".tif");
        fig_fname = strcat(join(verStrVect,"_"),"_MeanStDev",".tif");
        frameTmp=getframe(gcf);
        imwrite(frameTmp.cdata, fig_fname)
        %mean_sdfig = imread(strcat(fig_fname));
        cd(start_dir)
    end

    % to indicate that this quest session was successfully
    % completed, rename the quest directory/add the
    % "completed" tag...
    cd(directory); % change to data directory
    cd(subj)
    cd(figs_dir)
    curDir=pwd;
    save(strcat(fileid,".mat"),"-append"); % Save entire workspace at the end..
    rxivMatlabCode(main_dir,strcat(curDir,"/",todays_dir)); % Save archived copies of all matlab code in the main revCorrStim directory at end...
    movefile(diaryfile,todays_dir); % move command window diary to todays folder..
    movefile(todays_dir,strcat(todays_dir,"_cmpltd")) % we got to the end! mark this directory as completed...


%% Close things up.. (screens, etc.)
    Screen('FillRect', window, BlackIndex(window));
    ListenChar(0);
    ShowCursor;
    %     Screen('Resolution', 0, oldResolution.width, oldResolution.height, oldResolution.hz, oldResolution.pixelSize); %set it back
    Screen('CloseAll');
    cd(start_dir);
    if serial_test == 0
        clear;
    end
    close all

    
catch
%     Screen('Resolution', 0, oldResolution.width, oldResolution.height, oldResolution.hz, oldResolution.pixelSize); %set it back
    ListenChar(0);
    ShowCursor;
    Screen('CloseAll');
    rethrow(lasterror);
    cd(start_dir);
    save(strcat(fileid,".mat"),"-append"); % save entire workspace in the case of an error..
end
