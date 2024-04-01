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
clc
clear mex
clear all

% Serial Testing/Development Params..
serial_test = 1;
quest_testreps = 1;
test_initials = "XXX"; %phoney initials for testing...
print_figs = 1;
print_fig_gif = 1;
gif_delay = 0.1;
mean_sd_fitCrv_fig = 1;
%figs_dir = "/Users/eduwell/OneDrive - mcw.edu/Documents/SNAP/projects/Priyanka_Faces/quest_dev/quest_pics/test2_bin10_thrs70-75-84/";
% if serial_test == 0
%     quest_testreps = 1;
%     print_figs = 0;
% end

red = [255, 0, 0, 255];	% red
link = 'F'; % 'B'=boundary collinearity; 'L'=luminance; 'Z'=luminance/collinearity combined; 'T' = Texture (EJD ADDED); 'F' = Faces

%if serial_test == 0
subj = upper(input('Please type subject initials and press return. ','s'));
%end

% Ask if they want to run the practice trials..
run_pracTrls = input('Do you want to run practice trials before the experiment? Press "y" for yes or "n" for no. Then press return.  ','s');
% Don't let them pass until they enter either a 'y' or an 'n'
while (run_pracTrls ~= "y") && (run_pracTrls ~= "n")
run_pracTrls = input('Input not recognized. Please press either "y" for yes or "n" for no. Then press return.  ','s');
end


% if serial_test == 1
% subj = test_initials;
% end

db_mode = 0;
db_mode_screen = 0;
dbm_skip = 1;

quest_db = 1;
start_dir = pwd;

% if running "texture" version (ie link = "T") the image_in parameter below
% specifies the full path/filename pointing to the texture image
%image_in = '/Users/eduwell/OneDrive - mcw.edu/Documents/SNAP/projects/Rev_Corr/code/RevCorr/working_dir/fractal_tex1.png';

% image directory for face-ape images ...
%orig
image_dir = "/Users/testadmin/Documents/tDCS_stim_MATLAB/images/filtering_stim/"; %******
% output directory path 
out_dir = "/Users/testadmin/Documents/tDCS_stim_MATLAB/data_master"; % specifies root directory in which you want to save subjects data subdirectories..

% accuracy threshold based hard/easy parameter criteria params
% (use these ones if you want to select hard/easy pars (hPar & ePar)
% in on a subject by subject basis based on % accuracy thresholds measured
% via QUEST..
FMpar_qst_thr = 0.80; % accuracy threshold for determining max facemorph parameter images you want included in quest trials
FMpar_nqst_thr = 0.70; % accuracy threshold for determining highest facemorph parameter images you want included in non-quest trials
quest_thr = 0.75; % specify the threshold you want quest to home in on..

% Alternatively...
% manual select hard/easy parameter criteria
manual_ParSelect = 0; % Set equal to 1 to set parameters manually
if manual_ParSelect == 1
    FMpar_mat = [44,300]; % manually specify the lowest and highest facemorph parameter images you want included in the experiment (*aside from the "cannonical ape"*)
end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% If not using manual parameter selection mode, navigate to the subject's
% most recent quest session output directory and set easy/hard parameters
% based on percent accuracy thresholds (specified above) from psychometric 
% curve model:

if manual_ParSelect == 0
% Go into the subjects directory and load the .mat from running quest
% use the quest-generated psychometric function to pull desired thresholds 
% for "easy" and "hard" facemorph parameters
% Go into the subjects directory and load the .mat from running quest
% use the quest-generated psychometric function to pull desired thresholds 
% for "easy" and "hard" facemorph parameters
cd(out_dir)
cd(subj)

% Get a list of the subdirectories..
folder = pwd;
subDs = GetSubDirsFirstLevelOnly(folder);

%If there is only one subdirectory (as there probably should be..) enter
%that one.. if there are multiple, alert the user and make them pick which
%one they want at the command line
if size(subDs,2) == 1
    qst_parDir = subDs{1,1};
else
    disp('***** This subject has multiple subdirectories! *****')
    disp("              Subject's subdirectories:              ")
    disp("-----------------------------------------------------")
    for ii = 1:size(subDs,2)
        disp(subDs{1,ii})
    end
    disp("-----------------------------------------------------")
    disp("");
    
    qst_parDir = input('Copy and paste one of the the directories above that you wish to use. Then press return.  ','s');
end
cd(qst_parDir);

% get list of subdirectories and make list of previously completed quest output subdirectories

% Get a list of the subdirectories..
folder = pwd;
subDs = GetSubDirsFirstLevelOnly(folder);
subDs = subDs';
pattern = "_QST_cmpltd";
qst_subDs = {};
numFilesProcessed = 0;	% For fun, let's keep track of how many files we processed.
for k = 1 : size(subDs,1)
    % Get this folder.
    thisFolderName = subDs{k,1};
    % See if it contains our required pattern.
    if ~contains(thisFolderName, pattern, 'IgnoreCase', true)
        % Skip this file because the filename does not contain the required pattern.
        continue;
    end
    % The pattern is in the filename if you get here, so do something with it.
    %fprintf('Now processing %s\n', thisFileName);
    qst_subDs = vertcat(qst_subDs,thisFolderName);
    numFilesProcessed = numFilesProcessed + 1;	% For fun, let's keep track of how many files we processed.
end
fprintf('We found %d folders with %s in the name.\n',numFilesProcessed, pattern);
clear numFilesProcessed

% Select the most recent quest output folder to use..
qst_subDs = sortrows(qst_subDs,'descend'); % sort in descending order (folders contain timestamp in name)
qst_Dir = qst_subDs{1,1}; % grab the one on top..

% enter the quest directory and load in the .mat file
fprintf("Entering the most recent quest output directory.. \n")
cd(qst_Dir);
folder_tmp = pwd;
% Get a list of all files.
fileList = dir(fullfile(folder_tmp, '*.*'));
allFileNames = {fileList.name};
% Define what pattern you will need in your filenames.
pattern = ".mat";
numFilesProcessed = 0;	% For fun, let's keep track of how many files we processed.
matFiles = {};
for k = 1 : length(allFileNames)
    % Get this filename.
    thisFileName = {fullfile(fileList(k).folder, allFileNames{k})};
    % See if it contains our required pattern.
    if ~contains(thisFileName, pattern, 'IgnoreCase', true)
        % Skip this file because the filename does not contain the required pattern.
        continue;
    end
    % The pattern is in the filename if you get here, so do something with it.
    %fprintf('Now processing %s\n', thisFileName);
    matFiles = vertcat(matFiles,thisFileName);
    numFilesProcessed = numFilesProcessed + 1;	% For fun, let's keep track of how many files we processed.
end
if size(matFiles,1) <= 1
fprintf('We found %d file with %s in the name. \n',numFilesProcessed, pattern);
else
    fprintf("*************** !WARNING! ***************\n")
    fprintf("********** MULTIPLE .MAT FILES **********\n")
    fprintf("*********** should only be 1 ************\n")
    fprintf('We found %d files with %s in the name.\n',numFilesProcessed, pattern);
    fprintf("*****************************************\n")
end

fprintf("Loading in Data from the Most Recent QUEST Session.. \n")
load(matFiles{1,1})

nFuncns = size(functions_out,2); % check number of psychometric function fits from quest in functions_out
Funcn_thrMat = zeros(nFuncns,3); % preallocate space for threshold parameter values..
FMpar_mat = [0,0]; % preallocate for FMpar_mat
if nFuncns > 1
    fprintf('We found %d psychometric function fits from the previous QUEST session. \n',nFuncns);
    fprintf('Using the average threshold values across the %d to determine the range of parameters selectable on QUEST and non-QUEST trials.. \n', nFuncns);
    fprintf('The max (easiest) parameter possible on QUEST trials will be  the average %-2.2f threshold. The max (easiest) parameter possible on non-QUEST trials will be the average %-2.2f threshold. \n', FMpar_qst_thr, FMpar_nqst_thr);
    for ii = 1:nFuncns
        fcn_tmp = functions_out{1,ii};
        Funcn_thrMat(ii,1) = NRfcnThr_Comp(FMpar_qst_thr, fcn_tmp);
        Funcn_thrMat(ii,2) = NRfcnThr_Comp(FMpar_nqst_thr, fcn_tmp);
        Funcn_thrMat(ii,3) = NRfcnThr_Comp(quest_thr, fcn_tmp);
    end
    FMpar_mat(1) = mean(Funcn_thrMat(:,1));
    FMpar_mat(2) = mean(Funcn_thrMat(:,2));
    quest_thr_parVal = mean(Funcn_thrMat(:,3));
    fprintf('The max (easiest) parameter possible on QUEST trials (%-2.2f threshold) for this subject is %-2.2f. The max (easiest) parameter possible on non-QUEST trials (%-2.2f threshold) for this subject is %-2.2f. \n', FMpar_qst_thr,FMpar_mat(1),FMpar_nqst_thr,FMpar_mat(2));
    fprintf('On the "QUEST trials", the inline QUEST procedure will be homing in on the subject''s  %-2.2f threshold  \n For this subject, that threshold was computed to be a facemorph parameter of %-2.2f in the initial QUEST session. \n That will be the initial starting "guess" value fed to QUEST in this session.\n', quest_thr,quest_thr_parVal);
else
    fprintf('We found %d psychometric function fit from the previous QUEST session. \n',nFuncns);
    fprintf('Using the threshold values from that function to determine the range of parameters selectable on QUEST and non-QUEST trials.. \n');
    fprintf('The max (easiest) parameter possible on QUEST trials will be  the %-2.2f threshold. The max (easiest) parameter possible on non-QUEST trials will be the %-2.2f threshold. \n', FMpar_qst_thr, FMpar_nqst_thr);
    fcn_tmp = functions_out{1,1};
    Funcn_thrMat(1,1) = NRfcnThr_Comp(FMpar_qst_thr, fcn_tmp);
    Funcn_thrMat(1,2) = NRfcnThr_Comp(FMpar_nqst_thr, fcn_tmp);
    Funcn_thrMat(1,3) = NRfcnThr_Comp(quest_thr, fcn_tmp);
    FMpar_mat(1) = Funcn_thrMat(1,1);
    FMpar_mat(2) = Funcn_thrMat(1,2);
    quest_thr_parVal = Funcn_thrMat(1,3);
    fprintf('The max (easiest) parameter possible on QUEST trials (%-2.2f threshold) for this subject is %-2.2f. The max (easiest) parameter possible on non-QUEST trials  (%-2.2f threshold) for this subject is %-2.2f. \n', FMpar_qst_thr,FMpar_mat(1),FMpar_nqst_thr,FMpar_mat(2));
    fprintf('On the "QUEST trials", the inline QUEST procedure will be homing in on the subject''s  %-2.2f threshold  \nFor this subject, that threshold was computed to be a facemorph parameter of %-2.2f in the initial QUEST session. \nThat will be the initial starting "guess" value fed to QUEST in this session.\n', quest_thr,quest_thr_parVal);

end
else
    disp("************ Currently in Manual Parameter Section Mode ************")
    fprintf("This means that the range of parameter values for QUEST and \n" + ...
        "non-Quest trials (ie FMpar_minmax) is set manually as opposed to \n" + ...
        "being computed based on the subjects accuracy thresholds determined \n" + ...
        "in a previous QUEST proceedure. \n \n")
    fprintf('The max (easiest) parameter possible on QUEST trials will be %-2.2f. \nThe max (easiest) parameter possible on non-QUEST trials will be %-2.2f. \n', FMpar_mat(1), FMpar_mat(2));
end

disp("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
input('If the parameters above look reasonable, press return to continue.  ','s')
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%Ethan 100 im deck..
%image_dir ="/Users/eduwell/OneDrive - mcw.edu/Documents/SNAP/projects/Priyanka_Faces/quest_dev/faces/out/test3_frames_28-07-2022-01-40-29";

imtag = ".png"; % unique tag in filename iding the images you want... %******

% Quest Parameters
run_quests = 1; % specifies whether quests will be run (1) or not (0)
n_int = 1; % specify the number of interleved quests you want.. if 1, %******
% no interleaving will occur.. only threshold 1 will be used. Can range
% from 0-5.. 5 is the max # of quests you can currently interleave .. if
% you need more, we'll need to alter Facequest5_int.m ...

quests_ntrials = 120; % specify number of trials in each quest test rep  %****** (FOR REAL SESSION)
quests_ntrials_prac = 20; % specify number of trials in each quest test rep  %****** (FOR PRACTICE SESSION)

t_prob = 0.5; % probablility/proportion of ape-human trials (vs. ape only) %******
q_prob = 0.5;

maxmin = cell(1,2); % initialize maxmin
maxmin{1,1} = 800; % specifies the maximum parameter value allowed in quest..
maxmin{1,2} = 1; % specifies the minimum parameter value allowed in quest..
maxmin_seg = 0; % if 1, will segment the curve into n_int pieces and set
%separate maxs/mins for each segment to try to ensure even coverage..
% NOTE: IF YOU SELECT MAXMIN_SEG, MAKE SURE YOUR PTHRESHOLDS ARE IN
% ASCENDING ORDER!!!

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
% NOTE: make this your threshold of interest.. this is the one compared on the graphs..
pThreshold1=quest_thr; 

%Staircase2
% pThreshold2=0.75;    
% %Staircase3
% pThreshold3=0.85;
% %Staircase4
% pThreshold4=0.85;
% %Staircase5
% pThreshold5=0.85;

% Load the above into a structure...
pthrds = struct();

pthrds.thr1 = pThreshold1;

% pthrds.thr2 = pThreshold2;
% pthrds.thr3 = pThreshold3;
% pthrds.thr4 = pThreshold4;
% pthrds.thr5 = pThreshold5;

% Quest initial starting params.. 
% Set initial threshold guess to those computed in previous QUEST session
% Other start params? .. not sure whether we should change those..
qst_startParams = struct();
qst_startParams.tGuess = quest_thr_parVal; 
qst_startParams.tGuessSd = 20;  % std. dev. on initial guess of threshold
qst_startParams.grain = 0.01;
qst_startParams.range = 1000;

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
upper = 1000; 
lower = 0;
incrmt = 1; % controls the x axis resolution/smallest incriment..
c50fstart = 100; % specifies starting value in c50 optimization for NR function fit
c50_upper = 500; % specifies upper limit value of c50 for NR function fit
c50_lower = 0; % specifies lower limit value of c50 for NR function fit

% Optimizable NR parameters.. (control which params are optimized/which are fixed..)
nr_opt = 0; % if 0, only c50 is allowed to vary freely in the optimization/fitting process.., if 1, others specified below are also optimized..
nr_opt_shown2 = 1; % if 1, will fit the version of the nr curve without n optimized too for comparison..
if nr_opt == 1
nr_n_upper = 10; % specifies the upper limit of the "n" parameter (exponent) in the naka rushton model
nr_n_lower = 2; % specifies lower limit of nr_n..
end

% Bin threshold parameters for binning and including datapoints on model
% fitting...
n_bins = 1000; % specifies the number of bins..
excl_bin1 = 0; % if 1, excludes the first bin from the model fitting, if 0, includes all bins..
% above option was added when Ethan realized that the first bin (ape face
% only) would likely/potentially not conform with the psychometric curve..

%n_thresh = round((quests_ntrials*n_int)*0.04); % specifies the # of datapoints needed to be included as a "bin"
n_thresh = 4;

bin_weight = 1; % if 1, bin means will be "weighted" based on their N... ie if N=10, that bin's mean will be included 10x..
bweight_exp =2; % if>1, bin means will be weighted exponentially by N^bweight_exp...
%n_thresh = 2;

% Session Parameters
run_session = 0; % Specifies whether sessions will be run in addition to 
% quests.. (CURRENTLY A PLACEHOLDER.. ONLY DEVELOPED QUESTS AS OF NOW..)

% Build descriptive output directory name for serial tests
thr_str ="";
for kk = 1:n_int
    thr_tmp = strcat("pThreshold",num2str(kk));
    thr_tmp = eval(thr_tmp);
    %thr_tmp=num2str(round(thr_tmp,2));
    thr_tmp = num2str(round(thr_tmp*100));
    thr_tmp= thr_tmp(1:2);
    if kk == 1
        thr_str = strcat(thr_str,thr_tmp);
    else
        thr_str = strcat(thr_str,"-",thr_tmp);
    end
end
tguess_str = num2str(round(qst_startParams.tGuess,2));
tguess_str = tguess_str(end-1:end);
tGuessSd_str = num2str(qst_startParams.tGuessSd);% std. dev. on initial guess of threshold
grain_str = num2str(round(qst_startParams.grain,2));
grain_str = grain_str(end-1:end);
range_str = num2str(round(qst_startParams.range,2));
if bin_weight == 1
    bweight_str = "_binwtd";
else
    bweight_str = "";
end
figs_dir_base = "/Users/testadmin/Documents/tDCS_stim_MATLAB/data_master/";
figs_dir = strcat("bins",num2str(n_bins),"_binTh",num2str(n_thresh),"_nqint",num2str(n_int),"_ntrls",num2str(quests_ntrials),"_thrs",thr_str,"_grn",grain_str,"_tgst",tguess_str,"_rng",range_str,"_sdG",tGuessSd_str,"_tprob",num2str(t_prob*100),"_c50u",num2str(c50_upper),"_c50l",num2str(c50_lower),"_c50st",num2str(c50fstart),bweight_str,"_bwExp",num2str(bweight_exp),"_maxmin",num2str(maxmin{1,1}),"_",num2str(maxmin{1,2}),maxmin_seg_str);

% Build descriptive output directory 
cd(out_dir) % go to the base output dir..
cd(subj)
cd(qst_parDir)
currdir = pwd;

%get unique timestamp for directory
%time_vect = clock;
%unq_timeStamp = strcat(num2str(time_vect(3)),"-",num2str(time_vect(2)),"-",num2str(time_vect(1)),"-",num2str(time_vect(4)),"-",num2str(time_vect(5)),"-",num2str(round(time_vect(6))));
today = datestr(now,30);

% make unique directory name by concatenating subject's 3 digit id with
% unique time stamp..
out_dir_sub = strcat(subj,"_",today,"_TRN");
out_dir_final = strcat(currdir,"/",out_dir_sub,"/");
mkdir(out_dir_sub);
cd(start_dir);
%% Screen Setup
try
    %%% some standard things that we've got to set up
    Screen('Preference','Tick0Secs',nan);
    %ListenChar(2);  % stop throwing characters to matlab windows     %EJD commented
    %HideCursor;                                                      %EJD commented
    GetSecs;						% Pre-load and pre-allocate
    CharAvail;						% Pre-load and pre-allocate
    %KbName('UnifyKeyNames');
    %KbCheck(-1);
    
    %%% choose appropriate resolution
    critstimdisplay = 100;
    rwidth = 1024;	% requested resolution width
    rheight = 768;	% requested resolution height
    res1=Screen('Resolutions', 0); %list of all resolutions possible
    res2=find([res1(:).width]==rwidth & [res1(:).height]==rheight);
    res3=unique([res1(res2).hz]);
    for r=1:length(res3)
        rmod(r)=mod(critstimdisplay,1000*(1/res3(r)));  % match frame rate to value stored in critstimdisplay
    end
    [x,rmin]=min(rmod); % best frequency match
    res4=find([res1(:).width]==rwidth & [res1(:).height]==rheight & [res1(:).hz]==res3(rmin)); % best freq. & size
    [x,rmax]=max([res1(res4).pixelSize]);   % largest pixelSize value of good resolutions
    newres = res1(res4(rmax));
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
        directory = '/Users/eduwell/OneDrive - mcw.edu/Documents/SNAP/projects/Rev_Corr/code/RevCorr/working_dir/data/';
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
        %scale_f = 1;
        scale_f = 0.5;
        screenrect = screenrect * scale_f;
    end


    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'UseFastOffscreenWindows');
    window = PsychImaging('OpenWindow', 0, GrayIndex(0,0.5), screenrect);	% Open generic on-screen window
    
    Screen(window,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);    % allow for transparency!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    olddir = cd;
    fileid = strcat(out_dir_final,"/",subj);
    save(fileid,'today','link');

    day=1;  % index for determining which day of testing this is for this subj
%     if exist(strcat(subj,'1.mat'),'file')
%         %there already exists a day1 file with that Ss ID
%         load(strcat(subj,'1.mat'),'day','quest*','stair'); % load the order of conditions, day number, threshold, & stairs status for this subject
%         day = day + 1;  % advance the day counter
%         fileid1 = strcat(directory,subj,'1');
%         save(fileid1,'day','-append');  % save the new day number
%         fileid = strcat(directory,subj,num2str(day));
%         save(fileid,'today');
%     else
%         %this must be day1 for this subject
%         day = 1; %EJD UNCOMMENTED
%         stair = 'N';  %EJD UNCOMMENTED
%         fileid = strcat(directory,subj,'1');
%         save(fileid,'today','link');
%     end

    % save room location for each day
%     eval(['room' num2str(day) ' = room;']);
%     eval(['save(fileid,''room' num2str(day) ''',''-append'');']);
    
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
    
    if run_pracTrls == "y"
    % Run the Practice Session
    [trialtm_p,ts_out_p,qs_out_p,sd_out_p,nr_ideals_p, params_out_p,tdf_out_practice_TRN] = FACEquest5_int_dfTrain_practice_v2(room,kboard,quests_ntrials_prac,link,0,window,image_dir,db_mode,db_mode_screen,nr_params, imtag, pthrds,n_int,qst_startParams,t_prob,maxmin,FMpar_mat,q_prob);
    else
        tdf_out_practice_TRN = [];
    end

    % Run the Real Session
    for vv = 1:quest_testreps
                [trialtm,ts_out,qs_out,sd_out,nr_ideals, params_out,tdf_out] = FACEquest5_int_dfTrain_v2(room,kboard,quests_ntrials,link,0,window,image_dir,db_mode,db_mode_screen,nr_params, imtag, pthrds,n_int,qst_startParams,t_prob,maxmin,FMpar_mat,q_prob);
                
                
                % Unpack the outputs from the quests
                % quest estimated thresholds
                t1_1_qst = ts_out(1);
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
%                 nr_ideal_t2_qst = nr_ideals(2);
%                 nr_ideal_t3_qst = nr_ideals(3);
%                 nr_ideal_t4_qst = nr_ideals(4);
%                 nr_ideal_t5_qst = nr_ideals(5);
                
                paramq1 = params_out{1,1};
                if n_int>1
                    paramq2 = params_out{1,2};
                    if n_int>2
                        paramq3 = params_out{1,3};
                        if n_int>3
                            paramq4 = params_out{1,4};
                            if n_int>4
                                paramq5 = params_out{1,5};
                            end
                        end
                    end
                end 
                
                % Extract the raw response data from the interleaved quests
                % and fit naka rushton function to the raw data.. 
                % With this, then give the estimated psychometric thresholds 
                % from the fit NR model 
                q_list = who("q*_*_qst");
                for ii = 1:length(q_list)
                    cmd_str_qparam = strcat("q_param = ", "paramq",num2str(ii),"(1:(",q_list{ii},".trialCount));");
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

                for n = 1:n_bins
                    bin_means(:,n) = mean(q_respmat(bins==n,:))';
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
                
                if nr_opt == 0 % If nr_opt == 0, only optimize c50/fix the exponent param
                % "n" at 2...
                f = fit(x',y',"nr_rmax*((x^nr_n)/(nr_c50^nr_n + x^nr_n)) + nr_b","StartPoint",[nr_b,c50fstart, nr_n, nr_rmax],'Lower',[nr_b, c50_lower, nr_n, nr_rmax],'Upper',[nr_b, c50_upper, nr_n, nr_rmax],'Robust', 'Bisquare');
                

                end

                if nr_opt_shown2 == 1 % If nr_opt == 0, only optimize c50/fix the exponent param
                % "n" at 2...
                f2 = fit(x',y',"nr_rmax*((x^nr_n)/(nr_c50^nr_n + x^nr_n)) + nr_b","StartPoint",[nr_b,c50fstart, nr_n, nr_rmax],'Lower',[nr_b, c50_lower, nr_n, nr_rmax],'Upper',[nr_b, c50_upper, nr_n, nr_rmax],'Robust', 'Bisquare');
                end

                
                if nr_opt == 1 % If nr_opt == 1, optimize both c50 and nr_n..
                f = fit(x',y',"nr_rmax*((x^nr_n)/(nr_c50^nr_n + x^nr_n)) + nr_b","StartPoint",[nr_b,c50fstart, (nr_n_upper + nr_n_lower)/2, nr_rmax],'Lower',[nr_b, c50_lower, nr_n_lower, nr_rmax],'Upper',[nr_b, c50_upper, nr_n_upper, nr_rmax],'Robust', 'Bisquare');
                end
                
                figure()
                hold
                
                fitted_curve = double(f(interval));
                plot(x',y','o','Color','blue')
                plot(interval,fitted_curve,'Color','red')
                %plot(f,x',y')
                
                if nr_opt == 0
                fit_t1 = double(max(solve(f.nr_rmax*((c^f.nr_n)/(f.nr_c50^f.nr_n + c^f.nr_n)) + f.nr_b == pThreshold1,c)));
                end
                
                if nr_opt == 1
                    % define NR psychometric function for based on the parameters
                    % optimized in function f after the fitting above..
                    syms cf [1 1] real
                    %nr_fncn_f(cf) = nr_rmax*((cf^nr_nf)/(nr_c50f^nr_nf + cf^nr_nf)) + nr_b;
                    
                    fit_t1 = double(max(solve(f.nr_rmax*((cf^f.nr_n)/(f.nr_c50^f.nr_n + cf^f.nr_n)) + f.nr_b == pThreshold1,cf)));
                    
                end

                if db_mode == 1
                nr_ideal_t1 = double(max(solve(nr_rmax*((c^nr_n)/(nr_c50^nr_n + c^nr_n)) + nr_b == pThreshold1,c)));
                perct1_fit_err = ((abs(nr_ideal_t1-fit_t1))/nr_ideal_t1)*100; 
                end

                th_str = num2str(pThreshold1);
                plot(fit_t1,pThreshold1,'o','LineWidth',3,'Color','red')
                
                if db_mode == 1
                plot(interval,nr_verid,'blue')
                plot(nr_ideal_t1,pThreshold1,'o','LineWidth',3, 'Color','blue')
                end
                
                if nr_opt == 1 && nr_opt_shown2 == 1
                fitted_curve2 = double(f2(interval));
                plot(interval,fitted_curve2,'Color','green')
                end

                if db_mode == 1

                    if nr_opt == 0
                    legend("data","fitted curve", strcat("fitted curve"," ",th_str," threshold"), "veridical curve", strcat("veridical"," ",th_str," threshold"))
                    end
                    if nr_opt == 1
                    legend("data","fitted curve", strcat("fitted curve"," ",th_str," threshold"), "veridical curve", strcat("veridical"," ",th_str," threshold"), "fitted curve c50 opt. only..")

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
                functions_out{1,vv} = f;
                fit_t1_out{1,vv} = fit_t1;
                fit_curve_out{1,vv} = fitted_curve;
                ver_curve_out{1,vv} = nr_verid;
                

                % Print out figs if requested..
                if print_figs == 1
                    fig_fname = strcat(figs_dir,"_rep_",num2str(vv),".tif");
%                     cd(figs_dir_base)
%                     if not(isfolder(figs_dir))
%                         mkdir(figs_dir)
%                     end
%                     cd(figs_dir);
                    cd(out_dir_final)
                    imwrite(getframe(gcf).cdata, fig_fname)
                    fig_readback = imread(strcat(fig_fname));
                    cd(start_dir)
                end
                figs_out{1,vv} = fig_readback;
    end
            % Save the output parameters/data to file..
            if db_mode == 1
                save(fileid,'trialtm_out','ts_out_out','qs_out_out','sd_out_out','nr_ideals_out','nr_params','params_out_out','functions_out','fit_t1_out','ver_curve_out','figs_out','-append');
                %save(fileid,'trialtm','ts_out','qs_out','sd_out','nr_ideals','nr_params','params_out','f','fit_t1','idealNR_comparison_fig','-append');
            end
            if db_mode == 0
            save(fileid,'trialtm','ts_out','qs_out','sd_out','params_out','f','fit_t1','tdf_out','tdf_out_practice_TRN','-append');
            end
    end    
    
    % compile the images in figs_out into a .gif file for rapid viewing if requested...
    if print_fig_gif == 1
%         cd(figs_dir_base)
%         cd(figs_dir);
        cd(out_dir_final)
        nImages = size(figs_out,2);
        gif_filename = strcat(figs_dir,'_reps1-',num2str(nImages),'.gif'); % Specify the output file name
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
%         cd(figs_dir_base);
%         cd(figs_dir);
        cd(out_dir_final)
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
        %plot(interval,nr_verid,'LineWidth',2,'Color','blue')
        %alpha 0.5
        %plot(nr_ideal_t1,pThreshold1,'o','LineWidth',2, 'Color','blue')
        %alpha 0.5
        plot(interval,min_fitCrv,"--",'LineWidth',1,'Color','black')
        plot(interval,max_fitCrv,"--",'LineWidth',1,'Color','black')

        xlim([lower upper])
        ylim([0.5 1])
        fig_fname = strcat(figs_dir,"_MeanStDev",".tif");
        imwrite(getframe(gcf).cdata, fig_fname)
        %mean_sdfig = imread(strcat(fig_fname));
        cd(start_dir) 
    end
    
   
    
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
    jnk = 'pause';
end
cd(start_dir);
ListenChar(0);
ShowCursor;
Screen('CloseAll');
rethrow(lasterror);
%cd(start_dir);
%jnk = 'pause';
if serial_test == 0
clear;
end
cd(start_dir);
close all

%% Functions
function [subDirsNames] = GetSubDirsFirstLevelOnly(parentDir)
    % Get a list of all files and folders in this folder.
    files = dir(parentDir);
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subDirs = files(dirFlags);
    subDirsNames = cell(1, numel(subDirs) - 2);
    for i=3:numel(subDirs)
        subDirsNames{i-2} = subDirs(i).name;
    end
end
%end
%%%------------------------------------------------------------------------
%%%------------------------------------------------------------------------