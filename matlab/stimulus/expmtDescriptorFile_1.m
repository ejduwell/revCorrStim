function expParz = expmtDescriptorFile_1()

% This is an instance of a expmtDescriptorFile intended to be used for
% specifying key parameters used by reverse correlation experiments 
% run using the revCorrStim repository. 
%
% The expmtDescriptorFile is essentially a 'master descriptor file' that
% is intended to serve as the 'one-stop-shop' where key parameters
% regarding choice of base images, noise images, QUEST parameters, and trial
% timing parameters are set.
%
% parameters set here should all be stored within a unified struct variable
% called 'expParz'.
% 
% This struct is exported as the output variable such that other programs
% in the revCorrStim repository can access them. 

%% Initialize Output Struct
expParz=struct;

%% Set Base Image Parameters

% General to all
expParz.BI_parz.gen.imtag = ".png"; % unique tag in filename iding the base images you want...

% LUMINANCE ONLY VERSION:
% #########################################################################

% luminance bi directory:
expParz.BI_dirs.lum="test4_LumOnly-05-Mar-2024-10-12-48";

expParz.BI_parz.lum.qstParStr="lumDiff"; % string tag indicating quest should manipulate: diff between lum1/lum2
expParz.BI_parz.lum.parTagz = ["lum1","lum2"]; % file name tags for parameters manipulated
expParz.BI_parz.lum.obj1Par="lum1"; % object 1 parameter
expParz.BI_parz.lum.obj2Par="lum2"; % object 2 parameter
expParz.BI_parz.lum.nWt=0.75; % weighting of noise vs. bi
% #########################################################################

% TEXTURE ONLY VERSION:
% #########################################################################

% texture bi directory:
expParz.BI_dirs.tex="test6_TexOnly-19-Apr-2024-11-43-14";

expParz.BI_parz.tex.qstParStr="texDiff"; % string tag indicating quest should manipulate: diff between texture resolution1/texture resolution2
expParz.BI_parz.tex.parTagz = ["Tr1","Tr2"]; % file name tags for parameters manipulated
expParz.BI_parz.tex.obj1Par="Tr1"; % object 1 parameter
expParz.BI_parz.tex.obj2Par="Tr2"; % object 2 parameter
expParz.BI_parz.tex.nWt=0.75; % weighting of noise vs. bi
% #########################################################################

% COMMON REGION ONLY VERSION:
% #########################################################################

% common region bi directory:
expParz.BI_dirs.cr="test5_CROnly-11-Mar-2024-14-16-09";

expParz.BI_parz.cr.qstParStr="CRpl"; % string tag indicating quest should manipulate: common region boundary luminance
expParz.BI_parz.cr.parTagz = ["CRpl"]; % file name tags for parameters manipulated
expParz.BI_parz.cr.obj1Par="CRpl"; % object 1 parameter
expParz.BI_parz.cr.obj2Par="CRpl"; % object 2 parameter
expParz.BI_parz.cr.nWt=0.75; % weighting of noise vs. bi
% #########################################################################

%% Set Noise Image Parameters

% LUMINANCE ONLY VERSION:
% #########################################################################

% White Noise
expParz.noise.lum.white.noise_dir="512by512_whiteNoise_20000frms_smpl1";
expParz.noise.lum.white.nTag=".png";

% Kernel Noise
expParz.noise.lum.kernel.noise_dir="lumOnlyBI_20000frms_krnlNz_imzPerKrnl_111111100";
expParz.noise.lum.kernel.nTag=".png";

% #########################################################################

% TEXTURE ONLY VERSION:
% #########################################################################

% White Noise
expParz.noise.tex.white.noise_dir="512by512_whiteNoise_20000frms_smpl2";
expParz.noise.tex.white.nTag=".png";

% Kernel Noise
expParz.noise.tex.kernel.noise_dir="texOnlyBI_v2_20000frms_krnlNz_imzPerKrnl_111111100";
expParz.noise.tex.kernel.nTag=".png";

% #########################################################################

% COMMON REGION ONLY VERSION:
% #########################################################################

% White Noise
expParz.noise.cr.white.noise_dir="512by512_whiteNoise_20000frms_smpl3";
expParz.noise.cr.white.nTag=".png";

% Kernel Noise
expParz.noise.cr.kernel.noise_dir="crOnlyBI_20000frms_krnlNz_imzPerKrnl_111111100";
expParz.noise.cr.kernel.nTag=".png";

% #########################################################################

%% Set QUEST Parameters

% General to all:
% #########################################################################
expParz.questParz.n_int = 3; % specify the number of interleved quests you want.. if 1,
                             % no interleaving will occur.. only threshold 1 will be used. Can range
                             % from 0-5.. 5 is the max # of quests you can currently interleave ..

expParz.questParz.quests_ntrials_init = 50; % specify number of trials (per quest) in each quest test rep on the first quest block
expParz.questParz.quests_ntrials_reg = 100; % specify number of trials (per quest) in each quest test rep  on the quest block after the first block
expParz.questParz.quests_ntrials_practice = 6; % specify number of trials (per quest) in each quest test rep

expParz.questParz.t_prob = 0.5; % probablility/proportion of right vs. left trials (FOR THE ACTUAL EXPERIMENT)
expParz.questParz.t_prob_practice = 0.5; % probablility/proportion of right vs. left trials (FOR THE PRACTICE SESSION)

% Interleaved staircase threshold parameters..
% NOTE: DEPENDING ON THE # OF INTERLEAVED QUESTS SELECETED (WITH N_INT
% ABOVE) SOME OF THESE THRESHOLDS MAY OR MAY NOT BE USED.. IT USES THE
% PTHRESHOLDS 1-N_INT...
%Staircase1
expParz.questParz.pThreshold1=0.65; % NOTE: this is the one compared on the graphs..
%Staircase2
expParz.questParz.pThreshold2=0.75;
%Staircase3
expParz.questParz.pThreshold3=0.85;
%Staircase4
expParz.questParz.pThreshold4=0.85;
%Staircase5
expParz.questParz.pThreshold5=0.85;
% #########################################################################

% Set different max/min parameter limits based on experiment type...
% #########################################################################

% LUMINANCE:
expParz.questParz.lum.maxmin = cell(1,2); % initialize maxmin
expParz.questParz.lum.maxmin{1,1} = 255; % specifies the maximum parameter value allowed in quest..
expParz.questParz.lum.maxmin{1,2} = 0; % specifies the minimum parameter value allowed in quest..
expParz.questParz.lum.maxmin_seg = 0; % if 1, will segment the curve into n_int pieces and set
% separate maxs/mins for each segment to try to ensure even coverage..
% NOTE: IF YOU SELECT MAXMIN_SEG, MAKE SURE YOUR PTHRESHOLDS ARE IN
% ASCENDING ORDER!!!


% TEXTURE:
expParz.questParz.tex.maxmin = cell(1,2); % initialize maxmin
expParz.questParz.tex.maxmin{1,1} = 6.25; % specifies the maximum parameter value allowed in quest..
expParz.questParz.tex.maxmin{1,2} = 0; % specifies the minimum parameter value allowed in quest..
expParz.questParz.tex.maxmin_seg = 0; % if 1, will segment the curve into n_int pieces and set
% separate maxs/mins for each segment to try to ensure even coverage..
% NOTE: IF YOU SELECT MAXMIN_SEG, MAKE SURE YOUR PTHRESHOLDS ARE IN
% ASCENDING ORDER!!!

% COMMON REGION:
expParz.questParz.cr.maxmin = cell(1,2); % initialize maxmin
expParz.questParz.cr.maxmin{1,1} = 127.5; % specifies the maximum parameter value allowed in quest..
expParz.questParz.cr.maxmin{1,2} = 0; % specifies the minimum parameter value allowed in quest..
expParz.questParz.cr.maxmin_seg = 0; % if 1, will segment the curve into n_int pieces and set
%separate maxs/mins for each segment to try to ensure even coverage..
% NOTE: IF YOU SELECT MAXMIN_SEG, MAKE SURE YOUR PTHRESHOLDS ARE IN
% ASCENDING ORDER!!!
% #########################################################################

%% Set Trial Parameters

expParz.trialParz.presTime=0.25; % time after trial start that stimuli are presented ..(sec)
expParz.trialParz.presTime2=expParz.trialParz.presTime+0.5; % time after trial start that stimuli are removed ..(sec)
expParz.trialParz.resp_cut = 2; % Max time after trial start that responses are recorded
expParz.trialParz.ITI=0.5; % inter-trial interval time (sec).. (time btw trials)



%% Set Parameters for Example Images / Instructions Page

% for lum
expParz.exImz.lum.bi.occImgR="/images/test4_LumOnly-05-Mar-2024-10-12-48/occ/R/BaseIm_occ_0_ori_m__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3214_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_0_lum2_255_LWT_1_Lal_1_Ocon_1.png";
expParz.exImz.lum.bi.occImgL="/images/test4_LumOnly-05-Mar-2024-10-12-48/occ/L/BaseIm_occ_0_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3214_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_0_lum2_255_LWT_1_Lal_1_Ocon_1.png";
expParz.exImz.lum.bi.noccImgR="/images/test4_LumOnly-05-Mar-2024-10-12-48/nocc/R/BaseIm_occ_1_ori_m__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3214_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_0_lum2_255_LWT_1_Lal_1_Ocon_1.png";
expParz.exImz.lum.bi.noccImgL="/images/test4_LumOnly-05-Mar-2024-10-12-48/nocc/L/BaseIm_occ_1_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3214_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_0_lum2_255_LWT_1_Lal_1_Ocon_1.png";
expParz.exImz.lum.bi.fix_im="/images/test4_LumOnly-05-Mar-2024-10-12-48/fixation.png";
expParz.exImz.lum.noiseDir.kernel="/noise/lumOnlyBI_20000frms_krnlNz_imzPerKrnl_111111100";
expParz.exImz.lum.noiseDir.white="/noise/512by512_whiteNoise_20000frms_smpl1";
expParz.exImz.lum.nWt=0.5;


% for cr
expParz.exImz.cr.bi.occImgR="/images/test5_CROnly-11-Mar-2024-14-16-09/occ/R/BaseIm_occ_0_ori_m__CR_1_CRpw_4_CRpl_0_CRal_1_CRob_1_CRdm_3315_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_51_lum2_51_LWT_1_Lal_1_Ocon_1.png";
expParz.exImz.cr.bi.occImgL="/images/test5_CROnly-11-Mar-2024-14-16-09/occ/L/BaseIm_occ_0_ori_z__CR_1_CRpw_4_CRpl_0_CRal_1_CRob_1_CRdm_3315_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_51_lum2_51_LWT_1_Lal_1_Ocon_1.png";
expParz.exImz.cr.bi.noccImgR="/images/test5_CROnly-11-Mar-2024-14-16-09/nocc/R/BaseIm_occ_1_ori_m__CR_1_CRpw_4_CRpl_0_CRal_1_CRob_1_CRdm_3315_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_51_lum2_51_LWT_1_Lal_1_Ocon_1.png";
expParz.exImz.cr.bi.noccImgL="/images/test5_CROnly-11-Mar-2024-14-16-09/nocc/L/BaseIm_occ_1_ori_z__CR_1_CRpw_4_CRpl_0_CRal_1_CRob_1_CRdm_3315_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_51_lum2_51_LWT_1_Lal_1_Ocon_1.png";
expParz.exImz.cr.bi.fix_im="/images/test5_CROnly-11-Mar-2024-14-16-09/fixation.png";
expParz.exImz.cr.noiseDir.kernel="/noise/crOnlyBI_20000frms_krnlNz_imzPerKrnl_111111100";
expParz.exImz.cr.noiseDir.white="/noise/512by512_whiteNoise_20000frms_smpl2";
expParz.exImz.cr.nWt=0.5;


% for tex
expParz.exImz.tex.bi.occImgR="/images/test6_TexOnly-19-Apr-2024-11-43-14/occ/R/BaseIm_occ_0_ori_m__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3315_T_1checkrz_Tr1_2.5_Tr2_0.4_Tal_1_L_0_lum1_127_lum2_127_LWT_1_Lal_1_Ocon_1.png";
expParz.exImz.tex.bi.occImgL="/images/test6_TexOnly-19-Apr-2024-11-43-14/occ/L/BaseIm_occ_0_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3315_T_1checkrz_Tr1_2.5_Tr2_0.4_Tal_1_L_0_lum1_127_lum2_127_LWT_1_Lal_1_Ocon_1.png";
expParz.exImz.tex.bi.noccImgR="/images/test6_TexOnly-19-Apr-2024-11-43-14/nocc/R/BaseIm_occ_1_ori_m__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3315_T_1checkrz_Tr1_2.5_Tr2_0.4_Tal_1_L_0_lum1_127_lum2_127_LWT_1_Lal_1_Ocon_1.png";
expParz.exImz.tex.bi.noccImgL="/images/test6_TexOnly-19-Apr-2024-11-43-14/nocc/L/BaseIm_occ_1_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3315_T_1checkrz_Tr1_2.5_Tr2_0.4_Tal_1_L_0_lum1_127_lum2_127_LWT_1_Lal_1_Ocon_1.png";
expParz.exImz.tex.bi.fix_im="/images/test6_TexOnly-19-Apr-2024-11-43-14/fixation.png";
expParz.exImz.tex.noiseDir.kernel="/noise/texOnlyBI_v2_20000frms_krnlNz_imzPerKrnl_111111100";
expParz.exImz.tex.noiseDir.white="/noise/512by512_whiteNoise_20000frms_smpl3";
expParz.exImz.tex.nWt=0.5;
