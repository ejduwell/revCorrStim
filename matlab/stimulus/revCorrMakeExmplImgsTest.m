%% Set Parameters..

whichVersion="tex"; % 'lum', 'tex', or 'cr'

if whichVersion=="lum"

    % for lum
    parzIn=struct; % initialize
    parzIn.occImgR=imread("/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/images/test3_LumOnly-26-Feb-2024-19-00-35/occ/R/BaseIm_occ_0_ori_m__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3214_T_0clovers_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_0_lum2_255_LWT_1_Lal_1_Ocon_1.png");
    parzIn.occImgL=imread("/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/images/test3_LumOnly-26-Feb-2024-19-00-35/occ/L/BaseIm_occ_0_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3214_T_0clovers_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_0_lum2_255_LWT_1_Lal_1_Ocon_1.png");
    parzIn.noccImgR=imread("/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/images/test3_LumOnly-26-Feb-2024-19-00-35/nocc/R/BaseIm_occ_1_ori_m__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3214_T_0clovers_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_0_lum2_255_LWT_1_Lal_1_Ocon_1.png");
    parzIn.noccImgL=imread("/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/images/test3_LumOnly-26-Feb-2024-19-00-35/nocc/L/BaseIm_occ_1_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3214_T_0clovers_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_0_lum2_255_LWT_1_Lal_1_Ocon_1.png");
    parzIn.noiseDir="/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/noise/lumOnlyBI_20000frms_krnlNz_imzPerKrnl_111111100";
    parzIn.fix_im=imread("/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/images/test3_LumOnly-26-Feb-2024-19-00-35/fixation.png");
    parzIn.nWt=0.5;

end

if whichVersion=="cr"

    % for cr
    parzIn=struct; % initialize
    parzIn.occImgR=imread("/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/images/test3_CROnlyWlum51-26-Feb-2024-18-52-21/occ/R/BaseIm_occ_0_ori_m__CR_1_CRpw_4_CRpl_0_CRal_1_CRob_1_CRdm_3214_T_0clovers_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_51_lum2_51_LWT_1_Lal_1_Ocon_1.png");
    parzIn.occImgL=imread("/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/images/test3_CROnlyWlum51-26-Feb-2024-18-52-21/occ/L/BaseIm_occ_0_ori_z__CR_1_CRpw_4_CRpl_0_CRal_1_CRob_1_CRdm_3214_T_0clovers_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_51_lum2_51_LWT_1_Lal_1_Ocon_1.png");
    parzIn.noccImgR=imread("/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/images/test3_CROnlyWlum51-26-Feb-2024-18-52-21/nocc/R/BaseIm_occ_1_ori_m__CR_1_CRpw_4_CRpl_0_CRal_1_CRob_1_CRdm_3214_T_0clovers_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_51_lum2_51_LWT_1_Lal_1_Ocon_1.png");
    parzIn.noccImgL=imread("/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/images/test3_CROnlyWlum51-26-Feb-2024-18-52-21/nocc/L/BaseIm_occ_1_ori_z__CR_1_CRpw_4_CRpl_0_CRal_1_CRob_1_CRdm_3214_T_0clovers_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_51_lum2_51_LWT_1_Lal_1_Ocon_1.png");
    parzIn.noiseDir="/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/noise/crOnlyBI_20000frms_krnlNz_imzPerKrnl_111111100";
    parzIn.fix_im=imread("/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/images/test3_CROnlyWlum51-26-Feb-2024-18-52-21/fixation.png");
    parzIn.nWt=0.5;

end

if whichVersion=="tex"

    % for tex
    parzIn=struct; % initialize
    parzIn.occImgR=imread("/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/images/test5_TexOnly-12-Mar-2024-18-48-17/occ/R/BaseIm_occ_0_ori_m__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3315_T_1checkrz_Tr1_1.3_Tr2_0.76923_Tal_1_L_0_lum1_127_lum2_127_LWT_1_Lal_1_Ocon_1.png");
    parzIn.occImgL=imread("/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/images/test5_TexOnly-12-Mar-2024-18-48-17/occ/L/BaseIm_occ_0_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3315_T_1checkrz_Tr1_1.3_Tr2_0.76923_Tal_1_L_0_lum1_127_lum2_127_LWT_1_Lal_1_Ocon_1.png");
    parzIn.noccImgR=imread("/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/images/test5_TexOnly-12-Mar-2024-18-48-17/nocc/R/BaseIm_occ_1_ori_m__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3315_T_1checkrz_Tr1_1.3_Tr2_0.76923_Tal_1_L_0_lum1_127_lum2_127_LWT_1_Lal_1_Ocon_1.png");
    parzIn.noccImgL=imread("/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/images/test5_TexOnly-12-Mar-2024-18-48-17/nocc/L/BaseIm_occ_1_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3315_T_1checkrz_Tr1_1.3_Tr2_0.76923_Tal_1_L_0_lum1_127_lum2_127_LWT_1_Lal_1_Ocon_1.png");
    parzIn.noiseDir="/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/noise/texOnly4BI_krnlNz_imzPerKrnl_111111100";
    parzIn.fix_im=imread("/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/images/test5_TexOnly-12-Mar-2024-18-48-17/fixation.png");
    parzIn.nWt=0.5;

end

%% Test Call

% Without Noise
parzIn.useNoise=0;
exmplImgs = revCorrMakeExmplImgs(parzIn);

% With Noise
parzIn.useNoise=1;
exmplImgs_wNoise = revCorrMakeExmplImgs(parzIn);

comboImg=horzcat(exmplImgs,exmplImgs_wNoise);

figure;
imshow(comboImg);
axis image;