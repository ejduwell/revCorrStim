function revCorrInstructions_v1(window,screenrect,intro_txtHeaderSize,intro_txtSize,scrn_top,scrn_bot,scrn_1percY,scrn_1percX,scrn_left,scrn_right,black,gray,whichVersion)
Screen('Preference', 'SkipSyncTests', 1); %EJD: added to skip sync tests
Screen('Preference','SuppressAllWarnings', 1);


% Detect and set directory/path params:
% -------------------------------------
% get matlab directory path by "which-ing" for this file..
% NOTE: key assumption: there is only one copy of this on the path.. (that
% should always be the case to avoid other conflicts..)
mlab_dir = fileparts(which('RevCorr_QSTmain5.m'));

cd(mlab_dir);
% cd 2 directory to enter the main dir..
cd ..
cd ..
% save the path..
main_dir = pwd;

% go back to the matlab dir..
cd(mlab_dir)
% -------------------------------------

%% Generate the example images
% whichVersion="tex"; % 'lum', 'tex', or 'cr'

if whichVersion=="lum"

    % for lum
    parzIn=struct; % initialize
    parzIn.occImgR=imread(strcat(main_dir,"/images/test3_LumOnly-26-Feb-2024-19-00-35/occ/R/BaseIm_occ_0_ori_m__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3214_T_0clovers_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_0_lum2_255_LWT_1_Lal_1_Ocon_1.png"));
    parzIn.occImgL=imread(strcat(main_dir,"/images/test3_LumOnly-26-Feb-2024-19-00-35/occ/L/BaseIm_occ_0_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3214_T_0clovers_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_0_lum2_255_LWT_1_Lal_1_Ocon_1.png"));
    parzIn.noccImgR=imread(strcat(main_dir,"/images/test3_LumOnly-26-Feb-2024-19-00-35/nocc/R/BaseIm_occ_1_ori_m__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3214_T_0clovers_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_0_lum2_255_LWT_1_Lal_1_Ocon_1.png"));
    parzIn.noccImgL=imread(strcat(main_dir,"/images/test3_LumOnly-26-Feb-2024-19-00-35/nocc/L/BaseIm_occ_1_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3214_T_0clovers_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_0_lum2_255_LWT_1_Lal_1_Ocon_1.png"));
    parzIn.noiseDir=strcat(main_dir,"/noise/lumOnlyBI_20000frms_krnlNz_imzPerKrnl_111111100");
    parzIn.fix_im=imread(strcat(main_dir,"/images/test3_LumOnly-26-Feb-2024-19-00-35/fixation.png"));
    parzIn.nWt=0.5;

end

if whichVersion=="cr"

    % for cr
    parzIn=struct; % initialize
    parzIn.occImgR=imread(strcat(main_dir,"/images/test3_CROnlyWlum51-26-Feb-2024-18-52-21/occ/R/BaseIm_occ_0_ori_m__CR_1_CRpw_4_CRpl_0_CRal_1_CRob_1_CRdm_3214_T_0clovers_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_51_lum2_51_LWT_1_Lal_1_Ocon_1.png"));
    parzIn.occImgL=imread(strcat(main_dir,"/images/test3_CROnlyWlum51-26-Feb-2024-18-52-21/occ/L/BaseIm_occ_0_ori_z__CR_1_CRpw_4_CRpl_0_CRal_1_CRob_1_CRdm_3214_T_0clovers_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_51_lum2_51_LWT_1_Lal_1_Ocon_1.png"));
    parzIn.noccImgR=imread(strcat(main_dir,"/images/test3_CROnlyWlum51-26-Feb-2024-18-52-21/nocc/R/BaseIm_occ_1_ori_m__CR_1_CRpw_4_CRpl_0_CRal_1_CRob_1_CRdm_3214_T_0clovers_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_51_lum2_51_LWT_1_Lal_1_Ocon_1.png"));
    parzIn.noccImgL=imread(strcat(main_dir,"/images/test3_CROnlyWlum51-26-Feb-2024-18-52-21/nocc/L/BaseIm_occ_1_ori_z__CR_1_CRpw_4_CRpl_0_CRal_1_CRob_1_CRdm_3214_T_0clovers_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_51_lum2_51_LWT_1_Lal_1_Ocon_1.png"));
    parzIn.noiseDir=strcat(main_dir,"/noise/crOnlyBI_20000frms_krnlNz_imzPerKrnl_111111100");
    parzIn.fix_im=imread(strcat(main_dir,"/images/test3_CROnlyWlum51-26-Feb-2024-18-52-21/fixation.png"));
    parzIn.nWt=0.5;

end

if whichVersion=="tex"

    % for tex
    parzIn=struct; % initialize
    parzIn.occImgR=imread(strcat(main_dir,"images/test5_TexOnly-12-Mar-2024-18-48-17/occ/R/BaseIm_occ_0_ori_m__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3315_T_1checkrz_Tr1_1.3_Tr2_0.76923_Tal_1_L_0_lum1_127_lum2_127_LWT_1_Lal_1_Ocon_1.png"));
    parzIn.occImgL=imread(strcat(main_dir,"/images/test5_TexOnly-12-Mar-2024-18-48-17/occ/L/BaseIm_occ_0_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3315_T_1checkrz_Tr1_1.3_Tr2_0.76923_Tal_1_L_0_lum1_127_lum2_127_LWT_1_Lal_1_Ocon_1.png"));
    parzIn.noccImgR=imread(strcat(main_dir,"/images/test5_TexOnly-12-Mar-2024-18-48-17/nocc/R/BaseIm_occ_1_ori_m__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3315_T_1checkrz_Tr1_1.3_Tr2_0.76923_Tal_1_L_0_lum1_127_lum2_127_LWT_1_Lal_1_Ocon_1.png"));
    parzIn.noccImgL=imread(strcat(main_dir,"/images/test5_TexOnly-12-Mar-2024-18-48-17/nocc/L/BaseIm_occ_1_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3315_T_1checkrz_Tr1_1.3_Tr2_0.76923_Tal_1_L_0_lum1_127_lum2_127_LWT_1_Lal_1_Ocon_1.png"));
    parzIn.noiseDir=strcat(main_dir,"/noise/texOnly4BI_krnlNz_imzPerKrnl_111111100");
    parzIn.fix_im=imread(strcat(main_dir,"/images/test5_TexOnly-12-Mar-2024-18-48-17/fixation.png"));
    parzIn.nWt=0.5;

end

% Without Noise
parzIn.useNoise=0;
exmplImgs = revCorrMakeExmplImgs(parzIn);
% With Noise
parzIn.useNoise=1;
exmplImgs_wNoise = revCorrMakeExmplImgs(parzIn);
% Combination
comboImg=horzcat(exmplImgs,exmplImgs_wNoise);

%% Build instructions page
    %Start KbQueue
    KbQueueCreate();

    % INSTRUCTIONS PAGE 1: 
    % ---------------------------------------------------------------------
    % A) Create Page
    Screen('FillRect', window, gray); % gray out window
    
    speech1 = ['Instructions:'];
    Screen('TextSize',window, intro_txtHeaderSize);
    headerDistFrmTop=(4*scrn_1percY); % specify distance of header from top as percent of screen dimension..
    DrawFormattedText(window, speech1,'center', (scrn_top+headerDistFrmTop),black); % orig
    
    txtDistFrmTop1=(8*scrn_1percY); % specify distance of text from top as percent of screen dimension..
    speech2_1_1 = ['In this experiment, you''ll see sets of two rectangular objects with noise overlaid. \n' ...
        'You''ll also see a central fixation square. You must always look at the central square. \n' ...
        'On some trials the rectangles will be fully visible. \n' ...
        'On other trials, portions of the rectangles will be hidden behind occluder bars such that only parts of the objects are visible:'];
   

    scrn_Ctr=(scrn_left+scrn_right)/2;
    imLft=scrn_Ctr-(size(comboImg,2)/2);
    imRght=scrn_Ctr+(size(comboImg,2)/2);
    imTop=(40*scrn_1percY);
    imBot=(40*scrn_1percY)+size(comboImg,1);
    imgRect=[imTop,imLft,imBot,imRght];
    
    % draw the images
    comboImg = imresize(comboImg,0.85,"nearest");
    Screen('PutImage', window, comboImg);
    %Screen('PutImage', window, comboImg,imgRect);

    txtDistFrmTop2=85*scrn_1percY;
    speech2_1_2 = ['Your task will ALWAYS be to identify whether the objects are angled RIGHT or LEFT. \n' ...
        'If they''re angled RIGHT, press the ''m'' key on the keyboard, as quickly and as accurately as possible. \n' ... 
        'If they''re angled LEFT, press the ''z'' key on the keyboard, as quickly and as accurately as possible.'];
        
    txtDistFrmTop3=95*scrn_1percY; % specify distance of text from top as percent of screen dimension..
    speech2_2=['Press RETURN/ENTER to continue.'];
    
    Screen('TextSize',window, intro_txtSize);
    DrawFormattedText(window, speech2_1_1,'center', (scrn_top+txtDistFrmTop1),black);
    DrawFormattedText(window, speech2_1_2,'center', (scrn_top+txtDistFrmTop2),black);
    DrawFormattedText(window, speech2_2,'center', (scrn_top+txtDistFrmTop3),black);
   
    % B) Display page until button is pressed
    keyIsDown = 0;
    enter_pressed = 0;
    swtch = 0;
    respns = "";
    %t_trial0 = GetSecs;
    % While stimulus is on the screen and no response
    while enter_pressed == 0 %trialResp == 0  

        % Draw stimulus
        if swtch == 0
            Screen('Flip', window);% show linking display
            t_disp = GetSecs; % record time display came up..
            %start listening for responses
            KbQueueStart();
            swtch = 1;
        end

        [keyIsDown, firstPress, ~, ~, ~] = KbQueueCheck();
        if keyIsDown == 1
            if KbName(firstPress) == 'Return'
                enter_pressed = 1;
            end
        end
    end
    enter_pressed = 0;
    keyIsDown = 0;
    swtch = 0;
    % ---------------------------------------------------------------------
        
    % INSTRUCTIONS PAGE 2 Example Image..: 
    % ---------------------------------------------------------------------
    % A) Create Page
    Screen('FillRect', window, gray); % gray out window
    
    % add formatted text..
    txtDistFrmTop4=50*scrn_1percY; % specify distance of text from top as percent of screen dimension..
    Screen('TextSize',window, intro_txtSize);
    speech1 = ['Once you respond (or after X seconds, if you don''t respond), the experiment will go on to the next trial. \n' ...
        'Try to identify the target as quickly and accurately as possible. \n\n' ...
        'Press ENTER/RETURN to Continue...'];
    DrawFormattedText(window, speech1,'center', (scrn_top+txtDistFrmTop4), black);
    
    % B) Display page until button is pressed
    keyIsDown = 0;
    enter_pressed = 0;
    swtch = 0;
    respns = "";
    %t_trial0 = GetSecs;
    % While stimulus is on the screen and no response
    while enter_pressed == 0 %trialResp == 0  

        % Draw stimulus
        if swtch == 0
            Screen('Flip', window);% show linking display
            t_disp = GetSecs; % record time display came up..
            %start listening for responses
            KbQueueStart();
            swtch = 1;
        end

        [keyIsDown, firstPress, ~, ~, ~] = KbQueueCheck();
        if keyIsDown == 1
            if KbName(firstPress) == 'Return'
                enter_pressed = 1;
            end
        end
    end
    enter_pressed = 0;
    keyIsDown = 0;
    swtch = 0;
    % ---------------------------------------------------------------------
% 
%      % INSTRUCTIONS PAGE 3: (Task instructions)
%     % ---------------------------------------------------------------------
%     % A) Create Page
%     Screen('FillRect', window, gray); % gray out window
%     % add formatted text..
%     speech1 = ['INSERT TASK INSTRUCTIONS HERE... \n\n' ...
%         ' When you''re ready to begin, please\n\n' ...
%         'press RETURN/ENTER.'];
%     Screen('TextSize',window, intro_txtSize);
%     DrawFormattedText(window, speech1,'center','center',black);
% 
% 
%     % B) Display page until button is pressed
%     keyIsDown = 0;
%     enter_pressed = 0;
%     swtch = 0;
%     respns = "";
%     %t_trial0 = GetSecs;
%     % While stimulus is on the screen and no response
%     while enter_pressed == 0 
% 
%         % Draw stimulus
%         if swtch == 0
%             Screen('Flip', window);% show linking display
%             t_disp = GetSecs; % record time display came up..
%             %start listening for responses
%             KbQueueStart();
%             swtch = 1;
%         end
% 
%         [keyIsDown, firstPress, ~, ~, ~] = KbQueueCheck();
%         if keyIsDown == 1
%             if KbName(firstPress) == 'Return'
%                 enter_pressed = 1;
%             end
%         end
%     end
%     enter_pressed = 0;
%     keyIsDown = 0;
%     swtch = 0;


    Screen('FillRect', window, gray); % gray out window

end