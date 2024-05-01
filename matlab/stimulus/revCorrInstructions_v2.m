function revCorrInstructions_v2(window,screenrect,intro_txtHeaderSize,intro_txtSize,scrn_top,scrn_bot,scrn_1percY,scrn_1percX,scrn_left,scrn_right,black,gray,whichVersion,noizType,exImz)
Screen('Preference', 'SkipSyncTests', 1); %EJD: added to skip sync tests
Screen('Preference','SuppressAllWarnings', 1);


% Detect and set directory/path params:
% -------------------------------------
% get matlab directory path by "which-ing" for this file..
% NOTE: key assumption: there is only one copy of this on the path.. (that
% should always be the case to avoid other conflicts..)
mlab_dir = fileparts(which('RevCorr_QSTmain7.m'));

cd(mlab_dir);
% cd 2 directory to enter the main dir..
cd ..
cd ..
% save the path..
main_dir = pwd;

% go back to the matlab dir..
cd(mlab_dir)
% -------------------------------------

printImz=1; %if 1, will save images of instruction pages..
imOutDir=strcat(main_dir,"/matlab/stimulus/instructions");
%% Generate the example images
% whichVersion="tex"; % 'lum', 'tex', or 'cr'

if noizType=="krnlNz"
if whichVersion=="lum"

    % for lum
    parzIn=struct; % initialize
    parzIn.occImgR=imread(strcat(main_dir,exImz.lum.bi.occImgR));
    parzIn.occImgL=imread(strcat(main_dir,exImz.lum.bi.occImgL));
    parzIn.noccImgR=imread(strcat(main_dir,exImz.lum.bi.noccImgR));
    parzIn.noccImgL=imread(strcat(main_dir,exImz.lum.bi.noccImgL));
    parzIn.fix_im=imread(strcat(main_dir,exImz.lum.bi.fix_im));
    parzIn.noiseDir=strcat(main_dir,exImz.lum.noiseDir.kernel);
    parzIn.nWt=exImz.lum.nWt;

end

if whichVersion=="cr"

    parzIn=struct; % initialize
    parzIn.occImgR=imread(strcat(main_dir,exImz.cr.bi.occImgR));
    parzIn.occImgL=imread(strcat(main_dir,exImz.cr.bi.occImgL));
    parzIn.noccImgR=imread(strcat(main_dir,exImz.cr.bi.noccImgR));
    parzIn.noccImgL=imread(strcat(main_dir,exImz.cr.bi.noccImgL));
    parzIn.fix_im=imread(strcat(main_dir,exImz.cr.bi.fix_im));
    parzIn.noiseDir=strcat(main_dir,exImz.cr.noiseDir.kernel);
    parzIn.nWt=exImz.cr.nWt;

end

if whichVersion=="tex"

    parzIn=struct; % initialize
    parzIn.occImgR=imread(strcat(main_dir,exImz.tex.bi.occImgR));
    parzIn.occImgL=imread(strcat(main_dir,exImz.tex.bi.occImgL));
    parzIn.noccImgR=imread(strcat(main_dir,exImz.tex.bi.noccImgR));
    parzIn.noccImgL=imread(strcat(main_dir,exImz.tex.bi.noccImgL));
    parzIn.fix_im=imread(strcat(main_dir,exImz.tex.bi.fix_im));
    parzIn.noiseDir=strcat(main_dir,exImz.tex.noiseDir.kernel);
    parzIn.nWt=exImz.tex.nWt;

end
end

if noizType=="white"
if whichVersion=="lum"

    % for lum
    parzIn=struct; % initialize
    parzIn.occImgR=imread(strcat(main_dir,exImz.lum.bi.occImgR));
    parzIn.occImgL=imread(strcat(main_dir,exImz.lum.bi.occImgL));
    parzIn.noccImgR=imread(strcat(main_dir,exImz.lum.bi.noccImgR));
    parzIn.noccImgL=imread(strcat(main_dir,exImz.lum.bi.noccImgL));
    parzIn.fix_im=imread(strcat(main_dir,exImz.lum.bi.fix_im));
    parzIn.noiseDir=strcat(main_dir,exImz.lum.noiseDir.white);
    parzIn.nWt=exImz.lum.nWt;

end

if whichVersion=="cr"

    parzIn=struct; % initialize
    parzIn.occImgR=imread(strcat(main_dir,exImz.cr.bi.occImgR));
    parzIn.occImgL=imread(strcat(main_dir,exImz.cr.bi.occImgL));
    parzIn.noccImgR=imread(strcat(main_dir,exImz.cr.bi.noccImgR));
    parzIn.noccImgL=imread(strcat(main_dir,exImz.cr.bi.noccImgL));
    parzIn.fix_im=imread(strcat(main_dir,exImz.cr.bi.fix_im));
    parzIn.noiseDir=strcat(main_dir,exImz.cr.noiseDir.white);
    parzIn.nWt=exImz.cr.nWt;

end

if whichVersion=="tex"

    parzIn=struct; % initialize
    parzIn.occImgR=imread(strcat(main_dir,exImz.tex.bi.occImgR));
    parzIn.occImgL=imread(strcat(main_dir,exImz.tex.bi.occImgL));
    parzIn.noccImgR=imread(strcat(main_dir,exImz.tex.bi.noccImgR));
    parzIn.noccImgL=imread(strcat(main_dir,exImz.tex.bi.noccImgL));
    parzIn.fix_im=imread(strcat(main_dir,exImz.tex.bi.fix_im));
    parzIn.noiseDir=strcat(main_dir,exImz.tex.noiseDir.white);
    parzIn.nWt=exImz.tex.nWt;

end
end


% Without Noise
parzIn.useNoise=0;
exmplImgs = revCorrMakeExmplImgs(parzIn);
% With Noise
parzIn.useNoise=1;
exmplImgs_wNoise = revCorrMakeExmplImgs(parzIn);
% Combination
comboImg=horzcat(exmplImgs,exmplImgs_wNoise);

if printImz==1
if noizType=="krnlNz"
    % Capture window as image..
    out_name1 = strcat(imOutDir,"/","instrxns_","exmplImz","_",whichVersion,"_wKrnlNoise.png");
    out_name2 = strcat(imOutDir,"/","instrxns_","exmplImz","_",whichVersion,".png");
    imwrite(exmplImgs_wNoise,out_name1,"png");
    imwrite(exmplImgs,out_name2,"png");
end
if noizType=="white"
    % Capture window as image..
    out_name1 = strcat(imOutDir,"/","instrxns_","exmplImz","_",whichVersion,"_wWhiteNoise.png");
    out_name2 = strcat(imOutDir,"/","instrxns_","exmplImz","_",whichVersion,".png");
    imwrite(exmplImgs_wNoise,out_name1,"png");
    imwrite(exmplImgs,out_name2,"png");
end
end


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
        'On other trials, portions of the rectangles will be hidden behind occluder bars \n'...
        'such that only parts of the objects are visible:'];
   

    scrn_Ctr=(scrn_left+scrn_right)/2;
    imLft=scrn_Ctr-(size(comboImg,2)/2);
    imRght=scrn_Ctr+(size(comboImg,2)/2);
    imTop=(40*scrn_1percY);
    imBot=(40*scrn_1percY)+size(comboImg,1);
    imgRect=[imTop,imLft,imBot,imRght];
    
    % draw the images
    sclFctr=0.6/14.4;
    ImgSclFctr=sclFctr*scrn_1percY;
    exmplImgs_wNoise=imresize(exmplImgs_wNoise,ImgSclFctr,"nearest");

    Screen('PutImage', window, exmplImgs_wNoise);

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
            if string(KbName(firstPress)) == "Return"
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
    speech1 = ['Once you respond (or after 2 seconds if you don''t respond), the experiment will go on to the next trial. \n' ...
        'Try to identify the target as quickly and accurately as possible. \n\n' ...
        'Press ENTER/RETURN to Continue...'];
    DrawFormattedText(window, speech1,'center', (scrn_top+txtDistFrmTop4), black);
    
    % B) Display page until button is pressed
    keyIsDown = 0;
    enter_pressed = 0;
    swtch = 0;
    respns = "";

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
            if string(KbName(firstPress)) == "Return"
                enter_pressed = 1;
            end
        end
    end
    enter_pressed = 0;
    keyIsDown = 0;
    swtch = 0;
    % ---------------------------------------------------------------------

    Screen('FillRect', window, gray); % gray out window

end