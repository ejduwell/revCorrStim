function GenRevCorrBaseImsRF_v4
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
% should always be the case to avoid other conflicts..)
thisFile = "GenRevCorrBaseImsRF_v4.m"; % string should correspond to the name of this file's name..
mlab_dir = fileparts(which(thisFile));

% go to the matlab dir..
cd(mlab_dir);
% -------------------------------------
%% Parameters
    %% Output File Parameters
    %basename = "revcorrBI";
    basename = "BaseIm";
    out_fmt = "png";
    crop = 1;     % crop into square? 1=yes,0=no
    cropsz = 512; %output edge dimension after cropping into square
    method = "nearest"; % Interpolation Method
    
    %% System Parameters
    %--------------------------------------------------------------
    %outdir = "test3_CROnlyWlum51"; % Specify output directory (which will be created)
    outdir = "test5_TexOnly"; % Specify output directory (which will be created)

    %outdir = "LumTest5WhiteNoiseRan"; % Specify output directory (which will be created)
    cd("..");
    cd("..");
    cd('images');
    outdir_base = pwd; % specify base directory in which outdir will be created/saved.. 
    %outdir_base = mlab_dir; % specify base directory in which outdir will be created/saved.. (currently set to the base directory auto-detected above under "Set up path" to avoid path annoyances across machines accessing shared folder..)
    % go back to the matlab dir..
    cd(mlab_dir)
    %% Condition-General Parameters
    %--------------------------------------------------------------
    orientStrz = ['m','z']; % possible string values: 'm' or 'z'
    oriCon = [1,2]; % corresponding numeric values (for coding simplicity)
    orient = 'm'; % orientation condition
    objcon = [1,2]; % 1 or 2
    occluders = [0,1]; % include occluders? 0 = yes, 1 = no
    
    change_screensize = 1; % if 1, will resize the screen dimensions by scale factor "scale_f" defined below.., if 0, does nothing
    scale_f = 0.5; % Screen rescale factor (number between 0 and 1) 
    cmdLineOut =0; % if 1, will print the parameters used for creating the base image out on the command line
    
    norm=0; % normalization in PORClinkLTCR_v1 .. 1=normalize before output, 0=don't normalize.. #Probably only want to normalize on conditions where there are combos of lum/texture/etc that could result in out-of range values...
    
    % Fixation/Inducer Parameters
    indcrPenW=4; % sets pen width for the "inducer" lines
    fxnSqrSclr=50; % this par controls the scaling of the central fixation square.. Note: its proportion of the "main dimension" of the stimulus.. ie mainDim/fxnSqrSclr.. Bigger fxnSqrSclr-->smaller square
    fxnParz=[indcrPenW,fxnSqrSclr];
    
    %% Common Region Parameters
    %--------------------------------------------------------------
    CR = 0; %Include Common Region Box?
    CR_penWidth = 4; % If so, specify pen width for common region box..
    cr_CornerRad=60;% controls the radius of the rounded corners of the common region box
    %CR_penLum = 0; % specify the luminance value for the common region pen strokes too..
    %CR_penLum=linspace(63.75,191.25,56);
    
    %CR_penLum=linspace(127.5,255,101);
    CR_penLum=linspace(0,255,101);

    CR_al = 1; % Specify whether the common region should align with the objects (1) or be perpendicular (0)
    
    %CR_obj = 2; % Specify which "object" the common region should surround (2 for upper object or 1 for lower object)
    CR_obj = [1,2];
    
    cr_long = 33; % Long dimension of common region box expressed in % of screen height
    cr_short = 15; % Short dimension of common region box expressed in % of screen height
    
    if CR==0
        % If common region is not being used.. set defaults (even thogh they won't be used..) such that there are
        % not issues with linspace vectors breaking code..
        clear CR_penLum
        CR_penLum=0;
        cr_CornerRad=60;
        CR_obj = "NA";
    end
    %% Texture Parameters
    %--------------------------------------------------------------
    T = 1; % Run texture version? (1=yes, 0=no)
    T_al = 1; % Specify whether the texture should align with the objects (1) or be perpendicular (0)
    
    %image_in = 'checkerboard_rs10_nearest.jpg';
    image_in = 'checkerBrd_5000by5000Array_50by50_Checks.png';

    imdirBase = strcat(mlab_dir,"/","texImgs");

    image_in = strcat(imdirBase,"/",image_in); % add the image directory base path to the filename
    
    imName = "checkrz";
    %rscale_f1 = 1.5; % texture scaling factor 1
    %rscale_f1 = linspace(0.4,2,8); % texture scaling factor 
    %rscale_f2 = 0.5; % texture scaling factor 2
    %rscale_f2 = linspace(0.4,2,8); % texture scaling factor 2

    % mean zoom constant at 1.. (use inverses..)
%     rscale_f1 = linspace(1,2,101);  % texture scaling factor
%     rscale_f2 = 1./rscale_f1; % texture scaling factor 2 (make all values the inverse of rscale_f1 counterparts..)

    % HEADS UP: THE MAX POSSIBLE RESCALE VALUE FOR F1/MIN FOR F2 WILL
    % DEPEND ON THE SIZE OF THE "TEXTURE IMAGE" AND SIZE OF YOUR SCREEN
    % WINDOW YOU'RE DRAWING ON.. IF YOU GET AN ERROR SAYING:
    % Error using centerCropWindow2d
    % Target size may not be larger than input image size in any dimension.
    % ...
    % Then Thats whats going on..
    rscale_f1 = linspace(1,2.5,101);  % texture scaling factor
    rscale_f2 = 1./rscale_f1; % texture scaling factor 2 (make all values the inverse of rscale_f1 counterparts..)

    if T==0
        % If texture is not being used.. set defaults (even thogh they won't be used..) such that there are
        % not issues with linspace vectors breaking code..
        clear rscale_f1
        clear rscale_f2
        rscale_f1=1;
        rscale_f2=1;
    end
    %% Luminance Parameters
    %--------------------------------------------------------------
    L = 0; % Run luminance version? (1=yes, 0=no)
    
    %LUM VARY around background (Larger # of steps) (FOR LUM ONLY)
    %lum1 = linspace(0,127.5,127);  
    %lum2 = linspace(127.5,255,127); % num between 0-255 150
    
    %LUM1&2 BOTH VARY SAME RANGE 0-255 (FOR CR w/LUM)  
    lum1 = linspace(0,255,11);  
    lum2 = linspace(0,255,11); % num between 0-255 150

    lumWt = 1; % luminance weight (for texture/luminance combo.. note: texture weight will be 1-lumWt..) (num between 0-1)
    L_al = 1; % Specify whether the luminance should align with the objects (1) or be perpendicular (0)

    if L==0
        % If luminance is not being used.. set defaults (even though they won't be used..) such that there are
        % not issues with linspace vectors breaking code..
        clear lum1
        clear lum2
        lum1=127;
        lum2=127;
    end
    %% Parameters to Recursively Iterate Through Possible Combinations
    %Specify each of the parameter vectors you want to include
    
    % LUM ONLY NARROW RANGE, MEAN OBJECT LUM == BACKGROUND
    if L==1 && T==0 && CR==0
    itrPars   = {occluders,oriCon,objcon,lum1,lum2};           % list names of parameter variables you want to iterate
    itParNamz = {"occluders","oriCon","objcon","lum1","lum2"}; % list names of parameter variables you want to iterate (in itrParz above) as strings
    Eqlz = {"lum1~=lum2","((lum1+lum2)/2)==127.5"}; % specify equality conditions you want
    itrParsCom = combvec(itrPars{1},itrPars{2},itrPars{3},itrPars{4},itrPars{5}); % This builds a matrix of all the combos of pars in vectors listed in itrPars.
    end

    % TEXTURE ONLY
    if T==1 && L==0 && CR==0
    itrPars   = {occluders,oriCon,objcon,rscale_f1,rscale_f2,};           % list names of parameter variables you want to iterate
    itParNamz = {"occluders","oriCon","objcon","rscale_f1","rscale_f2"}; % list names of parameter variables you want to iterate (in itrParz above) as strings
    Eqlz = {"rscale_f1~=rscale_f2","(rscale_f1*rscale_f2)==1"}; % specify equality conditions you do now want
    itrParsCom = combvec(itrPars{1},itrPars{2},itrPars{3},itrPars{4},itrPars{5}); % This builds a matrix of all the combos of pars in vectors listed in itrPars.
    end

    % COMMON REGION WITH LUM ONLY W/IN OBJECTS BUT NO LUM DIFF BTW
    % OBJECTS.. (IE ONLY CR CUE CAN DISCRIMINATE..)
    if T==0 && L==1 && CR==1
    itrPars   = {occluders,oriCon,objcon,lum1,lum2,CR_penLum,CR_obj};           % list names of parameter variables you want to iterate
    itParNamz = {"occluders","oriCon","objcon","lum1","lum2","CR_penLum","CR_obj"}; % list names of parameter variables you want to iterate (in itrParz above) as strings
    Eqlz = {"abs(CR_penLum-127.5)>0","lum1==lum2","lum1==51.0"}; % specify equality conditions you do now want
    itrParsCom = combvec(itrPars{1},itrPars{2},itrPars{3},itrPars{4},itrPars{5},itrPars{6},itrPars{7}); % This builds a matrix of all the combos of pars in vectors listed in itrPars.    
    end
   
    %More complex version of above..
    %itrPars   = {occluders,oriCon,objcon,lum1,lum2,rscale_f1,rscale_f2,lumWt};           % list names of parameter variables you want to iterate
    %itParNamz = {"occluders","oriCon","objcon","lum1","lum2","rscale_f1","rscale_f2","lumWt"}; % list names of parameter variables you want to iterate (in itrParz above) as strings
    %avoidEqlz = {"rscale_f1==rscale_f2"}; % specify equality conditions you do now want
    %itrParsCom = combvec(itrPars{1},itrPars{2},itrPars{3},itrPars{4},itrPars{5},itrPars{6},itrPars{7},itrPars{8}); % This builds a matrix of all the combos of pars in vectors listed in itrPars.

    % Generate set of combos from the specifications above and
    % reduce to set matching the equality statements..
    itrParsCom_all=itrParsCom; % save copy of all combinations
    [itrParsCom_matches] = combvecLogicalIdx(itrPars,itParNamz,Eqlz,itrParsCom); % reduce to just the ones which match all the logical conditions in Eqlz
    itrParsCom = itrParsCom_matches; %reassign itrParsCom to itrParsCom_matches (was easier to update code this way.. lazy Ethan is lazy..)
    imOutMat = num2cell(zeros(2,size(itrParsCom,2)));

    SubDirStruc=struct;
    SubDirStruc.levels = {"occluders","oriCon"};
    SubDirStruc.levelVals = {{0,1},{1,2}};
    SubDirStruc.levelDirs = {{"occ","nocc"},{"R","L"}};
    
    % make unique timestamp
    timestamp=datestr((datetime('now')));
    timestamp= strrep(timestamp," ","-");
    timestamp= strrep(timestamp,":","-");
    baseOutDir = strcat(outdir_base,'/',outdir,"-",timestamp); % save the outdir base level ..

%% Set up output subdirectory structure
    mkdir(baseOutDir); % make base directory
    cd(baseOutDir); % enter it
    for nsubdir=1:length(SubDirStruc.levelDirs{1})
        mkdir(SubDirStruc.levelDirs{1}{(nsubdir)}); % make subdirectory
        cd(SubDirStruc.levelDirs{1}{(nsubdir)}); %enter it
        for nsubsubdir=1:length(SubDirStruc.levelDirs{2})
            mkdir(SubDirStruc.levelDirs{2}{(nsubsubdir)}); % make subdirectory
        end
        cd(baseOutDir); % go back to base output level
    end
    clear nsubdir
    clear nsubsubdir
    % go back to the matlab dir..
    cd(mlab_dir)

imNum = 0; %initialize imNum counter
disp("#####################################")
disp("Beginning Main Image Generation Loop:")
disp(strcat("Will generate a total of"," ",num2str(size(itrParsCom,2))," ","images."))
disp("#####################################")

firstPass = 1; % initialize first pass detector.
tic
%% Begin Main Image Generation Loop    
for ii = 1:size(itrParsCom,2) % Main for loop: loops through par combos (columns in itrParsCom)                                                            
imNum = imNum+1; % set imNum for this pass
abortPass=0; %initialize abortPass for this pass
disp(strcat("Image"," ",num2str(imNum)," ","of"," ",num2str(size(itrParsCom,2)),"..."))
%% Update Parameters for This Pass
for prz = 1:length(itParNamz)
    parStr    = itParNamz{prz};             % store par name to update as string
    parValStr = num2str(itrParsCom(prz,ii));         % get corresponding value from itrParsCom as string
    parCmd = strcat(parStr,"=",parValStr,";");  % build command string useing strcat to update parameter
    eval(parCmd);                           % evaluate string command with eval
end

% If "oriCon" is listed in itParNamz:
% Update "orient" string value based on current "oriCon" (set above as numeric proxy for orient)
if ~any(strcmp(itParNamz,"oriCon"))
orient = orientStrz(oriCon); 
end

% Update output directory to appropriate level/sublevel..
%----------------------------------------------------------
subDirTmp = {}; %initialize temporary cell
for nsubdir=1:length(SubDirStruc.levelVals)
    for zz = 1: length(SubDirStruc.levelVals{nsubdir})
    conChkCmd=strcat(SubDirStruc.levels{nsubdir},"==",num2str(SubDirStruc.levelVals{nsubdir}{(zz)}),";");
    conChkCmdOut=eval(conChkCmd);
    if conChkCmdOut == 1
        subDirTmp{end+1} = SubDirStruc.levelDirs{nsubdir}{(zz)};
    end
    end
    clear zz
end

strTmp="/"; % initialize
for xx = 1:length(subDirTmp)
    strTmp = strcat(strTmp,subDirTmp{xx},"/");
end
clear xx

outdir = strcat(baseOutDir,strTmp); % update output dir variable..
% Clean up
clear subDirTmp;
clear strTmp;
%----------------------------------------------------------

% Check whether this is an unwanted equality condition
%for vv = 1:length(avoidEqlz)
%    tstEqz = eval(avoidEqlz{vv});
%    if tstEqz==1
%        abortPass = abortPass+1;
%    end
%end
%clear vv


%% Set Up Windows
if firstPass == 1
screenrect = Screen('Rect',0);	% get the size of the display screen

if change_screensize == 1
    comp = Screen('Computer');
    if strcmp(comp.machineName,'tron')
        scale_f = 0.5;
        % ejd attempt to rescale whole double screen into a single screen..
        %screenrect = screenrect * scale_f; %ORIG
        screenrect(3) = screenrect(3) * scale_f;

        scl2crpFctr=cropsz/screenrect(4);
        screenrect(3)=(cropsz*2);
        screenrect(4)=(cropsz*2);
        %screenrect(3)=round(screenrect(3)*scl2crpFctr);
        %screenrect(4)=round(screenrect(4)*scl2crpFctr);
        pause = "";
    else
        scale_f = 0.5;
        screenrect = screenrect * scale_f;
    end
end
scrn_wid = screenrect(3);
scrn_hei = screenrect(4);
scrnWidHei = [scrn_wid,scrn_hei];

Screen('Preference', 'SkipSyncTests', 1); %EJD: added to skip sync tests
window = PsychImaging('OpenWindow', 0, GrayIndex(0,0.5),screenrect);

% PsychImaging('PrepareConfiguration');
% PsychImaging('AddTask', 'General', 'UseFastOffscreenWindows');
Screen(window,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);    % allow for transparency!

%set up par cell structure
imCons = {orient,objcon,occluders,change_screensize,scale_f,CR,CR_penWidth,CR_penLum,CR_al,CR_obj,cr_long,cr_short,T,T_al,image_in,imName,rscale_f1,rscale_f2,L,lum1,lum2,lumWt,L_al,window,outdir,outdir_base,basename,out_fmt,cmdLineOut,crop,cropsz,cr_CornerRad};
parsIn =cell(size(itrParsCom,2),length(imCons));

%SET UP NA ROW
naRow = cell(1,length(imCons));
for nn=1:length(naRow)
naRow{1,nn}="NA";
end

firstPass = 0; % flip off first pass switch..
end
%% Add rows to parsIn cell array
if abortPass == 0
imCons = {orient,objcon,occluders,change_screensize,scale_f,CR,CR_penWidth,CR_penLum,CR_al,CR_obj,cr_long,cr_short,T,T_al,image_in,imName,rscale_f1,rscale_f2,L,lum1,lum2,lumWt,L_al,window,outdir,outdir_base,basename,out_fmt,cmdLineOut,crop,cropsz,cr_CornerRad};
parsIn(ii,:)=imCons;

%EJD TEMP COMMENT..
%imgSnap = RevCorrBaseImGen_v1(orient,objcon,occluders,change_screensize,scale_f,CR,CR_penWidth,CR_penLum,CR_al,CR_obj,cr_long,cr_short,T,T_al,image_in,imName,rscale_f1,rscale_f2,L,lum1,lum2,lumWt,L_al,window,outdir,outdir_base,basename,out_fmt,cmdLineOut,crop,cropsz);
%imOutMat{1,ii} = imgSnap;
%imOutMat{2,ii} = imCons;

else
disp("Unwanted Parameter Combo Condition Met. Skipping this Image...")
parsIn(ii,1:end)=naRow;
% imOutMat{1,ii} = "NA";
% imOutMat{2,ii} = {};
end
% cd(strtDir); % make sure we end at the starting dir
end

%% Run the RevCorrBaseImGen_v3 recursively with the selected parameter in parsIn cell array
imOutMat = RevCorrBaseImGen_v3(parsIn,window,change_screensize,scale_f,imOutMat,method,baseOutDir,norm,fxnParz);

%% Bundle, Save Output .MAT, and Close Up
%sca
Screen('CloseAll');
Screen('Close');

cd(baseOutDir)
outmatname = strcat("RevCorBIs-",timestamp,".mat"); %make .mat filename
%bundle desired outputs into output structure

% CLEAN UP imOutMat such that it only includes images which are included in
% the output dataset.
% ========================================================================
% %separate the 2 rows..
% imOutMatT1=imOutMat(1,:);
% imOutMatT2=imOutMat(2,:);
% 
% % create a logical arrays indicating which cells are not string/"NA"
% logicalArray1 = cellfun(@(x) ~isstring(x), imOutMatT1);
% logicalArray2 = cellfun(@(x) ~isstring(x), imOutMatT2);
% 
% notNA1=imOutMatT1(logicalArray1);
% notNA2=imOutMatT2(logicalArray2);
% 
% imOutMatIncld=vertcat(notNA1,notNA2);
% ========================================================================
structOut = struct;
structOut.itParNamz = itParNamz;
structOut.itrPars = itrPars;
structOut.itrParsCom = itrParsCom;
structOut.imOutMat = imOutMat;
structOut.itrParsCom_all = itrParsCom_all;
structOut.imConParNames = "{orient,objcon,occluders,change_screensize,scale_f,CR,CR_penWidth,CR_penLum,CR_al,CR_obj,cr_long,cr_short,T,T_al,image_in,imName,rscale_f1,rscale_f2,L,lum1,lum2,lumWt,L_al,window,outdir,outdir_base,basename,out_fmt,cmdLineOut,crop,cropsz}";
structOut.Eqlz = Eqlz;

%save it
save(outmatname,"structOut",'-mat')
cd(strtDir); % make sure we end at the starting dir
toc
end