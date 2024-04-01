function GenRevCorrBaseImsRF_v2
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
thisFile = "GenRevCorrBaseImsRF_v2.m"; % string should correspond to the name of this file's name..
mlab_dir = fileparts(which(thisFile));

% go to the matlab dir..
cd(mlab_dir)
% -------------------------------------
%% Parameters
    %% Output File Parameters
    %basename = "revcorrBI";
    basename = "BaseIm";
    out_fmt = "png";
    crop = 1;     % crop into square? 1=yes,0=no
    cropsz = 450; %output edge dimension after cropping into square
    method = "nearest"; % Interpolation Method
    
    %% System Parameters
    %--------------------------------------------------------------
    outdir = "TexDiffLim2Inv_v1"; % Specify output directory (which will be created)
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
    %% Common Region Parameters
    %--------------------------------------------------------------
    CR = 0; %Include Common Region Box?
    CR_penWidth = 2; % If so, specify pen width for common region box..
    CR_penLum = 0; % specify the luminance value for the common region pen strokes too..
    %CR_penLum=linspace(63.75,191.25,56);
    
    %CR_penLum=linspace(127.5,191.25,28);
    %CR_penLum=linspace(63.75,127.5,28);

    CR_al = 1; % Specify whether the common region should align with the objects (1) or be perpendicular (0)
    
    CR_obj = 2; % Specify which "object" the common region should surround (2 for upper object or 1 for lower object)
    %CR_obj = [1,2];
    
    cr_long = 32; % Long dimension of common region box expressed in % of screen height
    cr_short = 14; % Short dimension of common region box expressed in % of screen height
    
    if CR==0
        % If common region is not being used.. set defaults (even thogh they won't be used..) such that there are
        % not issues with linspace vectors breaking code..
        clear CR_penLum
        CR_penLum=0;
    end
    %% Texture Parameters
    %--------------------------------------------------------------
    T = 1; % Run texture version? (1=yes, 0=no)
    T_al = 1; % Specify whether the texture should align with the objects (1) or be perpendicular (0)
    
    image_in = 'spheres_pic_rs2_sqr.jpg';

    imdirBase = strcat(mlab_dir,"/","texImgs");

    image_in = strcat(imdirBase,"/",image_in); % add the image directory base path to the filename
    
    imName = "spheres";
    %rscale_f1 = 1.5; % texture scaling factor 1
    %rscale_f1 = linspace(0.4,2,8); % texture scaling factor 
    %rscale_f2 = 0.5; % texture scaling factor 2
    %rscale_f2 = linspace(0.4,2,8); % texture scaling factor 2

    % mean zoom constant at 1.. (use inverses..)
    rscale_f1 = linspace(1,2,201);  % texture scaling factor
    rscale_f2 = 1./rscale_f1; % texture scaling factor 2 (make all values the inverse of rscale_f1 counterparts..)

    if T==0
        % If texture is not being used.. set defaults (even thogh they won't be used..) such that there are
        % not issues with linspace vectors breaking code..
        clear rscale_f1
        clear rscale_f2
        rscale_f1=127;
        rscale_f2=127;
    end
    %% Luminance Parameters
    %--------------------------------------------------------------
    L = 0; % Run luminance version? (1=yes, 0=no)


    %lum1 = 90; % num between 0-255 90
    %lum1 = linspace(0,255,16);
    
    %LUM VARY
    %lum1 = linspace(0,255,8);

    %LUM VARY around background
    %lum1 = linspace(107.5,127.5,41); 
    lum1 = linspace(0,127.5,101); 
    
    %lum1 = linspace(127.5,127.5,1);
    
    %LUM CONST SLIGHT BRIGHT
    %lum1 = linspace(170,170,1);

    %lum2 = linspace(0,255,16); % num between 0-255 150
    
    %LUM VARY
    %lum2 = linspace(0,255,8); % num between 0-255 150
    
    %LUM VARY around background
    %lum2 = linspace(127.5,147.5,41); % num between 0-255 150
    lum2 = linspace(127.5,255,101); % num between 0-255 150


    %lum2 = linspace(127.5,127.5,1); % num between 0-255 150

    %LUM CONST SLIGHT BRIGHT
    %lum2 = linspace(170,170,1);

    %lumWt = 0.7; % luminance weight (for texture/luminance combo.. note: texture weight will be 1-lumWt..) (num between 0-1)
    lumWt = 1; % luminance weight (for texture/luminance combo.. note: texture weight will be 1-lumWt..) (num between 0-1)
    L_al = 1; % Specify whether the luminance should align with the objects (1) or be perpendicular (0)

    if L==0
        % If luminance is not being used.. set defaults (even thogh they won't be used..) such that there are
        % not issues with linspace vectors breaking code..
        clear lum1
        clear lum2
        lum1=127;
        lum2=127;
    end
    %% Parameters to Recursively Iterate Through Possible Combinations
    %Specify each of the parameter vectors you want to include
    
    % LUM ONLY
%     itrPars   = {occluders,oriCon,objcon,lum1,lum2};           % list names of parameter variables you want to iterate
%     itParNamz = {"occluders","oriCon","objcon","lum1","lum2"}; % list names of parameter variables you want to iterate (in itrParz above) as strings
%     avoidEqlz = {"lum1==lum2"}; % specify equality conditions you do now want
%     itrParsCom = combvec(itrPars{1},itrPars{2},itrPars{3},itrPars{4},itrPars{5}); % This builds a matrix of all the combos of pars in vectors listed in itrPars.
    
    % LUM ONLY NARROW RANGE, MEAN OBJECT LUM == BACKGROUND
    %itrPars   = {occluders,oriCon,objcon,lum1,lum2};           % list names of parameter variables you want to iterate
    %itParNamz = {"occluders","oriCon","objcon","lum1","lum2"}; % list names of parameter variables you want to iterate (in itrParz above) as strings
    %avoidEqlz = {"lum1==lum2","((lum1+lum2)/2)~=127.5"}; % specify equality conditions you do now want
    %itrParsCom = combvec(itrPars{1},itrPars{2},itrPars{3},itrPars{4},itrPars{5}); % This builds a matrix of all the combos of pars in vectors listed in itrPars.

    % TEXTURE ONLY
%     itrPars   = {occluders,oriCon,objcon,rscale_f1,rscale_f2,};           % list names of parameter variables you want to iterate
%     itParNamz = {"occluders","oriCon","objcon","rscale_f1","rscale_f2"}; % list names of parameter variables you want to iterate (in itrParz above) as strings
%     avoidEqlz = {"rscale_f1==rscale_f2"}; % specify equality conditions you do now want
%     itrParsCom = combvec(itrPars{1},itrPars{2},itrPars{3},itrPars{4},itrPars{5}); % This builds a matrix of all the combos of pars in vectors listed in itrPars.
    
% TEXTURE ONLY LIMIT MEAN TEXTURE DIFF..
    itrPars   = {occluders,oriCon,objcon,rscale_f1,rscale_f2};           % list names of parameter variables you want to iterate
    itParNamz = {"occluders","oriCon","objcon","rscale_f1","rscale_f2"}; % list names of parameter variables you want to iterate (in itrParz above) as strings
    %avoidEqlz = {"rscale_f1==rscale_f2","(rscale_f1*rscale_f2)~=1"}; % specify equality conditions you do now want
    Eqlz = {"rscale_f1~=rscale_f2","(rscale_f1*rscale_f2)==1"}; % specify equality conditions you DO want
    
    itrParsCom = combvec(itrPars{1},itrPars{2},itrPars{3},itrPars{4},itrPars{5}); % This builds a matrix of all the combos of pars in vectors listed in itrPars.

   
    [itrParsCom_matches] = combvecLogicalIdx(itrPars,itParNamz,Eqlz,itrParsCom);

    % Compute a logical index for the columns that meet the condition
    idx = (itrParsCom(4,:) .* itrParsCom(5,:)) == 1;
    
    % Extract the columns that meet the condition
    result = itrParsCom(:, idx);



% CR Lum Const
%     itrPars   = {occluders,oriCon,objcon,lum1,lum2,CR_penLum,CR_obj};           % list names of parameter variables you want to iterate
%     itParNamz = {"occluders","oriCon","objcon","lum1","lum2","CR_penLum","CR_obj"}; % list names of parameter variables you want to iterate (in itrParz above) as strings
%     avoidEqlz = {}; % specify equality conditions you do now want
%     itrParsCom = combvec(itrPars{1},itrPars{2},itrPars{3},itrPars{4},itrPars{5},itrPars{6},itrPars{7}); % This builds a matrix of all the combos of pars in vectors listed in itrPars.    


    %More complex version of above..
    %itrPars   = {occluders,oriCon,objcon,lum1,lum2,rscale_f1,rscale_f2,lumWt};           % list names of parameter variables you want to iterate
    %itParNamz = {"occluders","oriCon","objcon","lum1","lum2","rscale_f1","rscale_f2","lumWt"}; % list names of parameter variables you want to iterate (in itrParz above) as strings
    %avoidEqlz = {"rscale_f1==rscale_f2"}; % specify equality conditions you do now want
    %itrParsCom = combvec(itrPars{1},itrPars{2},itrPars{3},itrPars{4},itrPars{5},itrPars{6},itrPars{7},itrPars{8}); % This builds a matrix of all the combos of pars in vectors listed in itrPars.

    %Build a matrix of all the combinations 
    % (NOTE: If you add a vector to itrPars, you must add a corresponding
    % itrPars{n} to itrParsCom).. 
    %itrParsCom = combvec(itrPars{1},itrPars{2},itrPars{3},itrPars{4},itrPars{5}); % This builds a matrix of all the combos of pars in vectors listed in itrPars.
    %itrParsCom = combvec(itrPars{1},itrPars{2},itrPars{3},itrPars{4},itrPars{5},itrPars{6},itrPars{7},itrPars{8}); % This builds a matrix of all the combos of pars in vectors listed in itrPars.

                                                                                 % Each column is a condition/base image. Each row corresponds to the parameter value of the parameter contained in the corresponding itrPars index.                                                                            % Each column is a condition/base image. Each row corresponds to the parameter value of the parameter contained in the corresponding itrPars index.
    % build a corresponding cell array to store
    % image matrix outputs for each condition
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
for vv = 1:length(avoidEqlz)
    tstEqz = eval(avoidEqlz{vv});
    if tstEqz==1
        abortPass = abortPass+1;
    end
end
clear vv


%% Set Up Windows
if firstPass == 1
screenrect = Screen('Rect',0);	% get the size of the display screen

if change_screensize == 1
    screenrect = screenrect * scale_f;
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
imCons = {orient,objcon,occluders,change_screensize,scale_f,CR,CR_penWidth,CR_penLum,CR_al,CR_obj,cr_long,cr_short,T,T_al,image_in,imName,rscale_f1,rscale_f2,L,lum1,lum2,lumWt,L_al,window,outdir,outdir_base,basename,out_fmt,cmdLineOut,crop,cropsz};
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
imCons = {orient,objcon,occluders,change_screensize,scale_f,CR,CR_penWidth,CR_penLum,CR_al,CR_obj,cr_long,cr_short,T,T_al,image_in,imName,rscale_f1,rscale_f2,L,lum1,lum2,lumWt,L_al,window,outdir,outdir_base,basename,out_fmt,cmdLineOut,crop,cropsz};
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

%% Run the RevCorrBaseImGen_v2 recursively with the selected parameter in parsIn cell array
imOutMat = RevCorrBaseImGen_v2(parsIn,window,change_screensize,scale_f,imOutMat,method,baseOutDir,norm);

%% Bundle, Save Output .MAT, and Close Up
sca
cd(baseOutDir)
outmatname = strcat("RevCorBIs-",timestamp,".mat"); %make .mat filename
%bundle desired outputs into output structure

% CLEAN UP imOutMat such that it only includes images which are included in
% the output dataset.
% ========================================================================
%separate the 2 rows..
imOutMatT1=imOutMat(1,:);
imOutMatT2=imOutMat(2,:);

% create a logical arrays indicating which cells are not string/"NA"
logicalArray1 = cellfun(@(x) ~isstring(x), imOutMatT1);
logicalArray2 = cellfun(@(x) ~isstring(x), imOutMatT2);

notNA1=imOutMatT1(logicalArray1);
notNA2=imOutMatT2(logicalArray2);

imOutMatIncld=vertcat(notNA1,notNA2);
% ========================================================================
structOut = struct;
structOut.itParNamz = itParNamz;
structOut.itrPars = itrPars;
structOut.itrParsCom = itrParsCom;
structOut.imOutMat = imOutMat;
structOut.imOutMatIncld = imOutMatIncld;
structOut.imConParNames = "{orient,objcon,occluders,change_screensize,scale_f,CR,CR_penWidth,CR_penLum,CR_al,CR_obj,cr_long,cr_short,T,T_al,image_in,imName,rscale_f1,rscale_f2,L,lum1,lum2,lumWt,L_al,window,outdir,outdir_base,basename,out_fmt,cmdLineOut,crop,cropsz}";
structOut.avoidEqlz = avoidEqlz;

%save it
save(outmatname,"structOut",'-mat')
cd(strtDir); % make sure we end at the starting dir
toc
end