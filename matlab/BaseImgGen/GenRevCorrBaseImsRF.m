function GenRevCorrBaseImsRF
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
thisFile = "GenRevCorrBaseImsRF.m"; % string should correspond to the name of this file's name..
mlab_dir = fileparts(which(thisFile));

% go to the matlab dir..
cd(mlab_dir)
% -------------------------------------
%% Parameters
    %% Output File Parameters
    basename = "revcorrBI";
    out_fmt = "png";
    crop = 1;     % crop into square? 1=yes,0=no
    cropsz = 512; %output edge dimension after cropping into square

    
    %% System Parameters
    %--------------------------------------------------------------
    outdir = "LumTest"; % Specify output directory (which will be created)
    outdir_base = mlab_dir; % specify base directory in which outdir will be created/saved.. (currently set to the base directory auto-detected above under "Set up path" to avoid path annoyances across machines accessing shared folder..)
    
    %% Condition-General Parameters
    %--------------------------------------------------------------
    orientStrz = ['m','z']; % possible string values: 'm' or 'z'
    oriCon = [1,2]; % corresponding numeric values (for coding simplicity)
    orient = 'm'; % orientation condition
    objcon = [1,2]; % 1 or 2
    occluders = [0,1]; % include occluders? 0 = yes, 1 = no
    
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
    %lum1 = 90; % num between 0-255 90
    
    %lum1 = linspace(0,255,16);
    lum1 = linspace(0,255,11);

    %lum2 = linspace(0,255,16); % num between 0-255 150
    lum2 = linspace(0,255,11); % num between 0-255 150

    %lumWt = 0.7; % luminance weight (for texture/luminance combo.. note: texture weight will be 1-lumWt..) (num between 0-1)
    lumWt = 1; % luminance weight (for texture/luminance combo.. note: texture weight will be 1-lumWt..) (num between 0-1)
    L_al = 1; % Specify whether the luminance should align with the objects (1) or be perpendicular (0)
    
    %% Parameters to Recursively Iterate Through Possible Combinations
    %Specify each of the parameter vectors you want to include
    itrPars   = {occluders,oriCon,objcon,lum1,lum2};           % list names of parameter variables you want to iterate
    itParNamz = {"occluders","oriCon","objcon","lum1","lum2"}; % list names of parameter variables you want to iterate (in itrParz above) as strings

    avoidEqlz = {"lum1==lum2","(lum1+lum2)/2~=127.5"}; % specify equality conditions you do now want

    %Build a matrix of all the combinations 
    % (NOTE: If you add a vector to itrPars, you must add a corresponding
    % itrPars{n} to itrParsCom).. 
    itrParsCom = combvec(itrPars{1},itrPars{2},itrPars{3},itrPars{4},itrPars{5}); % This builds a matrix of all the combos of pars in vectors listed in itrPars.
                                                                                 % Each column is a condition/base image. Each row corresponds to the parameter value of the parameter contained in the corresponding itrPars index.
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

if abortPass == 0
%% Set Up Windows
screenrect = Screen('Rect',0);	% get the size of the display screen

if change_screensize == 1
    screenrect = screenrect * scale_f;
end
scrn_wid = screenrect(3);
scrn_hei = screenrect(4);
scrnWidHei = [scrn_wid,scrn_hei];

Screen('Preference', 'SkipSyncTests', 1); %EJD: added to skip sync tests
window = PsychImaging('OpenWindow', 0, GrayIndex(0,0.5),[]);

% PsychImaging('PrepareConfiguration');
% PsychImaging('AddTask', 'General', 'UseFastOffscreenWindows');
Screen(window,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);    % allow for transparency!

%% Run the RevCorrBaseImGen_v1 with the selected parameters..
imCons = {orient,objcon,occluders,change_screensize,scale_f,CR,CR_penWidth,CR_penLum,CR_al,CR_obj,cr_long,cr_short,T,T_al,image_in,imName,rscale_f1,rscale_f2,L,lum1,lum2,lumWt,L_al,window,outdir,outdir_base,basename,out_fmt,cmdLineOut,crop,cropsz};
imgSnap = RevCorrBaseImGen_v1(orient,objcon,occluders,change_screensize,scale_f,CR,CR_penWidth,CR_penLum,CR_al,CR_obj,cr_long,cr_short,T,T_al,image_in,imName,rscale_f1,rscale_f2,L,lum1,lum2,lumWt,L_al,window,outdir,outdir_base,basename,out_fmt,cmdLineOut,crop,cropsz);
imOutMat{1,ii} = imgSnap;
imOutMat{2,ii} = imCons;
else
disp("Unwanted Parameter Combo Condition Met. Skipping this Image...")
imOutMat{1,ii} = "NA";
imOutMat{2,ii} = {};
end
cd(strtDir); % make sure we end at the starting dir
end
%% Bundle, Save Output .MAT, and Close Up
sca
cd(baseOutDir)
outmatname = strcat("RevCorBIs-",timestamp,".mat"); %make .mat filename
%bundle desired outputs into output structure
structOut = struct;
structOut.itParNamz = itParNamz;
structOut.itrPars = itrPars;
structOut.itrParsCom = itrParsCom;
structOut.imOutMat = imOutMat;
structOut.imConParNames = "{orient,objcon,occluders,change_screensize,scale_f,CR,CR_penWidth,CR_penLum,CR_al,CR_obj,cr_long,cr_short,T,T_al,image_in,imName,rscale_f1,rscale_f2,L,lum1,lum2,lumWt,L_al,window,outdir,outdir_base,basename,out_fmt,cmdLineOut,crop,cropsz}";
structOut.avoidEqlz = avoidEqlz;

%save it
save(outmatname,"structOut",'-mat')
cd(strtDir); % make sure we end at the starting dir
toc
end