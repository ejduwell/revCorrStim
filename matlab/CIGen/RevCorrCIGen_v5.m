%function RevCorrCIGen_v4()

% EJD Tried moving the bp filter operation to only on the mean images.. 
% Thinks/testing whether this is equivilant to applying to each individual
% noise image then taking mean (as originally done in v3)
%
% Also Added the ability to "flip" the noise frames from one orientation
% condition prior to the construction of CIs
%% Auto-detect and set package-level directory/path params:
% -------------------------------------
% get matlab directory path by "which-ing" for this file..
% NOTE: key assumption: there is only one copy of this on the path.. (that
% should always be the case to avoid other conflicts..)
mlab_dir = fileparts(which('RevCorrCIGen_v5.m'));
startDir=pwd;
cd(mlab_dir);
% cd 2 directory to enter the main dir..
cd ..
cd ..
% save the path..
main_dir = pwd;

% go back to the matlab dir..
cd(mlab_dir)
% ------------------------------------

%% Set Parameters

% SUBJECT AND SESSION PARS
subjID="DVH"; % 3 letter subject ID
noiseType="krnlNz"; % which noise (white or krnlNz)
stimType="lum"; % which grouping parameter type (lum,cr, or tex)
stimVer="both"; % which version (occ, nocc, or both)

% RESPONSE BUTTONS FOR THE 2 ORI CONDITIONS
conRespz=['m','z']; % response options for the two conditions (right and left angle here..)

oconSplit=1; % if 1, will seperate off only one of the 2 object conditions..
oconVal=2; % object condition to use oconSplit is on if  1 or 2

% Specify CI overlay color parameters
clrMap="parula"; % colormap used on colorized/thresholded CI overlays for real and ideal CIs over the base image..
alphVal=0.5; % alpha/transparency value used on thresholded/colorized CI overlay..
% specify weights and thresholds
ciThr=75; % percent threshold used for thesholded/colorized CI and ideal CI overlays..
useSmthVer4clrd=0; % if 1, will colorize/threshold smoothed version instead.

% BASE IMAGES:
path2BIz=cell(1,2); % base image 1 and 2 path array..
path2BIz{1,1}="/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/images/test4_LumOnly-05-Mar-2024-10-12-48/occ/L/BaseIm_occ_0_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3214_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_0_lum2_255_LWT_1_Lal_1_Ocon_1.png"; %path to base image file
path2BIz{1,2}="/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/images/test4_LumOnly-05-Mar-2024-10-12-48/nocc/L/BaseIm_occ_1_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3214_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_0_lum2_255_LWT_1_Lal_1_Ocon_1.png"; %path to base image file
desiredDimz=[512,512]; % bi/ci dims
biIn1 = imgReformater(imread(path2BIz{1,1}),desiredDimz); % read in and format base image to desired dimensions
biIn2 = imgReformater(imread(path2BIz{1,2}),desiredDimz); % read in and format base image to desired dimensions

% Flip Pars
oriConFlp=1; % if 1, will flip all noise from one of the two conditions all on axis specified by before making CI.. (adam's request)
flpAxis=2; % axis along which flip will occur if oriConFlp is on/=1. (1=vertical, 2=horiz)
FlpdOriCon=1; % which orientation condition gets flipped? (1 or 2)

% SMOOTHING PARS
%--------------------------------------------------------------------------
% gaussian smoothing kernel size/sigma
smthKrnlSigma=5;
% filtering/smoothing domain : 'auto', 'spatial', or 'frequency'
fDomain='auto'; 
%--------------------------------------------------------------------------

% ------------                 BP Filter Pars:                 ------------
%==========================================================================
% Band Pass Filter On/Off
bpFilter=1; % 0=don't bandpass filter; 1=apply bandpass filter to noise images..

% ------------ LOWER & UPPER BOUNDS OF BAND PASS FILTER WINDOW ------------
lowerBnd=2; % lower frequency bound
upperBnd=8; % upper frequency bound
filterOrder=4; % filter order
pltFig=0;

% LOW PASS #1
%--------------------------------------------------------------------------
lowFreq=0;  % lower frequency bound of bandpass filter..
highFreq=5; % upper frequency bound of bandpass filter..
%--------------------------------------------------------------------------
% LOW PASS #2
%--------------------------------------------------------------------------
%lowFreq=0;  % lower frequency bound of bandpass filter..
%highFreq=32; % upper frequency bound of bandpass filter..
%--------------------------------------------------------------------------
% LOW PASS #3
%--------------------------------------------------------------------------
%lowFreq=0;  % lower frequency bound of bandpass filter..
%highFreq=16; % upper frequency bound of bandpass filter..
%--------------------------------------------------------------------------
% HIGH PASS
%--------------------------------------------------------------------------
% lowFreq=16;
% highFreq=256;
%--------------------------------------------------------------------------
%==========================================================================

dataDir=strcat(main_dir,"/data_master/",subjID,"/",noiseType,"/",stimType,"/",stimVer);

%% Read in all data from subject's tdfs
disp(" ");
disp("Reading in subject's TDFs...")
disp(" ");
tpoint11=datetime;
% Enter subject's data directory
% =========================================================================
try
    cd(dataDir);
catch
    disp(" ")
    disp("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    disp("!!!!!!!!!!!!! ERROR ENTERING/FINDING DIRECTORY !!!!!!!!!!!!!!!!!")
    disp("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    disp(" ")
    assert(exist(dataDir, 'dir'),strcat("Directory: '",dataDir,"' does not exist.."));
end
% =========================================================================

% Get list of subdirectories and make list of previously completed 
% quest output subdirectories
% =========================================================================
% 1. Get a list of the subdirectories..
folder = pwd;
subDs = GetSubDirsFirstLevelOnly(folder);

% 2. make list of previously completed quest output subdirectories
pattern = "_cmpltd";
qst_subDs = {};
numFolders = 0;	% keep track of how many folders we find.
for k = 1 : size(subDs,1)
    % Get this folder.% Example setup

    thisFolderName = subDs{k,1};
    % See if it contains our required pattern.
    if ~contains(thisFolderName, pattern, 'IgnoreCase', true)
        % Skip this file because the filename does not contain the required pattern.
        continue;
    end
    % The pattern is in the filename if you get here, so do something with it.
    %fprintf('Now processing %s\n', thisFileName);
    qst_subDs = vertcat(qst_subDs,thisFolderName);
    numFolders = numFolders + 1;	% For fun, let's keep track of how many files we processed.
end
disp(" ");
fprintf('We found %d folders with %s in the name.\n',numFolders, pattern);
disp("Entering each completed session folder and extracting the TDFs...");
disp(" ");

% =========================================================================

% sort in ascending order from oldest to latest
% (folders contain timestamp in name)
qst_subDs = sortrows(qst_subDs,'ascend'); 

% Enter completed session directory and load the tdfs into a struct
% =========================================================================
allTDFs=struct; % pre-allocate
strtDirTmp=pwd; % save starting location
for bb = 1: size(qst_subDs,1)
    dirTmp=qst_subDs{bb,1};
    cd(dirTmp);
    tmpStrct=load(strcat(subjID,".mat"));
    tdfs=tmpStrct.tdfs;
    strTmp=strcat("block",num2str(bb));
    allTDFs.(strTmp).tdfs=tdfs;
    clear tdfs
    clear tmpStrct;
    cd(strtDirTmp); % return to starting location..
end
clear tmpStrct
clear strTmp
% =========================================================================

% Combine the data from all TDFs into a single array..
% =========================================================================

% grab the top row of the first tdf in block1 and save as header
% NOTE: assuming that all tdfs have same headers
% This should be the case for all tdfs within a particular stimulus
% type/version/noise combo for each subject..
headerz=allTDFs.block1.tdfs{1,1}(1,:);

%preallocate/create fullTDF cell variable to contain data from all TDFs
% (Starting with just the headers.. will concatenate on from there..)
fullTDF=headerz;
nColz=size(fullTDF,2);
fullTDF{1,nColz+1}="blockNumber";
nColz=nColz+1;
blkCountr=1; % initialize countr for blocks..
for bb=1:numFolders
    strTmp=strcat("block",num2str(bb));
    tdfs=allTDFs.(strTmp).tdfs;
    for ii=1:size(tdfs,2)
        thisTDF=tdfs{1,ii};
        thisTDF(:,end+1)={blkCountr}; % add block number in final additional column
        fullTDF=vertcat(fullTDF,thisTDF(2:end,:)); % concatenate all but the top header row onto the fullTDF array
        blkCountr=blkCountr+1; % update block countr
    end
    allTDFs.(strTmp).tdfs={};% in the interest of memory clear these out as we go..
end
% clear un-needed variables..
clear tdfs
clear thisTDF
clear blkCountr
clear allTDFs

fullTDF(2:end,5)={"removed for RAM"}; % REMOVE BASE IMAGES FROM ARRAY FOR RAM PURPOSES
% =========================================================================
tpoint12=datetime;
time1=tpoint12-tpoint11;
disp(strcat("Reading in subject's TDFs took ",num2str(seconds(time1)), " seconds..."));
disp(" ");

%% Catagorize the Noise Images Based on Condition and Correctness
disp(" ");
disp("Catagorizing noise images based on condition and correctness...")
disp(" ");
tpoint21=datetime;


% Add object condition column (ocon) by using file names..
fullTDF{1,end+1}="ObjectCon";
for hh=2:size(fullTDF,1)
    parTagzTmp=["Ocon"];
    [parValzTmp, parTagzTmp] = xtractParsFrmFilename(fullTDF{hh,1},parTagzTmp,"num");
    fullTDF{hh,end}=parValzTmp;
end



% temporarily separate headers from data..
headerz=fullTDF(1,:);
fullTDF=fullTDF(2:end,:);


% seperate off just one object condition if specified..
if oconSplit==1
fullTDF_oconIdx=cell2mat(fullTDF(:,end))==oconVal;
fullTDF = fullTDF(fullTDF_oconIdx,:);
end

% Use logical indexing to separate data into the 4 constituent 
% condition/correctness combination conditions ..
% =========================================================================
BIConColNum = 7; % ground truth column
subjRespCol= 11; % subject response column
noiseImgCol=18; % column containing noise images..
occConCol=3; % Column with the occlusion condition info

% Create a logical indices for rows corresponding to the 4 possible
% condition/correctness combos..
%-------------------------------------------------------------------------

biCon1Idx = cell2mat(fullTDF(:, BIConColNum)) == conRespz(1); % Base image condition1
biCon2Idx = cell2mat(fullTDF(:, BIConColNum)) == conRespz(2); % Base image condition2
%rspCon1Idx = cell2mat(fullTDF(:, subjRespCol)) == conRespz(1);% response condition1

respCol=cellstr(fullTDF(:, subjRespCol));
respCol=strrep(respCol,'NA', '0'); % replace 'NA's with '0' to keep dims same

% replace any other responses that aren't an 'm' or 'z' with a '0' too
uniqResp=unique(respCol);
for zz=1:size(uniqResp,1)
    rspTmp=uniqResp{zz,1};
    if ((string(rspTmp) ~= "z") && (string(rspTmp) ~= "m")) % if its not a 'z' and its also not an 'm'
        respCol=strrep(respCol,rspTmp, '0'); % then replace all instances of it with a '0'
    end
end

respCol=cell2mat(respCol);
rspCon1Idx = respCol == cell2mat({conRespz(1)});% response condition1
rspCon2Idx = respCol == conRespz(2);% response condition2
%rspCon1Idx = cell2mat(cellstr(fullTDF(:, subjRespCol))) == cell2mat({conRespz(1)});% response condition1
%rspCon2Idx = cell2mat(cellstr(fullTDF(:, subjRespCol))) == conRespz(2);% response condition2

occCon1Idx = cell2mat(fullTDF(:, occConCol)) == 0;% occlusion condition 1 (non occluded)
occCon2Idx = cell2mat(fullTDF(:, occConCol)) == 1;% occlusion condition 2 (occluded)

% FOR ALL OCCLUSION CONDITIONS COMBINED
% Base image condition1 (conRespz(1)) and correct
both_logicalIndex1_1=logical(biCon1Idx.*rspCon1Idx); % condition1 resp1
% Base image condition1 (conRespz(1)) and incorrect
both_logicalIndex1_2=logical(biCon1Idx.*rspCon2Idx); % condition1 resp2
% Base image condition2 (conRespz(2)) and correct
both_logicalIndex2_1=logical(biCon2Idx.*rspCon2Idx); % condition2 resp2
% Base image condition2 (conRespz(2)) and incorrect
both_logicalIndex2_2=logical(biCon2Idx.*rspCon1Idx); % condition2 resp1

% FOR OCCLUDED ONLY
% Base image condition1 (conRespz(1)) and correct
occ_logicalIndex1_1=logical(occCon1Idx.*(biCon1Idx.*rspCon1Idx)); % condition1 resp1
% Base image condition1 (conRespz(1)) and incorrect
occ_logicalIndex1_2=logical(occCon1Idx.*(biCon1Idx.*rspCon2Idx)); % condition1 resp2
% Base image condition2 (conRespz(2)) and correct
occ_logicalIndex2_1=logical(occCon1Idx.*(biCon2Idx.*rspCon2Idx)); % condition2 resp2
% Base image condition2 (conRespz(2)) and incorrect
occ_logicalIndex2_2=logical(occCon1Idx.*(biCon2Idx.*rspCon1Idx)); % condition2 resp1

% FOR UN-OCCLUDED ONLY
% Base image condition1 (conRespz(1)) and correct
nocc_logicalIndex1_1=logical(occCon2Idx.*(biCon1Idx.*rspCon1Idx)); % condition1 resp1
% Base image condition1 (conRespz(1)) and incorrect
nocc_logicalIndex1_2=logical(occCon2Idx.*(biCon1Idx.*rspCon2Idx)); % condition1 resp2
% Base image condition2 (conRespz(2)) and correct
nocc_logicalIndex2_1=logical(occCon2Idx.*(biCon2Idx.*rspCon2Idx)); % condition2 resp2
% Base image condition2 (conRespz(2)) and incorrect
nocc_logicalIndex2_2=logical(occCon2Idx.*(biCon2Idx.*rspCon1Idx)); % condition2 resp1
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% Use the logical index variables to filter rows
%-------------------------------------------------------------------------
CI_11both_nz=fullTDF(both_logicalIndex1_1, noiseImgCol)'; % con1 resp correct
CI_12both_nz=fullTDF(both_logicalIndex1_2, noiseImgCol)'; % con2 resp incorrect
CI_21both_nz=fullTDF(both_logicalIndex2_1, noiseImgCol)'; % con2 resp correct
CI_22both_nz=fullTDF(both_logicalIndex2_2, noiseImgCol)'; % con2 resp incorrect

CI_11occ_nz=fullTDF(occ_logicalIndex1_1, noiseImgCol)'; % con1 resp correct
CI_12occ_nz=fullTDF(occ_logicalIndex1_2, noiseImgCol)'; % con2 resp incorrect
CI_21occ_nz=fullTDF(occ_logicalIndex2_1, noiseImgCol)'; % con2 resp correct
CI_22occ_nz=fullTDF(occ_logicalIndex2_2, noiseImgCol)'; % con2 resp incorrect

CI_11nocc_nz=fullTDF(nocc_logicalIndex1_1, noiseImgCol)'; % con1 resp correct
CI_12nocc_nz=fullTDF(nocc_logicalIndex1_2, noiseImgCol)'; % con2 resp incorrect
CI_21nocc_nz=fullTDF(nocc_logicalIndex2_1, noiseImgCol)'; % con2 resp correct
CI_22nocc_nz=fullTDF(nocc_logicalIndex2_2, noiseImgCol)'; % con2 resp incorrect

clear fullTDF; % clear the fullTDF array for RAM space concerns..

%-------------------------------------------------------------------------
% Transfer data from cell arrays to mats
%-------------------------------------------------------------------------

% FOR ALL OCCLUSION CONDITIONS COMBINED
%--------------------------------------------------------------------------
CI_11bothmat=zeros(size(CI_11both_nz{1,1},1),size(CI_11both_nz{1,1},2),length(CI_11both_nz));
for jj = 1:length(CI_11both_nz)
%     if bpFilter==1
%         CI_11bothmat(:,:,jj) = applyBandPassFilter(CI_11both_nz{1,jj}, @myBandPassFilter,lowFreq,highFreq);
%     else
        CI_11bothmat(:,:,jj)=CI_11both_nz{1,jj};
%     end
end
clear jj;
clear CI_11both_nz;

CI_21bothmat=zeros(size(CI_21both_nz{1,1},1),size(CI_21both_nz{1,1},2),length(CI_21both_nz));
for jj = 1:length(CI_21both_nz)
%     if bpFilter==1
%         CI_21bothmat(:,:,jj) = applyBandPassFilter(CI_21both_nz{1,jj}, @myBandPassFilter,lowFreq,highFreq);
%     else
        CI_21bothmat(:,:,jj)=CI_21both_nz{1,jj};
%     end
end
clear jj;
clear CI_21both_nz;


CI_12bothmat=zeros(size(CI_12both_nz{1,1},1),size(CI_12both_nz{1,1},2),length(CI_12both_nz));
for jj = 1:length(CI_12both_nz)
%     if bpFilter==1
%         CI_12bothmat(:,:,jj) = applyBandPassFilter(CI_12both_nz{1,jj}, @myBandPassFilter,lowFreq,highFreq);
%     else
        CI_12bothmat(:,:,jj)=CI_12both_nz{1,jj};
%     end
end
clear jj;
clear CI_12both_nz;

CI_22bothmat=zeros(size(CI_22both_nz{1,1},1),size(CI_22both_nz{1,1},2),length(CI_22both_nz));
for jj = 1:length(CI_22both_nz)
%     if bpFilter==1
%         CI_22bothmat(:,:,jj) = applyBandPassFilter(CI_22both_nz{1,jj}, @myBandPassFilter,lowFreq,highFreq);
%     else
        CI_22bothmat(:,:,jj)=CI_22both_nz{1,jj};
%     end
end
clear jj;
clear CI_22both_nz;
%--------------------------------------------------------------------------

% FOR OCCLUDED ONLY
%--------------------------------------------------------------------------
CI_11occmat=zeros(size(CI_11occ_nz{1,1},1),size(CI_11occ_nz{1,1},2),length(CI_11occ_nz));
for jj = 1:length(CI_11occ_nz)
%     if bpFilter==1
%         CI_11occmat(:,:,jj) = applyBandPassFilter(CI_11occ_nz{1,jj}, @myBandPassFilter,lowFreq,highFreq);
%     else
        CI_11occmat(:,:,jj)=CI_11occ_nz{1,jj};
%     end

end
clear jj;
clear CI_11occ_nz;

CI_21occmat=zeros(size(CI_21occ_nz{1,1},1),size(CI_21occ_nz{1,1},2),length(CI_21occ_nz));
for jj = 1:length(CI_21occ_nz)
%     if bpFilter==1
%         CI_21occmat(:,:,jj) = applyBandPassFilter(CI_21occ_nz{1,jj}, @myBandPassFilter,lowFreq,highFreq);
%     else
        CI_21occmat(:,:,jj)=CI_21occ_nz{1,jj};
%     end
end
clear jj;
clear CI_21occ_nz;

CI_12occmat=zeros(size(CI_12occ_nz{1,1},1),size(CI_12occ_nz{1,1},2),length(CI_12occ_nz));
for jj = 1:length(CI_12occ_nz)

%     if bpFilter==1
%         CI_12occmat(:,:,jj) = applyBandPassFilter(CI_12occ_nz{1,jj}, @myBandPassFilter,lowFreq,highFreq);
%     else
        CI_12occmat(:,:,jj)=CI_12occ_nz{1,jj};
%     end
end
clear jj;
clear CI_12occ_nz;

CI_22occmat=zeros(size(CI_22occ_nz{1,1},1),size(CI_22occ_nz{1,1},2),length(CI_22occ_nz));
for jj = 1:length(CI_22occ_nz)
%     if bpFilter==1
%         CI_22occmat(:,:,jj) = applyBandPassFilter(CI_22occ_nz{1,jj}, @myBandPassFilter,lowFreq,highFreq);
%     else
        CI_22occmat(:,:,jj)=CI_22occ_nz{1,jj};
%     end
end
clear jj;
clear CI_22occ_nz;
%--------------------------------------------------------------------------


% FOR UN-OCCLUDED ONLY
%--------------------------------------------------------------------------
CI_11noccmat=zeros(size(CI_11nocc_nz{1,1},1),size(CI_11nocc_nz{1,1},2),length(CI_11nocc_nz));
for jj = 1:length(CI_11nocc_nz)
%     if bpFilter==1
%         CI_11noccmat(:,:,jj) = applyBandPassFilter(CI_11nocc_nz{1,jj}, @myBandPassFilter,lowFreq,highFreq);
%     else
        CI_11noccmat(:,:,jj)=CI_11nocc_nz{1,jj};
%     end
end
clear jj;
clear CI_11nocc_nz;

CI_21noccmat=zeros(size(CI_21nocc_nz{1,1},1),size(CI_21nocc_nz{1,1},2),length(CI_21nocc_nz));
for jj = 1:length(CI_21nocc_nz)
%     if bpFilter==1
%         CI_21noccmat(:,:,jj) = applyBandPassFilter(CI_21nocc_nz{1,jj}, @myBandPassFilter,lowFreq,highFreq);
%     else
        CI_21noccmat(:,:,jj)=CI_21nocc_nz{1,jj};
%     end

end
clear jj;
clear CI_21nocc_nz;

CI_12noccmat=zeros(size(CI_12nocc_nz{1,1},1),size(CI_12nocc_nz{1,1},2),length(CI_12nocc_nz));
for jj = 1:length(CI_12nocc_nz)
%     if bpFilter==1
%         CI_12noccmat(:,:,jj) = applyBandPassFilter(CI_12nocc_nz{1,jj}, @myBandPassFilter,lowFreq,highFreq);
%     else
        CI_12noccmat(:,:,jj)=CI_12nocc_nz{1,jj};
%     end
end
clear jj;
clear CI_12nocc_nz;

CI_22noccmat=zeros(size(CI_22nocc_nz{1,1},1),size(CI_22nocc_nz{1,1},2),length(CI_22nocc_nz));
for jj = 1:length(CI_22nocc_nz)
%     if bpFilter==1
%         CI_22noccmat(:,:,jj) = applyBandPassFilter(CI_22nocc_nz{1,jj}, @myBandPassFilter,lowFreq,highFreq);
%     else
        CI_22noccmat(:,:,jj)=CI_22nocc_nz{1,jj};
%     end

end
clear jj;
clear CI_22nocc_nz;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% =========================================================================
tpoint22=datetime;
time2=tpoint22-tpoint21;
disp(strcat("Catagorizing noise images took ",num2str(seconds(time2)), " seconds..."));
disp(" ");

%% Apply specified flip to noise matrices from 1 of the 2 orientations if requested

if oriConFlp==1

    if FlpdOriCon==1
    CI_11bothmat=flip(CI_11bothmat,flpAxis);
    CI_21bothmat=flip(CI_21bothmat,flpAxis);

    CI_11occmat=flip(CI_11occmat,flpAxis);
    CI_21occmat=flip(CI_21occmat,flpAxis);

    CI_11noccmat=flip(CI_11noccmat,flpAxis);
    CI_21noccmat=flip(CI_21noccmat,flpAxis);
    end

    if FlpdOriCon==2
    CI_12bothmat=flip(CI_12bothmat,flpAxis);
    CI_22bothmat=flip(CI_22bothmat,flpAxis);

    CI_12occmat=flip(CI_12occmat,flpAxis);
    CI_22occmat=flip(CI_22occmat,flpAxis);

    CI_12noccmat=flip(CI_12noccmat,flpAxis);
    CI_22noccmat=flip(CI_22noccmat,flpAxis);
    end

end

%% --- Compute CIs --- %
disp(" ");
disp("Computing CIs...")
disp(" ");
tpoint31=datetime;

if bpFilter==0
% Both Occ and Nocc Combined
%==========================================================================
CI_results1both = (mean(CI_11bothmat,3) + mean(CI_21bothmat,3)) - (mean(CI_12bothmat,3) + mean(CI_22bothmat,3));
CI_results2both = (mean(CI_12bothmat,3) + mean(CI_22bothmat,3)) - (mean(CI_11bothmat,3) + mean(CI_21bothmat,3));
%make smoothed version..
CI_results1bothSM=imgaussfilt(CI_results1both,smthKrnlSigma,'FilterDomain',fDomain);
CI_results2bothSM=imgaussfilt(CI_results2both,smthKrnlSigma,'FilterDomain',fDomain);
%==========================================================================
% Clean up
clear CI_11bothmat;
clear CI_21bothmat;
clear CI_12bothmat;
clear CI_22bothmat;

% Occ Only
%==========================================================================
CI_results1occ = (mean(CI_11occmat,3) + mean(CI_21occmat,3)) - (mean(CI_12occmat,3) + mean(CI_22occmat,3));
CI_results2occ = (mean(CI_12occmat,3) + mean(CI_22occmat,3)) - (mean(CI_11occmat,3) + mean(CI_21occmat,3));
%make smoothed version..
CI_results1occSM=imgaussfilt(CI_results1occ,smthKrnlSigma,'FilterDomain',fDomain);
CI_results2occSM=imgaussfilt(CI_results2occ,smthKrnlSigma,'FilterDomain',fDomain);
%==========================================================================
% Clean up
clear CI_11occmat;
clear CI_21occmat;
clear CI_12occmat;
clear CI_22occmat;


% Nocc Only
%==========================================================================
CI_results1nocc = (mean(CI_11noccmat,3) + mean(CI_21noccmat,3)) - (mean(CI_12noccmat,3) + mean(CI_22noccmat,3));
CI_results2nocc = (mean(CI_12noccmat,3) + mean(CI_22noccmat,3)) - (mean(CI_11noccmat,3) + mean(CI_21noccmat,3));
%make smoothed version..
CI_results1noccSM=imgaussfilt(CI_results1nocc,smthKrnlSigma,'FilterDomain',fDomain);
CI_results2noccSM=imgaussfilt(CI_results2nocc,smthKrnlSigma,'FilterDomain',fDomain);
%==========================================================================
% Clean up
clear CI_11noccmat;
clear CI_21noccmat;
clear CI_12noccmat;
clear CI_22noccmat;
end


if bpFilter==1
% Both Occ and Nocc Combined
%==========================================================================
%CI_results1both = (applyBandPassFilter(mean(CI_11bothmat,3), @myBandPassFilter,lowFreq,highFreq) + applyBandPassFilter(mean(CI_21bothmat,3), @myBandPassFilter,lowFreq,highFreq)) - (applyBandPassFilter(mean(CI_12bothmat,3), @myBandPassFilter,lowFreq,highFreq) + applyBandPassFilter(mean(CI_22bothmat,3), @myBandPassFilter,lowFreq,highFreq));
%CI_results2both = (applyBandPassFilter(mean(CI_12bothmat,3), @myBandPassFilter,lowFreq,highFreq) + applyBandPassFilter(mean(CI_22bothmat,3), @myBandPassFilter,lowFreq,highFreq)) - (applyBandPassFilter(mean(CI_11bothmat,3), @myBandPassFilter,lowFreq,highFreq) + applyBandPassFilter(mean(CI_21bothmat,3), @myBandPassFilter,lowFreq,highFreq));

CI_results1both = (butterworthbpf_v2(mean(CI_11bothmat,3), lowerBnd,upperBnd,filterOrder,pltFig) + butterworthbpf_v2(mean(CI_21bothmat,3), lowerBnd,upperBnd,filterOrder,pltFig)) - (butterworthbpf_v2(mean(CI_12bothmat,3), lowerBnd,upperBnd,filterOrder,pltFig) + butterworthbpf_v2(mean(CI_22bothmat,3), lowerBnd,upperBnd,filterOrder,pltFig));
CI_results2both = (butterworthbpf_v2(mean(CI_12bothmat,3), lowerBnd,upperBnd,filterOrder,pltFig) + butterworthbpf_v2(mean(CI_22bothmat,3), lowerBnd,upperBnd,filterOrder,pltFig)) - (butterworthbpf_v2(mean(CI_11bothmat,3), lowerBnd,upperBnd,filterOrder,pltFig) + butterworthbpf_v2(mean(CI_21bothmat,3), lowerBnd,upperBnd,filterOrder,pltFig));

%make smoothed version..
CI_results1bothSM=imgaussfilt(CI_results1both,smthKrnlSigma,'FilterDomain',fDomain);
CI_results2bothSM=imgaussfilt(CI_results2both,smthKrnlSigma,'FilterDomain',fDomain);
%==========================================================================
% Clean up
clear CI_11bothmat;
clear CI_21bothmat;
clear CI_12bothmat;
clear CI_22bothmat;

% Occ Only
%==========================================================================
%CI_results1occ = (applyBandPassFilter(mean(CI_11occmat,3), @myBandPassFilter,lowFreq,highFreq) + applyBandPassFilter(mean(CI_21occmat,3), @myBandPassFilter,lowFreq,highFreq)) - (applyBandPassFilter(mean(CI_12occmat,3), @myBandPassFilter,lowFreq,highFreq) + applyBandPassFilter(mean(CI_22occmat,3), @myBandPassFilter,lowFreq,highFreq));
%CI_results2occ = (applyBandPassFilter(mean(CI_12occmat,3), @myBandPassFilter,lowFreq,highFreq) + applyBandPassFilter(mean(CI_22occmat,3), @myBandPassFilter,lowFreq,highFreq)) - (applyBandPassFilter(mean(CI_11occmat,3), @myBandPassFilter,lowFreq,highFreq) + applyBandPassFilter(mean(CI_21occmat,3), @myBandPassFilter,lowFreq,highFreq));
CI_results1occ = (butterworthbpf_v2(mean(CI_11occmat,3), lowerBnd,upperBnd,filterOrder,pltFig) + butterworthbpf_v2(mean(CI_21occmat,3), lowerBnd,upperBnd,filterOrder,pltFig)) - (butterworthbpf_v2(mean(CI_12occmat,3), lowerBnd,upperBnd,filterOrder,pltFig) + butterworthbpf_v2(mean(CI_22occmat,3), lowerBnd,upperBnd,filterOrder,pltFig));
CI_results2occ = (butterworthbpf_v2(mean(CI_12occmat,3), lowerBnd,upperBnd,filterOrder,pltFig) + butterworthbpf_v2(mean(CI_22occmat,3), lowerBnd,upperBnd,filterOrder,pltFig)) - (butterworthbpf_v2(mean(CI_11occmat,3), lowerBnd,upperBnd,filterOrder,pltFig) + butterworthbpf_v2(mean(CI_21occmat,3), lowerBnd,upperBnd,filterOrder,pltFig));

%make smoothed version..
CI_results1occSM=imgaussfilt(CI_results1occ,smthKrnlSigma,'FilterDomain',fDomain);
CI_results2occSM=imgaussfilt(CI_results2occ,smthKrnlSigma,'FilterDomain',fDomain);
%==========================================================================
% Clean up
clear CI_11occmat;
clear CI_21occmat;
clear CI_12occmat;
clear CI_22occmat;


% Nocc Only
%==========================================================================
%CI_results1nocc = (applyBandPassFilter(mean(CI_11noccmat,3), @myBandPassFilter,lowFreq,highFreq) + applyBandPassFilter(mean(CI_21noccmat,3), @myBandPassFilter,lowFreq,highFreq)) - (applyBandPassFilter(mean(CI_12noccmat,3), @myBandPassFilter,lowFreq,highFreq) + applyBandPassFilter(mean(CI_22noccmat,3), @myBandPassFilter,lowFreq,highFreq));
%CI_results2nocc = (applyBandPassFilter(mean(CI_12noccmat,3), @myBandPassFilter,lowFreq,highFreq) + applyBandPassFilter(mean(CI_22noccmat,3), @myBandPassFilter,lowFreq,highFreq)) - (applyBandPassFilter(mean(CI_11noccmat,3), @myBandPassFilter,lowFreq,highFreq) + applyBandPassFilter(mean(CI_21noccmat,3), @myBandPassFilter,lowFreq,highFreq));
CI_results1nocc = (butterworthbpf_v2(mean(CI_11noccmat,3), lowerBnd,upperBnd,filterOrder,pltFig) + butterworthbpf_v2(mean(CI_21noccmat,3), lowerBnd,upperBnd,filterOrder,pltFig)) - (butterworthbpf_v2(mean(CI_12noccmat,3), lowerBnd,upperBnd,filterOrder,pltFig) + butterworthbpf_v2(mean(CI_22noccmat,3), lowerBnd,upperBnd,filterOrder,pltFig));
CI_results2nocc = (butterworthbpf_v2(mean(CI_12noccmat,3), lowerBnd,upperBnd,filterOrder,pltFig) + butterworthbpf_v2(mean(CI_22noccmat,3), lowerBnd,upperBnd,filterOrder,pltFig)) - (butterworthbpf_v2(mean(CI_11noccmat,3), lowerBnd,upperBnd,filterOrder,pltFig) + butterworthbpf_v2(mean(CI_21noccmat,3), lowerBnd,upperBnd,filterOrder,pltFig));

%make smoothed version..
CI_results1noccSM=imgaussfilt(CI_results1nocc,smthKrnlSigma,'FilterDomain',fDomain);
CI_results2noccSM=imgaussfilt(CI_results2nocc,smthKrnlSigma,'FilterDomain',fDomain);
%==========================================================================
% Clean up
clear CI_11noccmat;
clear CI_21noccmat;
clear CI_12noccmat;
clear CI_22noccmat;
end

tpoint32=datetime;
time3=tpoint32-tpoint31;
disp(strcat("Computing CIs took ",num2str(seconds(time3)), " seconds..."));
disp(" ");

%% Make Thresholded Versions..

disp(" ");
disp("Making/Plotting Thresholded CIs...")
disp(" ");

tpoint41=datetime;


% noise type/configuration used in sim
whichNtype=noiseType; % either kernel noise: "krnlNz" or white noise: "white"


% specify ideal CI to plot..
whichCI="";% either rawCI1, rawCI2, smthCI1, or smthCI2
whichConfig=""; %specify noise configuration in struct you want to visualize
whichIdealCI=""; % specify which ideal CI version you want (either idealCI1 or idealCI2).. should probably match the ci version you choose for whichCI

% specify figure dimensions..
figDims= [1000,500]; % figure dimensions ([xdim,ydim])

% auto-set noiseConfigStr describing noise used in this configuration based
% on whichNtype specified above..
if whichNtype == "krnlNz"
noiseConfigStr=strcat("Kernel Noise");
elseif whichNtype == "white"
noiseConfigStr="White Noise Control";
else
noiseConfigStr="";
end

% compile all CIs into input structure
cisIn=struct; % initialize
if useSmthVer4clrd == 0
cisIn.CI_results1both=CI_results1both;
cisIn.CI_results2both=CI_results2both;
cisIn.CI_results1occ=CI_results1occ;
cisIn.CI_results2occ=CI_results2occ;
cisIn.CI_results1nocc=CI_results1nocc;
cisIn.CI_results2nocc=CI_results2nocc;
end
if useSmthVer4clrd == 1
cisIn.CI_results1both=CI_results1bothSM;
cisIn.CI_results2both=CI_results2bothSM;
cisIn.CI_results1occ=CI_results1occSM;
cisIn.CI_results2occ=CI_results2occSM;
cisIn.CI_results1nocc=CI_results1noccSM;
cisIn.CI_results2nocc=CI_results2noccSM;
end

% Call plotCIasOverlay_v2
plotCIasOverlay_v2(biIn1,biIn2,cisIn,ciThr, figDims, clrMap,alphVal,whichConfig,noiseConfigStr);

tpoint42=datetime;
time4=tpoint42-tpoint41;
disp(strcat("Making/Plotting Thresholded CIs took ",num2str(seconds(time4)), " seconds..."));
disp(" ");

%% Plot Non-Thresholded CIs

disp(" ");
disp("Plotting non-thresholded CIs...")
disp(" ");

tpoint51=datetime;

pltCIs=1;

if pltCIs==1

% PLOT RAW CIs
%==========================================================================
f = figure;
f.Position = [100 100 900 750];
subplot(2,3,1);
imagesc(CI_results1both);
colormap("gray");
axis equal
axis tight
colorbar
title("Occ and Nocc Combined CI1 Raw")

subplot(2,3,4);
imagesc(CI_results2both);
colormap("gray");
axis equal
axis tight
colorbar
title("Occ and Nocc Combined CI2 Raw")

subplot(2,3,2);
imagesc(CI_results1occ);
colormap("gray");
axis equal
axis tight
colorbar
title("Occ Only CI1 Raw")

subplot(2,3,5);
imagesc(CI_results2occ);
colormap("gray");
axis equal
axis tight
colorbar
title("Occ Only CI2 Raw")

subplot(2,3,3);
imagesc(CI_results1nocc);
colormap("gray");
axis equal
axis tight
colorbar
title("Nocc Only CI1 Raw")

subplot(2,3,6);
imagesc(CI_results2nocc);
colormap("gray");
axis equal
axis tight
colorbar
title("Nocc Only CI2 Raw")
%==========================================================================

% PLOT SMOOTHED CIs
%==========================================================================
f2 = figure;
f2.Position = [100 100 900 750];
subplot(2,3,1);
imagesc(CI_results1bothSM);
colormap("gray");
axis equal
axis tight
colorbar
title("Occ and Nocc Combined CI1 Smooth")

subplot(2,3,4);
imagesc(CI_results2bothSM);
colormap("gray");
axis equal
axis tight
colorbar
title("Occ and Nocc Combined CI2 Smooth")

subplot(2,3,2);
imagesc(CI_results1occSM);
colormap("gray");
axis equal
axis tight
colorbar
title("Occ Only CI1 Smooth")

subplot(2,3,5);
imagesc(CI_results2occSM);
colormap("gray");
axis equal
axis tight
colorbar
title("Occ Only CI2 Smooth")

subplot(2,3,3);
imagesc(CI_results1noccSM);
colormap("gray");
axis equal
axis tight
colorbar
title("Nocc Only CI1 Smooth")

subplot(2,3,6);
imagesc(CI_results2noccSM);
colormap("gray");
axis equal
axis tight
colorbar
title("Nocc Only CI2 Smooth")
end

tpoint52=datetime;
time5=tpoint52-tpoint51;
disp(strcat("Plotting non-thresholded CIs took ",num2str(seconds(time5)), " seconds..."));
disp(" ");

cd(startDir);
%==========================================================================
%end

