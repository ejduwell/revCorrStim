%function RevCorrCIGen_v1()


%% Auto-detect and set directory/path params:
% -------------------------------------
% get matlab directory path by "which-ing" for this file..
% NOTE: key assumption: there is only one copy of this on the path.. (that
% should always be the case to avoid other conflicts..)
mlab_dir = fileparts(which('RevCorrCIGen_v1.m'));
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
subjID="EJD"; % 3 letter subject ID
noiseType="krnlNz"; % which noise (white or krnlNz)
stimType="lum"; % which grouping parameter type (lum,cr, or tex)
stimVer="both"; % which version (occ, nocc, or both)

conRespz=['m','z']; % response options for the two conditions (right and left angle here..)

dataDir=strcat(main_dir,"/data_master/",subjID,"/",noiseType,"/",stimType,"/",stimVer);

%% Read in all data from subject's tdfs

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
array1 = {
    'A', 1, 'Yes';
    'B', 2, 'No';
    'C', 1, 'Yes';
    'D', 3, 'No'
};
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
% =========================================================================
pause ="";
%% Catagorize the Noise Images Based on Condition and Correctness

% temporarily separate headers from data..
headerz=fullTDF(1,:);
fullTDF=fullTDF(2:end,:);

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
respCol=strrep(respCol,'NA', '0');
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

% Use the logical index variables to filter rows
%-------------------------------------------------------------------------
% FOR ALL OCCLUSION CONDITIONS COMBINED
% Full TDF Rows..
% CI_11both=fullTDF(both_logicalIndex1_1, :); % con1 resp correct
% CI_12both=fullTDF(both_logicalIndex1_2, :); % con2 resp incorrect
% CI_21both=fullTDF(both_logicalIndex2_1, :); % con2 resp correct
% CI_22both=fullTDF(both_logicalIndex2_2, :); % con2 resp incorrect
% Just the noise images
CI_11both_nz=fullTDF(both_logicalIndex1_1, noiseImgCol)'; % con1 resp correct
CI_12both_nz=fullTDF(both_logicalIndex1_2, noiseImgCol)'; % con2 resp incorrect
CI_21both_nz=fullTDF(both_logicalIndex2_1, noiseImgCol)'; % con2 resp correct
CI_22both_nz=fullTDF(both_logicalIndex2_2, noiseImgCol)'; % con2 resp incorrect


% FOR OCCLUDED ONLY
% Full TDF Rows..
% CI_11occ=fullTDF(occ_logicalIndex1_1, :); % con1 resp correct
% CI_12occ=fullTDF(occ_logicalIndex1_2, :); % con2 resp incorrect
% CI_21occ=fullTDF(occ_logicalIndex2_1, :); % con2 resp correct
% CI_22occ=fullTDF(occ_logicalIndex2_2, :); % con2 resp incorrect
% Just the noise images
CI_11occ_nz=fullTDF(occ_logicalIndex1_1, noiseImgCol)'; % con1 resp correct
CI_12occ_nz=fullTDF(occ_logicalIndex1_2, noiseImgCol)'; % con2 resp incorrect
CI_21occ_nz=fullTDF(occ_logicalIndex2_1, noiseImgCol)'; % con2 resp correct
CI_22occ_nz=fullTDF(occ_logicalIndex2_2, noiseImgCol)'; % con2 resp incorrect


% FOR UN-OCCLUDED ONLY
% CI_11nocc=fullTDF(nocc_logicalIndex1_1, :); % con1 resp correct
% CI_12nocc=fullTDF(nocc_logicalIndex1_2, :); % con2 resp incorrect
% CI_21nocc=fullTDF(nocc_logicalIndex2_1, :); % con2 resp correct
% CI_22nocc=fullTDF(nocc_logicalIndex2_2, :); % con2 resp incorrect
% Just the noise images
CI_11nocc_nz=fullTDF(nocc_logicalIndex1_1, noiseImgCol)'; % con1 resp correct
CI_12nocc_nz=fullTDF(nocc_logicalIndex1_2, noiseImgCol)'; % con2 resp incorrect
CI_21nocc_nz=fullTDF(nocc_logicalIndex2_1, noiseImgCol)'; % con2 resp correct
CI_22nocc_nz=fullTDF(nocc_logicalIndex2_2, noiseImgCol)'; % con2 resp incorrect
%-------------------------------------------------------------------------

% =========================================================================


% BOTH OCC AND NOCC COMBINED
%--------------------------------------------------------------------------
CI_11bothmat=zeros(size(CI_11both_nz{1,1},1),size(CI_11both_nz{1,1},2),length(CI_11both_nz));
for jj = 1:length(CI_11both_nz)
CI_11bothmat(:,:,jj)=CI_11both_nz{1,jj};
end
clear jj;

CI_21bothmat=zeros(size(CI_21both_nz{1,1},1),size(CI_21both_nz{1,1},2),length(CI_21both_nz));
for jj = 1:length(CI_21both_nz)
CI_21bothmat(:,:,jj)=CI_21both_nz{1,jj};
end
clear jj;

CI_12bothmat=zeros(size(CI_12both_nz{1,1},1),size(CI_12both_nz{1,1},2),length(CI_12both_nz));
for jj = 1:length(CI_12both_nz)
CI_12bothmat(:,:,jj)=CI_12both_nz{1,jj};
end
clear jj;

CI_22bothmat=zeros(size(CI_22both_nz{1,1},1),size(CI_22both_nz{1,1},2),length(CI_22both_nz));
for jj = 1:length(CI_22both_nz)
CI_22bothmat(:,:,jj)=CI_22both_nz{1,jj};
end
clear jj;
%--------------------------------------------------------------------------

% OCC Only
%--------------------------------------------------------------------------
CI_11occmat=zeros(size(CI_11occ_nz{1,1},1),size(CI_11occ_nz{1,1},2),length(CI_11occ_nz));
for jj = 1:length(CI_11occ_nz)
CI_11occmat(:,:,jj)=CI_11occ_nz{1,jj};
end
clear jj;

CI_21occmat=zeros(size(CI_21occ_nz{1,1},1),size(CI_21occ_nz{1,1},2),length(CI_21occ_nz));
for jj = 1:length(CI_21occ_nz)
CI_21occmat(:,:,jj)=CI_21occ_nz{1,jj};
end
clear jj;

CI_12occmat=zeros(size(CI_12occ_nz{1,1},1),size(CI_12occ_nz{1,1},2),length(CI_12occ_nz));
for jj = 1:length(CI_12occ_nz)
CI_12occmat(:,:,jj)=CI_12occ_nz{1,jj};
end
clear jj;

CI_22occmat=zeros(size(CI_22occ_nz{1,1},1),size(CI_22occ_nz{1,1},2),length(CI_22occ_nz));
for jj = 1:length(CI_22occ_nz)
CI_22occmat(:,:,jj)=CI_22occ_nz{1,jj};
end
clear jj;
%--------------------------------------------------------------------------

% NOCC Only
%--------------------------------------------------------------------------
CI_11noccmat=zeros(size(CI_11nocc_nz{1,1},1),size(CI_11nocc_nz{1,1},2),length(CI_11nocc_nz));
for jj = 1:length(CI_11nocc_nz)
CI_11noccmat(:,:,jj)=CI_11nocc_nz{1,jj};
end
clear jj;

CI_21noccmat=zeros(size(CI_21nocc_nz{1,1},1),size(CI_21nocc_nz{1,1},2),length(CI_21nocc_nz));
for jj = 1:length(CI_21nocc_nz)
CI_21noccmat(:,:,jj)=CI_21nocc_nz{1,jj};
end
clear jj;

CI_12noccmat=zeros(size(CI_12nocc_nz{1,1},1),size(CI_12nocc_nz{1,1},2),length(CI_12nocc_nz));
for jj = 1:length(CI_12nocc_nz)
CI_12noccmat(:,:,jj)=CI_12nocc_nz{1,jj};
end
clear jj;

CI_22noccmat=zeros(size(CI_22nocc_nz{1,1},1),size(CI_22nocc_nz{1,1},2),length(CI_22nocc_nz));
for jj = 1:length(CI_22nocc_nz)
CI_22noccmat(:,:,jj)=CI_22nocc_nz{1,jj};
end
clear jj;
%--------------------------------------------------------------------------

%% --- Compute CIs --- %

% Both Occ and Nocc Combined
%==========================================================================
CI_results1both = (mean(CI_11bothmat,3) + mean(CI_21bothmat,3)) - (mean(CI_12bothmat,3) + mean(CI_22bothmat,3));
CI_results2both = (mean(CI_12bothmat,3) + mean(CI_22bothmat,3)) - (mean(CI_11bothmat,3) + mean(CI_21bothmat,3));
%make smoothed version..
CI_results1bothSM=imgaussfilt(CI_results1both,2);
CI_results2bothSM=imgaussfilt(CI_results2both,2);
%==========================================================================

% Occ Only
%==========================================================================
CI_results1occ = (mean(CI_11occmat,3) + mean(CI_21occmat,3)) - (mean(CI_12occmat,3) + mean(CI_22occmat,3));
CI_results2occ = (mean(CI_12occmat,3) + mean(CI_22occmat,3)) - (mean(CI_11occmat,3) + mean(CI_21occmat,3));
%make smoothed version..
CI_results1occSM=imgaussfilt(CI_results1occ,2);
CI_results2occSM=imgaussfilt(CI_results2occ,2);
%==========================================================================

% Nocc Only
%==========================================================================
CI_results1nocc = (mean(CI_11noccmat,3) + mean(CI_21noccmat,3)) - (mean(CI_12noccmat,3) + mean(CI_22noccmat,3));
CI_results2nocc = (mean(CI_12noccmat,3) + mean(CI_22noccmat,3)) - (mean(CI_11noccmat,3) + mean(CI_21noccmat,3));
%make smoothed version..
CI_results1noccSM=imgaussfilt(CI_results1nocc,2);
CI_results2noccSM=imgaussfilt(CI_results2nocc,2);
%==========================================================================


%% Plot CIs
pltCIs=1;

if pltCIs==1

% PLOT RAW CIs
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

% PLOT SMOOTHED CIs
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


cd(startDir);
%end

