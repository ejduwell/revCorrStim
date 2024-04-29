function qparsOut = getPrevQuestPars(out_dir,subj,n_int,pthrds,stimVer)
%--------------------------------------------------------------------------
% Input Variables/Pars:
%--------------------------------------------------------------------------
% out_dir : parent/highest level data output directory (path string)
% subject : 3 letter subject ID (string)
% n_int   : number of interleaved quests (integer)
% pthrds  : desired thresholds for respective interleaved quests 
%          (vector of n_int acuracy threshold values ranging from 0-1.. 
%           ie 0.5-->"50 percent accuracy", etc)
% stimVer : vector of 1-->N number strings specifying stimulus version
%           being run. These are interpreted/used as subdiectories to keep the data
%           particular stimulus versions, versions of that particular version,
%           versions of versions of versions..etc. separated hierarchically in the
%           output directory structure. First string in list is the highest
%           level, next is subdir w/in that, and so forth.
%--------------------------------------------------------------------------
% Output Variables
%--------------------------------------------------------------------------
% qparsOut             : struct variable containing one or more of the following fields:
% qparsOut.isPrevQuest : indicates whether or not a previous QUEST session 
%                        was detected for this stimulus version/subjectID 
%                        combo (1 or 0).
% If a previous QUEST session is detected (qparsOut.isPrevQuest=1), the
% following fields will also be included:
% qparsOut.pthrds      : threshold estimates from previous QUEST session
%                      ( solutions from output from psychometric curve fit 
%                       to data in that session at the input theshold values)
% qparsOut.sdevs       : corresponding sd estimates from QUESTs 
%                       (ptb quest model)
%--------------------------------------------------------------------------

strtDir=pwd; % save start path location..
qparsOut=struct(); % initialize output struct..

% Go into the subjects directory and load the .mat from running quest
% use the quest-generated psychometric function to pull desired thresholds
cd(out_dir)

if ~exist(strcat(pwd,"/",subj), 'dir')
    disp(" ");
    disp(strcat("No output directory found for subject ID ",subj,". Creating it."));
    disp(" ");
    mkdir(subj);
end
cd(subj)

for ii = 1:length(stimVer)
    if ~exist(strcat(pwd,"/",stimVer(ii)), 'dir')
        disp(" ");
        disp(strcat("No ",stimVer(ii), " output subdirectory found. Creating it."));
        disp(" ");
        mkdir(stimVer(ii));
    end
    cd(stimVer(ii));
end

% Get a list of the subdirectories..
folder = pwd;
subDs = GetSubDirsFirstLevelOnly(folder);

if size(subDs,1) == 0
    % if there are no subdirectories, there is no previous quest data.. in
    % this case set qparsOut.isPrevQuest equal to 0 and skip the rest of this function;
    qparsOut.isPrevQuest=0;
    disp(" ");
    disp("No previous QUEST session detected.");
    disp(" ");
else

    


    % get list of subdirectories and make list of previously completed quest output subdirectories

    %subDs = subDs';
    pattern = "_cmpltd";
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
        % if there are subdirectories, set qparsOut.isPrevQuest=1;
        qparsOut.isPrevQuest=1;
        qst_subDs = vertcat(qst_subDs,thisFolderName);
        numFilesProcessed = numFilesProcessed + 1;	% For fun, let's keep track of how many files we processed.
    end
    fprintf('We found %d folders with %s in the name.\n',numFilesProcessed, pattern);
    
    if numFilesProcessed == 0
       qparsOut.isPrevQuest=0; 
    end
    clear numFilesProcessed

    if qparsOut.isPrevQuest == 1
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
    data = load(matFiles{1,1});

    nFuncns = size(data.functions_out,2); % check number of psychometric function fits from quest in functions_out
    Funcn_thrMat = zeros(nFuncns,n_int); % preallocate space for threshold parameter values..
    par_mat = [2,n_int]; % preallocate for par_mat

    if nFuncns > 1
        fprintf('We found %d psychometric function fits from the previous QUEST session. \n',nFuncns);
        fprintf('Using the average threshold values across the %d to determine the QUEST starting parameters for this session.. \n', nFuncns);
    else
        fprintf('We found %d psychometric function fit from the previous QUEST session. \n',nFuncns);
        fprintf('Using the threshold values from this function to determine the QUEST starting parameters for this session.. \n');
    end


    for ii = 1:nFuncns
        fcn_tmp = data.functions_out{1,ii};
        for nn=1:n_int
            thrldStr=strcat("thr",num2str(nn));
            thrTmp=pthrds.(thrldStr);
            par_mat(1,nn) = thrTmp; % save percent-accuracy threshold target value in first row of par_mat(these are to label)..
            Funcn_thrMat(ii,nn) = NRfcnThr_Comp(thrTmp, fcn_tmp);
        end
    end

    for nn=1:n_int
        par_mat(2,nn) = mean(Funcn_thrMat(:,nn)); % save the corresponding computed mean accuracy value in second row of par_mat
    end

    % pull the sd estimate data too.
    sdevsRaw=zeros(size(data.sd_out_out,2),n_int);%preallocate
    sdevs=zeros(1,n_int);
    % first pull in raw values..
    for vv=1:size(sdevsRaw,1)
        sdevsRaw(vv,1:n_int)=data.sd_out_out{1,vv}(1:n_int);
    end
    % Then compute means (this is only really relevant if there are multiple
    % runs.. if not, this step will not do anything)
    for vv=1:n_int
        sdevs(1,vv)=mean(sdevsRaw(:,vv));
    end

    disp(strcat("Below are the threshold estimates of the ",num2str(n_int)," interleaved QUEST thresholds from the last session:"));

    for nn=1:n_int
        thrldStr=strcat("thr",num2str(nn));
        disp(strcat(thrldStr, " (",num2str(pthrds.(thrldStr)),"): ",num2str(par_mat(2,nn))));
    end

    % Pack things up into output structure..
    qparsOut.pthrds=par_mat; % add threshold estimates
    qparsOut.sdevs=sdevs; % add the sd estimates..
    end
end

% return to start path location before closeing up..
cd(strtDir);
end