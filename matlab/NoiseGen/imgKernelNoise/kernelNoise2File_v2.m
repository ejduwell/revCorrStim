function kernelNoise2File_v2()

%% Parameters

%imgsPerKrnl=[1,1,1,1,1,1,0,0,0]; % adjusts the total number of subtile images from each respective kernel incorporated into the noise
% imgsPerKrnl_set = {...
%     [1,1,1,1,1,1,0,0,0],...
%     [6,0,0,0,0,0,0,0,0],...
%     [0,6,0,0,0,0,0,0,0],...
%     [0,0,6,0,0,0,0,0,0],...
%     [0,0,0,6,0,0,0,0,0],...
%     [0,0,0,0,6,0,0,0,0],...
%     [0,0,0,0,0,6,0,0,0]...
%     };

imgsPerKrnl_set = {...
    [1,1,1,1,1,1,1,0,0],...
    };

% imgsPerKrnl_set = {...
%     [1,2,3,4,5,6,7,0,0],...
%     };

Repz=10; % sets number of image reps/trials used for each pass..
%Repz=100; % sets number of image reps/trials used for each pass..
jobChunks=10; % sets the number of chunks you want to break the job 
              % (total reps) into..

krnlImgWgts=[1,1,1,1,1,1,1,1,1]; % adjusts relative weights applied to subtile images from each respective kernel
smoothTiles=[0,0,0,0,0,0,0,0,0]; % apply smoothing kernel to intermediate subtile images?
gSmthK1sz=[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]; % size of gaussian smoothing kernel for subtiles..
smoothFinal=0; % apply smoothing kernel to final image?
gSmthK2sz=0.5; % size of gaussian smoothing kernel for final image..
BI_WT=0.35;
N_WT=1-BI_WT;
nReps=Repz/jobChunks; % number of noise image copies/reps you want to run per pass..
nWorkers=8; % number of cpu nodes used in the parallelized portions..

%Specify desired size?
selectSize=1; % if 1, this means that we will use/resize the input image the selected size below
desiredSize=[512,512]; % (MUST BE SQUARE IMAGE AND DIMZ MUST BE POWER OF 2)

%LumOnly
%krnlNzBI="/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/images/test4_LumOnly-05-Mar-2024-10-12-48/nocc/L/BaseIm_occ_1_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3214_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_0_lum2_255_LWT_1_Lal_1_Ocon_1.png";

%TexOnly
%krnlNzBI="/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/images/test5_TexOnly-12-Mar-2024-18-48-17/nocc/L/BaseIm_occ_1_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3315_T_1checkrz_Tr1_1.075_Tr2_0.93023_Tal_1_L_0_lum1_127_lum2_127_LWT_1_Lal_1_Ocon_1.png";
%krnlNzBI="/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/images/test6_TexOnly-19-Apr-2024-11-43-14/nocc/L/BaseIm_occ_1_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3315_T_1checkrz_Tr1_2.5_Tr2_0.4_Tal_1_L_0_lum1_127_lum2_127_LWT_1_Lal_1_Ocon_1.png";

krnlNzBI="/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/images/test6_TexOnly-19-Apr-2024-11-43-14/nocc/L/BaseIm_occ_1_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3315_T_1checkrz_Tr1_2.5_Tr2_0.4_Tal_1_L_0_lum1_127_lum2_127_LWT_1_Lal_1_Ocon_1.png";

%krnlNzBI="/home/eduwell/SynologyDrive/SNAP/projects/imgKernelNoise/baseImgs/baboon512.png";

%CRonly w/Lum
%krnlNzBI="/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/images/test5_CROnly-11-Mar-2024-14-16-09/nocc/L/BaseIm_occ_1_ori_z__CR_1_CRpw_4_CRpl_0_CRal_1_CRob_1_CRdm_3315_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_51_lum2_51_LWT_1_Lal_1_Ocon_2.png";

% Out Directory Base (parent directory where you want your output dirs
% created)
outDirMain="/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/noise";

% Base for output directories created..
%baseOutDirName="lumOnlyBI_20000frms_krnlNz_imzPerKrnl";
%baseOutDirName="texOnlyBI_20000frms_krnlNz_imzPerKrnl";
%baseOutDirName="crOnlyBI_20000frms_krnlNz_imzPerKrnl";

%baseOutDirName="texOnlyBI_v2_20000frms_krnlNz_imzPerKrnl";
baseOutDirName="baboonTest";

%parProfile='Threads';
%parProfile='Processes';
parProfile='HPC Cluster';
ClstrTimeStr='--time=00-02:00:00';
mlabComDirTmp="/home/eduwell/Matlab_Sandbox/imgKernelNoise/tmp";

%% Call Function

% read in krnlNzBI
krnlNzBI=imread(krnlNzBI);
krnlNzBI = imgReformater(krnlNzBI,desiredSize);
% Go to output parent directory...
strtDir=pwd; % save start location..
cd(outDirMain);

% Set up local parallel computing stuff..
%--------------------------------------------------------------------------
parpool("Processes",nWorkers);
%--------------------------------------------------------------------------
for mm=1:length(imgsPerKrnl_set)
imgsPerKrnl=imgsPerKrnl_set{1,mm}; % set imgsPerKrnl to specify kernels used to build the tailored noise for this pass...

% Build Out Directory Name..
imgsPerKrnlStr=string(imgsPerKrnl);
imPerK_str="_"; % initialize
for xx=1:length(imgsPerKrnlStr)
    strTmp=imgsPerKrnlStr(xx);
    imPerK_str=strcat(imPerK_str,strTmp);
end
OutDirTmp=strcat(baseOutDirName,imPerK_str);


frmCountr=1; % initialize frame countr
for qq=1:jobChunks
% Generate the noise images
[nzImz] = kernelNoise_v1_1(imgsPerKrnl,desiredSize,krnlNzBI,nReps,krnlImgWgts,smoothTiles,gSmthK1sz,smoothFinal,gSmthK2sz);

if ~exist(OutDirTmp, 'dir')
    mkdir(OutDirTmp) % make the output directory..
end
cd(OutDirTmp); % enter it

% write out the noise files
nImzTemp=size(nzImz,3);
for ii=1:nImzTemp
    imNameStr=strcat("noiseSample",num2str(frmCountr,'%05.f'),".png");
    frmCountr=frmCountr+1; % update frame countr for next pass..
    imwrite(nzImz(:,:,ii),imNameStr);
end
cd .. %go back up a level..

end
end

% Close down local parpool
delete(gcp('nocreate'));

% return to start location
cd(strtDir);

end