
%% Compute ssim(reconBI,actualBI) and ssim(computedCI,idealCI) for the various kernel noises and white noise
structOut = KrnlSmplNoise_v8_Demo_multi();
pause="";

%% Extract Data from StructOut and Load Into Arrays for Plotting Purposes

% get the list of fieldnames..
allFields=fieldnames(structOut); % all field names
configFields={}; % initialize for just the "config" fields
% load only fields beginning with "config" into configFields
for ii=1:length(allFields)
    fldTmp=allFields{ii,1};
    if strcmp(fldTmp(1:6),'config')==1
        configFields{end+1,1}=fldTmp;
    end
end

% first get the number of config fields
nConfigs=length(configFields);

%preallocate arrays
nKrnlReconCfgs=size(structOut.config1.reconBIs_ssim2BI_valz,2);
krnlNzArray=cell(8,nConfigs*nKrnlReconCfgs);
whiteNzArray=cell(8,nConfigs*nKrnlReconCfgs);

% loop through and load data from struct into krnlNzArray
countr=1; %initialize counter
for ii = 1:nConfigs

    cfgPass=configFields{ii,1}; % set the structure sub-field for this pass
    
    % set column indices data loaded on this pass..
    startIdx=countr;
    endIdx=(countr+nKrnlReconCfgs)-1;

    % Load in the data from this structure sub-field/config
    %----------------------------------------------------------------------
    % (1) recon kernel labels for overall ssim btw reconBIs made with each krnl and the real BI: structOut.(cfgPass).reconBIs_ssim2BI_lbl
    krnlNzArray(1,startIdx:endIdx)=num2cell(structOut.(cfgPass).reconBIs_ssim2BI_lbl);

    % (2) overall ssim btw reconBIs made with each krnl and the real BI: structOut.(cfgPass).reconBIs_ssim2BI_valz (lower row is krnlNz)
    krnlNzArray(2,startIdx:endIdx)=num2cell(structOut.(cfgPass).reconBIs_ssim2BI_valz(2,:));

    % (3) stdev of overall ssim btw reconBIs made with each krnl and the real BI: structOut.(cfgPass).reconBIs_ssim2BI_std (lower row is krnlNz)
    krnlNzArray(3,startIdx:endIdx)=num2cell(structOut.(cfgPass).reconBIs_ssim2BI_std(2,:));

    % (4) overall ssim btw raw CI1 generated with this noise and the "ideal" CI1: structOut.(cfgPass).CIs.krnlNz.rawCI1_ssim2Ideal1.ssimval
    krnlNzArray(4,startIdx:endIdx)={structOut.(cfgPass).CIs.krnlNz.rawCI1_ssim2Ideal1.ssimval};

    % (5) overall ssim btw raw CI2 generated with this noise and the "ideal" CI2: structOut.(cfgPass).CIs.krnlNz.rawCI2_ssim2Ideal2.ssimval
    krnlNzArray(5,startIdx:endIdx)={structOut.(cfgPass).CIs.krnlNz.rawCI2_ssim2Ideal2.ssimval};

    % (6) overall ssim btw smoothed CI1 generated with this noise and the "ideal" CI1: structOut.(cfgPass).CIs.krnlNz.smthCI1_ssim2Ideal1.ssimval
    krnlNzArray(6,startIdx:endIdx)={structOut.(cfgPass).CIs.krnlNz.smthCI1_ssim2Ideal1.ssimval};

    % (7) overall ssim btw smoothed CI2 generated with this noise and the "ideal" CI2: structOut.(cfgPass).CIs.krnlNz.smthCI2_ssim2Ideal2.ssimval
    krnlNzArray(7,startIdx:endIdx)={structOut.(cfgPass).CIs.krnlNz.smthCI2_ssim2Ideal2.ssimval};

    % (8) vector indicating kernels used in the construction of this kernel noise configuration: structOut.(cfgPass).krnlNz_krnlsUsed
    krnlNzArray(8,startIdx:endIdx)={num2str(structOut.(cfgPass).krnlNz_krnlsUsed)};
    %----------------------------------------------------------------------
    countr=countr+nKrnlReconCfgs; % update countr for next pass..
end

% loop through and load data from struct into whiteNzArray
countr=1; %initialize counter
for ii = 1:nConfigs

    cfgPass=configFields{ii,1}; % set the structure sub-field for this pass
    
    % set column indices data loaded on this pass..
    startIdx=countr;
    endIdx=(countr+nKrnlReconCfgs)-1;

    % Load in the data from this structure sub-field/config
    %----------------------------------------------------------------------
    % (1) recon kernel labels for overall ssim btw reconBIs made with each krnl and the real BI: structOut.(cfgPass).reconBIs_ssim2BI_lbl
    whiteNzArray(1,startIdx:endIdx)=num2cell(structOut.(cfgPass).reconBIs_ssim2BI_lbl);

    % (2) overall ssim btw reconBIs made with each krnl and the real BI: structOut.(cfgPass).reconBIs_ssim2BI_valz (upper row is whiteNz)
    whiteNzArray(2,startIdx:endIdx)=num2cell(structOut.(cfgPass).reconBIs_ssim2BI_valz(1,:));

    % (3) stdev of overall ssim btw reconBIs made with each krnl and the real BI: structOut.(cfgPass).reconBIs_ssim2BI_std (upper row is whiteNz)
    whiteNzArray(3,startIdx:endIdx)=num2cell(structOut.(cfgPass).reconBIs_ssim2BI_std(1,:));

    % (4) overall ssim btw raw CI1 generated with this noise and the "ideal" CI1: structOut.(cfgPass).CIs.white.rawCI1_ssim2Ideal1.ssimval
    whiteNzArray(4,startIdx:endIdx)={structOut.(cfgPass).CIs.white.rawCI1_ssim2Ideal1.ssimval};

    % (5) overall ssim btw raw CI2 generated with this noise and the "ideal" CI2: structOut.(cfgPass).CIs.white.rawCI2_ssim2Ideal2.ssimval
    whiteNzArray(5,startIdx:endIdx)={structOut.(cfgPass).CIs.white.rawCI2_ssim2Ideal2.ssimval};

    % (6) overall ssim btw smoothed CI1 generated with this noise and the "ideal" CI1: structOut.(cfgPass).CIs.white.smthCI1_ssim2Ideal1.ssimval
    whiteNzArray(6,startIdx:endIdx)={structOut.(cfgPass).CIs.white.smthCI1_ssim2Ideal1.ssimval};

    % (7) overall ssim btw smoothed CI2 generated with this noise and the "ideal" CI2: structOut.(cfgPass).CIs.white.smthCI2_ssim2Ideal2.ssimval
    whiteNzArray(7,startIdx:endIdx)={structOut.(cfgPass).CIs.white.smthCI2_ssim2Ideal2.ssimval};

    % (8) vector indicating kernels used in the construction of this kernel
    % noise configuration: structOut.(cfgPass).krnlNz_krnlsUsed ("NA" for
    % white noise as its not made out of kernels..)
    whiteNzArray(8,startIdx:endIdx)={"NA"};
    %----------------------------------------------------------------------
    countr=countr+nKrnlReconCfgs; % update countr for next pass..
end

%% Reorganize Data by BI-Recon Kernel for Plotting Ease.. 
krnlNzArrayByKrnl={}; % initialize krnlNzArrayByKrnl for storing reorganized verion of krnlNzArray
whiteNzArrayByKrnl={}; % for white noise too..
nRwz=size(krnlNzArray,1); % save number of rows per column in krnlNzArray

% First get the set of unique BI-Recon Kernel conditions
biKrnlConRow=krnlNzArray(1,:); % grab the top row (row with kernel condition labels)
unqBIrecon=unique(string(biKrnlConRow)); % get unique set.

for ii = 1:length(unqBIrecon)
    lblTemp=unqBIrecon(ii); % get label for this pass
    logIdxTmp=ismember(unqBIrecon, lblTemp); % get logical index set of columns of interest matching the label for this pass
    nTemp=sum(logIdxTmp); % count the number of instances
    
    tmpArray={}; % initialize array for this pass
    for yy=1:size(krnlNzArray,2)
        if strcmp(krnlNzArray{1,yy},lblTemp)==1
            tmpArray(1:nRwz,end+1)=krnlNzArray(:,yy);
        end
    end    
    krnlNzArrayByKrnl{1,end+1}=lblTemp;
    krnlNzArrayByKrnl{2,end}=tmpArray;
end

for ii = 1:length(unqBIrecon)
    lblTemp=unqBIrecon(ii); % get label for this pass
    logIdxTmp=ismember(unqBIrecon, lblTemp); % get logical index set of columns of interest matching the label for this pass
    nTemp=sum(logIdxTmp); % count the number of instances
    
    tmpArray={}; % initialize array for this pass
    for yy=1:size(whiteNzArray,2)
        if strcmp(whiteNzArray{1,yy},lblTemp)==1
            tmpArray(1:nRwz,end+1)=whiteNzArray(:,yy);
        end
    end    
    whiteNzArrayByKrnl{1,end+1}=lblTemp;
    whiteNzArrayByKrnl{2,end}=tmpArray;
end

%% Reorganize Kernel Noise by Kernels Used to Build the Noise too..
krnlNzArrayByNzKrnl={}; % initialize krnlNzArrayByKrnl for storing reorganized verion of krnlNzArray
nRwz=size(krnlNzArray,1); % save number of rows per column in krnlNzArray

% First get the set of unique BI-Recon Kernel conditions
nzKrnlConRow=krnlNzArray(8,:); % grab the bottom row (row with kernel noise condition labels)
unqNzKrnl=unique(string(nzKrnlConRow)); % get unique set.

for ii = 1:length(unqNzKrnl)
    lblTemp=unqNzKrnl(ii); % get label for this pass
    logIdxTmp=ismember(unqNzKrnl, lblTemp); % get logical index set of columns of interest matching the label for this pass
    nTemp=sum(logIdxTmp); % count the number of instances
    
    tmpArray={}; % initialize array for this pass
    for yy=1:size(krnlNzArray,2)
        if strcmp(krnlNzArray{8,yy},lblTemp)==1
            tmpArray(1:nRwz,end+1)=krnlNzArray(:,yy);
        end
    end    
    krnlNzArrayByNzKrnl{1,end+1}=lblTemp;
    krnlNzArrayByNzKrnl{2,end}=tmpArray;
end

%% Plot Data

% Specify Axis Range
xrng=[0 0.5];
yrng=[0 0.2];
figSizeX=700;
figSizeY=600;

% ALL DATA COMBINED FOR WHITE AND KERNEL NOISE for all recon kernels....
figure('Position', [10 10 figSizeX figSizeY])
%errorbar(cell2mat(krnlNzArray(6,:)),cell2mat(krnlNzArray(2,:)),cell2mat(krnlNzArray(3,:)),"vertical","o")
errorbar(cell2mat(krnlNzArray(2,:)),cell2mat(krnlNzArray(6,:)),cell2mat(krnlNzArray(3,:)),"horizontal","o",'CapSize',1)
hold on
%errorbar(cell2mat(whiteNzArray(6,:)),cell2mat(whiteNzArray(2,:)),cell2mat(whiteNzArray(3,:)),"vertical","o")
errorbar(cell2mat(whiteNzArray(2,:)),cell2mat(whiteNzArray(6,:)),cell2mat(whiteNzArray(3,:)),"horizontal","o",'CapSize',1)
xlabel("mean ssim between noise-reconBI and realBI")
ylabel("ssim between computed CI and ideal CI")
xlim(xrng);
ylim(yrng);
legend('Kernel Noise','White Noise')
title("FOR WHITE vs. KERNEL NOISE FOR ALL RECON KERNELS")
hold off

% FOR KERNEL NOISE BY RECON-BI KERNEL
figure('Position', [10 10 figSizeX figSizeY])
hold on

for ii=1:size(krnlNzArrayByKrnl,2)
arrayTmp=krnlNzArrayByKrnl{2,ii};
%errorbar(cell2mat(arrayTmp(6,:)),cell2mat(arrayTmp(2,:)),cell2mat(arrayTmp(3,:)),"vertical","o")
errorbar(cell2mat(arrayTmp(2,:)),cell2mat(arrayTmp(6,:)),cell2mat(arrayTmp(3,:)),"horizontal","o",'CapSize',1)
end
xlabel("mean ssim between noise-reconBI and realBI")
ylabel("ssim between computed CI and ideal CI")
xlim(xrng);
ylim(yrng);
title("FOR KERNEL NOISE BY RECON-BI KERNEL")
legend(string(krnlNzArrayByKrnl(1,:)))
hold off

% FOR WHITE NOISE BY RECON-BI KERNEL
figure('Position', [10 10 figSizeX figSizeY])
hold on
for ii=1:size(whiteNzArrayByKrnl,2)
arrayTmp=whiteNzArrayByKrnl{2,ii};
%errorbar(cell2mat(arrayTmp(6,:)),cell2mat(arrayTmp(2,:)),cell2mat(arrayTmp(3,:)),"vertical","o")
errorbar(cell2mat(arrayTmp(2,:)),cell2mat(arrayTmp(6,:)),cell2mat(arrayTmp(3,:)),"horizontal","o",'CapSize',1)
end
xlabel("mean ssim between noise-reconBI and realBI")
ylabel("ssim between computed CI and ideal CI")
title("FOR WHITE NOISE BY RECON-BI KERNEL")
legend(string(whiteNzArrayByKrnl(1,:)))
xlim(xrng);
ylim(yrng);
hold off

% FOR KERNEL NOISE BY KERNEL NOISE TYPE
figure('Position', [10 10 figSizeX figSizeY])
hold on
for ii=1:size(krnlNzArrayByNzKrnl,2)
arrayTmp=krnlNzArrayByNzKrnl{2,ii};
%errorbar(cell2mat(arrayTmp(6,:)),cell2mat(arrayTmp(2,:)),cell2mat(arrayTmp(3,:)),"vertical","o")
errorbar(cell2mat(arrayTmp(2,:)),cell2mat(arrayTmp(6,:)),cell2mat(arrayTmp(3,:)),"horizontal","o",'CapSize',1)
end
xlabel("mean ssim between noise-reconBI and realBI")
ylabel("ssim between computed CI and ideal CI")
title("FOR KERNEL NOISE BY KERNEL NOISE TYPE")
legend(string(krnlNzArrayByNzKrnl(1,:)))
xlim(xrng);
ylim(yrng);
hold off


% FOR THE CIS + ssim(ci,idealci) + ssim(reconBI,realBI) (Kernel noise)
figure('Position', [10 10 5*figSizeX 5*figSizeY]);
hold on
for ii = 1:size(configFields,1)

    subplot(9,size(configFields,1),ii);
    imagesc(structOut.(configFields{ii,1}).CIs.krnlNz.smthCI1);
    colormap("gray");
    axis equal;
    axis tight;   
    colorbar;
    title(strcat(num2str(structOut.(configFields{ii,1}).krnlNz_krnlsUsed)," (CI)"));

    subplot(9,size(configFields,1),ii+size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).CIs.krnlNz.smthCI1_ssim2Ideal1.ssimmap);
    title("ssim(smthCI1,idealCI1)");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    subplot(9,size(configFields,1),ii+2*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.krnlNz.ssim2BI_allKsComb);
    title("ReconBIssim2BI-allKs");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    % NOTE!!! : EJD NOTICED THAT THE WHITE AND TAILORED NOISE DATA FOR
    % ssim2BI_Kimz were flip-flopped in the struct! NEED TO SWITCH! IN THE
    % SCRIPT THAT CREATES THESE
    % SIMPLY REFERENCING THE "WHITE" FIELD FOR TAILORED FOR THE MOMENT!
    ntyp="krnlNz";
    subplot(9,size(configFields,1),ii+3*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.(ntyp).ssim2BI_Kimz(:,:,1));
    title("ReconBIssim2BI-K1");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    subplot(9,size(configFields,1),ii+4*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.(ntyp).ssim2BI_Kimz(:,:,2));
    title("ReconBIssim2BI-K2");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    subplot(9,size(configFields,1),ii+5*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.(ntyp).ssim2BI_Kimz(:,:,3));
    title("ReconBIssim2BI-K3");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    subplot(9,size(configFields,1),ii+6*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.(ntyp).ssim2BI_Kimz(:,:,4));
    title("ReconBIssim2BI-K4");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    subplot(9,size(configFields,1),ii+7*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.(ntyp).ssim2BI_Kimz(:,:,5));
    title("ReconBIssim2BI-K5");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    subplot(9,size(configFields,1),ii+8*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.white.ssim2BI_Kimz(:,:,6));
    title("ReconBIssim2BI-K6");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

end
hold off

% FOR THE CIS + ssim(ci,idealci) + ssim(reconBI,realBI) (White)
figure('Position', [10 10 5*figSizeX 5*figSizeY]);
hold on
for ii = 1:size(configFields,1)

    subplot(9,size(configFields,1),ii);
    imagesc(structOut.(configFields{ii,1}).CIs.white.smthCI1);
    colormap("gray");
    axis equal;
    axis tight;   
    colorbar;
    title(strcat("White", "(CI)"));

    subplot(9,size(configFields,1),ii+size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).CIs.white.smthCI1_ssim2Ideal1.ssimmap);
    title("ssim(smthCI1,idealCI1)");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    subplot(9,size(configFields,1),ii+2*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.white.ssim2BI_allKsComb);
    title("ReconBIssim2BI-allKs");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    % NOTE!!! : EJD NOTICED THAT THE WHITE AND TAILORED NOISE DATA FOR
    % ssim2BI_Kimz were flip-flopped in the struct! NEED TO SWITCH! IN THE
    % SCRIPT THAT CREATES THESE
    % SIMPLY REFERENCING THE "WHITE" FIELD FOR TAILORED FOR THE MOMENT!
    ntyp="white";
    subplot(9,size(configFields,1),ii+3*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.(ntyp).ssim2BI_Kimz(:,:,1));
    title("ReconBIssim2BI-K1");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    subplot(9,size(configFields,1),ii+4*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.(ntyp).ssim2BI_Kimz(:,:,2));
    title("ReconBIssim2BI-K2");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    subplot(9,size(configFields,1),ii+5*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.(ntyp).ssim2BI_Kimz(:,:,3));
    title("ReconBIssim2BI-K3");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    subplot(9,size(configFields,1),ii+6*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.(ntyp).ssim2BI_Kimz(:,:,4));
    title("ReconBIssim2BI-K4");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    subplot(9,size(configFields,1),ii+7*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.(ntyp).ssim2BI_Kimz(:,:,5));
    title("ReconBIssim2BI-K5");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    subplot(9,size(configFields,1),ii+8*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.white.ssim2BI_Kimz(:,:,6));
    title("ReconBIssim2BI-K6");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

end
hold off



% FOR THE BIS and RECON BIs(Kernel noise)
figure('Position', [10 10 5*figSizeX 5*figSizeY]);
hold on
for ii = 1:size(configFields,1)

    subplot(9,size(configFields,1),ii);
    imagesc(structOut.(configFields{ii,1}).CIs.krnlNz.smthCI1);
    colormap("gray");
    axis equal;
    axis tight;   
    colorbar;
    title(strcat(num2str(structOut.(configFields{ii,1}).krnlNz_krnlsUsed)," (CI)"));

    subplot(9,size(configFields,1),ii+size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).CIs.krnlNz.smthCI1_ssim2Ideal1.ssimmap);
    title("ssim(smthCI1,idealCI1)");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    subplot(9,size(configFields,1),ii+2*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.krnlNz.reconBI_allKsComb);
    title("ReconBI-allKs");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    % NOTE!!! : EJD NOTICED THAT THE WHITE AND TAILORED NOISE DATA FOR
    % ssim2BI_Kimz were flip-flopped in the struct! NEED TO SWITCH! IN THE
    % SCRIPT THAT CREATES THESE
    % SIMPLY REFERENCING THE "WHITE" FIELD FOR TAILORED FOR THE MOMENT!
    ntyp="krnlNz";
    subplot(9,size(configFields,1),ii+3*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.(ntyp).reconBI_KImz(:,:,1));
    title("ReconBI-K1");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    subplot(9,size(configFields,1),ii+4*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.(ntyp).reconBI_KImz(:,:,2));
    title("ReconBI-K2");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    subplot(9,size(configFields,1),ii+5*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.(ntyp).reconBI_KImz(:,:,3));
    title("ReconBI-K3");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    subplot(9,size(configFields,1),ii+6*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.(ntyp).reconBI_KImz(:,:,4));
    title("ReconBI-K4");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    subplot(9,size(configFields,1),ii+7*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.(ntyp).reconBI_KImz(:,:,5));
    title("ReconBI-K5");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    subplot(9,size(configFields,1),ii+8*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.white.reconBI_KImz(:,:,6));
    title("ReconBI-K6");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

end
hold off


% FOR THE BIS and RECON BIs(white noise)
figure('Position', [10 10 5*figSizeX 5*figSizeY]);
hold on
for ii = 1:size(configFields,1)

    subplot(9,size(configFields,1),ii);
    imagesc(structOut.(configFields{ii,1}).CIs.white.smthCI1);
    colormap("gray");
    axis equal;
    axis tight;   
    colorbar;
    title(strcat("White", "(CI)"));

    subplot(9,size(configFields,1),ii+size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).CIs.white.smthCI1_ssim2Ideal1.ssimmap);
    title("ssim(smthCI1,idealCI1)");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    subplot(9,size(configFields,1),ii+2*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.white.reconBI_allKsComb);
    title("ReconBI-allKs");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    % NOTE!!! : EJD NOTICED THAT THE WHITE AND TAILORED NOISE DATA FOR
    % ssim2BI_Kimz were flip-flopped in the struct! NEED TO SWITCH! IN THE
    % SCRIPT THAT CREATES THESE
    % SIMPLY REFERENCING THE "WHITE" FIELD FOR TAILORED FOR THE MOMENT!
    ntyp="white";
    subplot(9,size(configFields,1),ii+3*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.(ntyp).reconBI_KImz(:,:,1));
    title("ReconBI-K1");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    subplot(9,size(configFields,1),ii+4*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.(ntyp).reconBI_KImz(:,:,2));
    title("ReconBI-K2");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    subplot(9,size(configFields,1),ii+5*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.(ntyp).reconBI_KImz(:,:,3));
    title("ReconBI-K3");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    subplot(9,size(configFields,1),ii+6*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.(ntyp).reconBI_KImz(:,:,4));
    title("ReconBI-K4");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    subplot(9,size(configFields,1),ii+7*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.(ntyp).reconBI_KImz(:,:,5));
    title("ReconBI-K5");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

    subplot(9,size(configFields,1),ii+8*size(configFields,1));
    imagesc(structOut.(configFields{ii,1}).reconBIs.(ntyp).reconBI_KImz(:,:,6));
    title("ReconBI-K6");
    colormap("gray");
    axis equal;
    axis tight;  
    colorbar;

end
hold off