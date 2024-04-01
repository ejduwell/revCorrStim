%% Parameters..

outDir="/home/eduwell/Matlab_Sandbox/imgKernelNoise/exp2Test/individualFaceBIs/1000";
outName="";
%% Compute ssim(reconBI,actualBI) and ssim(computedCI,idealCI) for the various kernel noises and white noise
structOut = KrnlSmplNoise_v10_Demo_multi();

%% Save Data
strtDir=pwd; % save start position
cd(outDir); % go to output dir..
save(outName,structOut,'-mat'); % save struct..
cd(strtDir); % return to start location