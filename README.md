# revCorrStim

<p align="center">

</p>

<p align="center">
  <img src="https://github.com/ejduwell/revCorrStim/blob/main/images/githubExamples/asciiSnapshot.png" alt="revCorrStim">
  <br>
</p>
<p align="center">

# Overview:
<p align="center">
  <div style="text-align: center;">

  RevCorrStim is a repository of Matlab code for optimizing and running reverse correlation experiments aimed specifically at examining various object grouping principles.
    
  This repositiory uses a combination of Psychtoolbox functions and independently developed code to address the four most basic aspects of the psychophysical reverse correlation method:

  1) Generating base images
  2) Generating noise images
  3) Stimulus presentation / task design
  4) Analysis software for generating classification and significance images (CIs and SIs)
  </div>
</p>

# Requirements/Dependencies:
- Must have a relatively recent installation of MATLAB (2018 or newer)
- Must have Psychtoolbox3 installed ( [link to PTB3 install](http://psychtoolbox.org/download.html) )
- Must be running either MacOS or a Linux-based operating system such as Ubuntu, Fedora, CentOs, etc.
- Must have git package installed to install via "git clone" method


# Installation:

1. Open a terminal window on your local machine and navigate to desired install location

2. Clone the revCorrStim repository to that location by running:

   ```
   git clone https://github.com/ejduwell/revCorrStim.git
   ```
   Alternatively, you can also download the package by hitting the green "<> Code" button above and then hitting the "Download ZIP" button to interactively download the package through your browser (clunkier, but works too..). You'll then need to unzip the folder in your desired location.

4. Add revCorrStim folder to your Matlab path. This can be done one of two ways:
   
   a. by manually right-clicking on the folder each time you open Matlab and selecting "Add to Path > Selected Folders and Subfolders"
   
   b. by permenantly adding the folder to your path by adding a line to your startup.m file such that this operation is done automatically each time 
 Matlab starts. (startup.m should be in '/home/{YourUserName}/Documents/MATLAB'.. if not present create it. then add line below ):

   ```
   % add revCorrStim to path
   addpath(genpath("full/path/to/your/revCorrStim/folder"));
   ```

    
# Usage:

This section explains how to use various aspects of the revCorrStim repository. Namely how to:
- Generate Base Images ( [jump to section](https://github.com/ejduwell/revCorrStim/blob/main/README.md#base-image-generation) )
- Generate Noise Images ( [jump to section](https://github.com/ejduwell/revCorrStim/blob/main/README.md#noise-image-generation) )
- Run the Stimulus ( [jump to section](https://github.com/ejduwell/revCorrStim/blob/main/README.md#stimulus-presentation) )
- Generate CIs/SIs ( [jump to section](https://github.com/ejduwell/revCorrStim/blob/main/README.md#ci--si-generation) )

In general, each of these aspects of the revCorrStim repository were intentionally created to stand alone as seperate "modules."

Base images are pre-generated and saved in an output subdirectory located within '/revCorrStim/images'. 
You'll notice that the auto-generated base images have obnoxiously long names. This is because they are named systematically such that all pertinent parameter values are encoded in the file name. These parameters are then read out when they are loaded by the stimulus presentation program.

Similarly, noise images are also pre-generated and saved in a subdirectory within /revCorrStim/noise. 

This approach allows users to create multiple different versions of an experiment with various base image / noise combinations, but run them all with the same centralized stimulus presentation program and analyze the output data with the same analysis software.

<p align="center">
  <div style="text-align: center;">
  **NOTE: Consequently, before you can run the stimulus your first time, you will first need to:
    
  1) generate base & noise images
  2) update a number of parameters to point the stimulus presentation program to these images. 
  
  Also, perhaps it goes without saying, but before you can generate CIs/SIs, you must first run an experiment and have data to analyze..
  </div>
</p>


## Base Image Generation

The Matlab code for generating base images is located in the '/revCorrStim/matlab/BaseImgGen' subdirectory.

The most recent copy of the 'main function' for generating base images is GenRevCorrBaseImsRF_v4.m

Detailed instructions on how to use GenRevCorrBaseImsRF_v4.m are provide below, but very briefly, here is what GenRevCorrBaseImsRF_v4.m is / what it does:
- GenRevCorrBaseImsRF is short for "Generate Reverse Correlation Base Images Recursively"
- Simply put, it pre-generates a set of base images which will later be read-in by the stimulus presentation program.
- However, it does this in a powerfully generalized manner:
    - The user specifies a set of vectors in the parameter section.
    - Each vector defines the desired set of possible values for a parameter that controls a particular aspect of the base image.
    - The user also specifies a list of logical constraints which limit which parameter combinations will be allowed in the output image set.
    - The full set of all possible parameter combinations is generated by feeding this list of vectors into "combvec"
    - The set of all possible parameter combos is then reduced via logical indexing to only those which meet the specified logical constraints.
    - This set of parameter combos is then used to generate the corresponding base images.
- Effectively: GenRevCorrBaseImsRF_v4.m allows users to pregenerate a set of base images which tile an entire n-dimensional parameter space in a highly flexible manner.

To generate base images you will need to update a few parameters in GenRevCorrBaseImsRF_v4.m and then run it.

Follow the instructions below:
1) Open GenRevCorrBaseImsRF_v4.m in your Matlab editor window.
2) Scroll down to the '%% Parameters' section
3) Under the '%% System Parameters' subsection, change the 'outdir' to the name you want to assign to the output parent directory.
4) Select which grouping parameters you want to use:
     - revCorrStim currently supports 3 different "object grouping parameters" out of the box (common region, texture, and luminance)
     - each of these hase their own parameter subsection ('%% Common Region Parameters', '%% Texture Parameters', and '%% Luminance Parameters')
     - to turn these on/off set the first parameter (T, CR, or L respectively) equal to either 1 or 0 (on or off).
     - You can, in fact, select any combination of these that you want. However, for the sake of simplicity, lets stick to luminance only for this example (set L=1, T=0, and CR=0)
5) Set the range of possibile parameter values for object 1 and object 2:
     - In this case, because we're just manipulating luminance, that means we will be adjusting 'lum1 and lum2' (in '%% Luminance Parameters' section)
     - Adjust the linspace command to specify the vectors of possible values for each. For example:
       ```
       lum1 = linspace(0,127.5,127);  
       lum2 = linspace(127.5,255,127); % num between 0-255 150
       ```
       This will make lum1 range from 0-127.5 in 127 steps and lum2 range from 127.5-255 in 127 steps.
6) Finally set up the parameters whose value vectors you want included in the set 'mega set' of all possible value combinations.
   - this is done in the '%% Parameters to Recursively Iterate Through Possible Combinations' subsection
   - Scroll down to the "if statement" matching your use case of interest. In our case, for this example its:
   ```
   % LUM ONLY NARROW RANGE, MEAN OBJECT LUM == BACKGROUND
    if L==1 && T==0 && CR==0
    itrPars   = {occluders,oriCon,objcon,lum1,lum2};           % list names of parameter variables you want to iterate
    itParNamz = {"occluders","oriCon","objcon","lum1","lum2"}; % list names of parameter variables you want to iterate (in itrParz above) as strings
    Eqlz = {"lum1~=lum2","((lum1+lum2)/2)==127.5"}; % specify equality conditions you want
    itrParsCom = combvec(itrPars{1},itrPars{2},itrPars{3},itrPars{4},itrPars{5}); % This builds a matrix of all the combos of pars in vectors listed     in itrPars.
    end
   ```
   - The current setup specifies (as noted in the comment header), a configuration that will generate images where lum1 and lum2 vary about the background luminance of 127.5 while maintaining a mean luminance of 127.5. We will not adjust the values in this example, but for reference here is what each of the parameters do:
       - itrPars -- specfies the set of parameter vectors (defined in earlier parameter sections) whose values will be included in generating the "set of all possible value combinations"
       - itParNamz -- specifies a set of 'name tag' strings corresponding to the respective parameters in itrPars. (These strings are used in the systematic file naming)
       - Eqlz -- specifies a list of strings containing equality statements. These logical statements basically allow you to specify conditions/cases you want included in the output image set. In this case we speficied that we only want conditions where lum1 and lum2 are different values/not equal ("lum1~=lum2") and conditions where the mean of lum1 and lum2 is equal to 127.5 ("((lum1+lum2)/2)==127.5"). These logical statements are used in a subsequent logical indexing step to constrain the set of all possible parameter combos to the specific subset which meet all the specified conditions.
         
         NOTE: Think carefully when setting this parameter in other scenarios! Without proper constraints, the number of possible image combos can quickly balloon into an enormous, unworkable number of images, take a huge amount of time to generate, and lots of disk space. However if used judicioulsy/creatively this parameter is extremely powerful! It allows you to quickly and intuitively create sets of images with any imaginable parameter combos and constraints all with a single strategy which essentially amounts to: "generate the set of all possible parameter combos, but only create images for the subset that meet these conditions.." Generating the set of all possible parameter combinations may seem extreme, but it is actually really fast/easy to compute using "combvec" (even for HUGE combo sets). Drawing/saving the images is the slow step. 

       - itrParsCom -- stores thes set of all possible combinations of the parameter values for the parameters in itrPars. itrParsCom is computed using the "combvec" command in this one-liner call. You feed it a list of vectors, it spits out the set of all possible combinations of values within those vectors (and mind-bogglingly fast..). Most importantly: If you make a new configuration with a different number of parameter vectors in itrPars, you will simply need to add or remove values such that the list of vectors in 'combvec(itrPars{1},itrPars{2},itrPars{3}, ... itrPars{n})' mirrors the number of parameter vectors (n) in itrPars.
     
    - You'll notice that the other configurations in the '%% Parameters to Recursively Iterate Through Possible Combinations' subsection all define the same 4 parameters explained above under their respective "if statements". All that changes in each different configuration is the list of particular parameter vectors included in the itrPars set, the corresponding strings in itParNamz, and the equality conditions specified. Once understand how one configuration works, you should be able to understand the rest / be able make your own new configurations if you wish.

7) Finally, run GenRevCorrBaseImsRF_v4 either by hitting the green "Run" button at the top of the matlab editor window or by running 'GenRevCorrBaseImsRF_v4' in the command window.

   Note: This may take some time to run. It will also take over a portion of your screen as many Psychtoolbox functions used to generate the stimuli write on screen windows.

   When its finished there should be a new subdirectory in the revCorrStim/images folder corresponding to the name assigned to the 'outdir' parameter earlier plus a unique time/date string.
   
   The output directory will contain two subdirectories: 'nocc' and 'occ' which contain the non-occluded and occluded versions of the base images respectively.

   The 'nocc' and 'occ' subdirectories will be further subdivided into "L" and "R" subdirectories which contain the left and right angled base images respectively.

   The main output directory will also always contain an image named 'fixation.png' This image contains only the fixation mark and inducers without the objects present.
   
   Finally, in the interest of transparency and reproducability, the main output directory will also contain a .mat file which contains copies of the images and workspace parameters from when the images were made such that all relevant aspects of how the set was made can be referenced later.
     - Data within the .mat file will all be saved within a structure called 'structOut' which contains the following:
         - structOut.itParNamz (see itParNamz in 6 above)
         - structOut.itrPars (see itrPars in 6 above)
         - structOut.itrParsCom : In this case itrParsCom is the set of parameter combos which met the logical constraints/were used in the output image set.
         - structOut.imOutMat : The output image set matching parameters in itrParsCom
         - structOut.itrParsCom_all : the set of all possible parameter combos
         - structOut.imConParNames : set of parameter names which correspond to /allow you to decipher the parameter tag abbreviations in the descriptive output file names
         - structOut.Eqlz (see Eqlz in 6 above)

## Noise Image Generation

The Matlab code for generating noise images is located in the 'revCorrStim/matlab/NoiseGen/imgKernelNoise' subdirectory.

The two "main functions" of interest for generating noise are currently:
- kernelNoise2File_v2.m : For generating "kernel noise" images
- whiteNoise2File_v1.m : For generating white noise images

Very briefly, here is how you use each:

### kernelNoise2File_v2.m:

This function generates a specified number of Ethan's 'kernel noise' images in an output directory.

In short this noise is essentially a weighted combination of sub-images which are each essentially a "patchwork quilt" tiled with various sized "kernel window" samples from random locations of an input base image. These are each also flipped at random orientations (vertical,horizontal, invert, or any combo of these). The set of sub-images are combined linearly with specified weights to form the output noise image.

How to use:
1) Adjust parameters in the '%% Parameters' section
   
   Comment is provided next to each parameter, but we also include a list of the ones you'll likely want to adjust/know below:
   ```
    imgsPerKrnl_set = {...
        [1,1,1,1,1,1,1,0,0],... % adjusts the total number of subtile images from each respective kernel incorporated into the noise
        };

    Repz=100; % sets total number of output images

    jobChunks=10; % sets the number of chunks you want to break the job 
                  % (total reps) into..

    krnlImgWgts=[1,1,1,1,1,1,1,1,1]; % adjusts relative weights applied to subtile images from each respective kernel
    smoothTiles=[0,0,0,0,0,0,0,0,0]; % apply smoothing kernel to intermediate subtile images?
    gSmthK1sz=[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]; % size of gaussian smoothing kernel for subtiles..
    smoothFinal=0; % apply smoothing kernel to final image?
    gSmthK2sz=0.5; % size of gaussian smoothing kernel for final image..

    nReps=Repz/jobChunks; % number of noise image copies/reps you want to run per pass..
    nWorkers=8; % number of cpu nodes used in the parallelized portions..

    %Specify desired size?
    selectSize=1; % if 1, this means that we will use/resize the input image the selected size below
    desiredSize=[512,512]; % (MUST BE SQUARE IMAGE AND DIMZ MUST BE POWER OF 2)

    % Path to base image used to generate kernel noise
    krnlNzBI="/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/images/test6_TexOnly-19-Apr-2024-11-43-14/nocc/L/BaseIm_occ_1_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3315_T_1checkrz_Tr1_2.5_Tr2_0.4_Tal_1_L_0_lum1_127_lum2_127_LWT_1_Lal_1_Ocon_1.png";

    % Base for output directories created..
    baseOutDirName="texOnlyBI_100frms_krnlNz_imzPerKrnl";
   
    % Out Directory Base (parent directory where you want your output dirs
    % created)
    outDirMain="/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/noise";
   ```
 3)   Either hit "Run" or run kernelNoise2File_v2 in the command window.

      This may take a while to finish.

      When its done, there should be a directory full of noise images corresponding to the name string defined by 'baseOutDirName' in the parent directory defined by 'outDirMain'

      Noise frames should be named "noiseSample00001.png" --> "noiseSample#####.png"

      So, make sure you pick a descriptive name for baseOutDirName
      

### whiteNoise2File_v1.m

This function generates a specified number of white noise images in an output directory.

It works very similarly to kernelNoise2File_v2 above ...

How to use:
1) Adjust parameters in the '%% Parameters' section
   
   Comment is provided next to each parameter, but we also include a list of the ones you'll likely want to adjust/know below:
   ```
    Repz=100; % sets total number of output images
    jobChunks=10; % sets the number of chunks you want to break the job 
                  % (total reps) into..
    
    nReps=Repz/jobChunks; % number of noise image copies/reps you want to run per pass..
    nWorkers=8; % number of cpu nodes used in the parallelized portions..
    
    %Specify desired size?
    selectSize=1; % if 1, this means that we will use/resize the input image the selected size below
    desiredSize=[512,512]; % (MUST BE SQUARE IMAGE AND DIMZ MUST BE POWER OF 2)
    
    % Out Directory Base (parent directory where you want your output dirs
    % created)
    outDirMain="/home/eduwell/SynologyDrive/SNAP/projects/revCorrStim/noise";
    baseOutDirName="512by512_whiteNoise_100frms_smpl";
   ```
   
3) Either hit "Run" or run whiteNoise2File_v1 in the command window.
   
   This may take a while to finish.

   When its done, there should be a directory full of noise images corresponding to the name string defined by 'baseOutDirName' in the parent directory defined by 'outDirMain'

   Noise frames should be named "noiseSample00001.png" --> "noiseSample#####.png"

   So, make sure you pick a descriptive name for baseOutDirName


## Machine-Specific Tweaks That Need To Be Done Before Running Stimulus For First Time:

There are still a few kinks to work out to make things work more seamlessly across machines.. 

I know, they are annoying..

Ethan intends to iron these out / come up with a more slick solution to handle these in the future.

However, as of now, there are still a number of "machine specific" things that will need to be done before the first time you run the stimuli:

1) Adjust hard-coded paths in revCorrInstructions_v1.m (this function generates example images/instructions before the task)
   - under the section titled '%% Generate the example images' you'll see section in which a large number of paths are hard coded:
     ```
      %% Generate the example images
      % whichVersion="tex"; % 'lum', 'tex', or 'cr'
      
      if noizType=="krnlNz"
      if whichVersion=="lum"
      
          % for lum
          parzIn=struct; % initialize
          parzIn.occImgR=imread(strcat(main_dir,"/images/test4_LumOnly-05-Mar-2024-10-12-48/occ/R/BaseIm_occ_0_ori_m__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3214_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_0_lum2_255_LWT_1_Lal_1_Ocon_1.png"));
          parzIn.occImgL=imread(strcat(main_dir,"/images/test4_LumOnly-05-Mar-2024-10-12-48/occ/L/BaseIm_occ_0_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3214_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_0_lum2_255_LWT_1_Lal_1_Ocon_1.png"));
          parzIn.noccImgR=imread(strcat(main_dir,"/images/test4_LumOnly-05-Mar-2024-10-12-48/nocc/R/BaseIm_occ_1_ori_m__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3214_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_0_lum2_255_LWT_1_Lal_1_Ocon_1.png"));
          parzIn.noccImgL=imread(strcat(main_dir,"/images/test4_LumOnly-05-Mar-2024-10-12-48/nocc/L/BaseIm_occ_1_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3214_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_0_lum2_255_LWT_1_Lal_1_Ocon_1.png"));
          parzIn.fix_im=imread(strcat(main_dir,"/images/test4_LumOnly-05-Mar-2024-10-12-48/fixation.png"));
          parzIn.noiseDir=strcat(main_dir,"/noise/lumOnlyBI_20000frms_krnlNz_imzPerKrnl_111111100");
          parzIn.nWt=0.5;
      
      end
      
      if whichVersion=="cr"
      
          % for cr
          parzIn=struct; % initialize
          parzIn.occImgR=imread(strcat(main_dir,"/images/test5_CROnly-11-Mar-2024-14-16-09/occ/R/BaseIm_occ_0_ori_m__CR_1_CRpw_4_CRpl_0_CRal_1_CRob_1_CRdm_3315_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_51_lum2_51_LWT_1_Lal_1_Ocon_1.png"));
          parzIn.occImgL=imread(strcat(main_dir,"/images/test5_CROnly-11-Mar-2024-14-16-09/occ/L/BaseIm_occ_0_ori_z__CR_1_CRpw_4_CRpl_0_CRal_1_CRob_1_CRdm_3315_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_51_lum2_51_LWT_1_Lal_1_Ocon_1.png"));
          parzIn.noccImgR=imread(strcat(main_dir,"/images/test5_CROnly-11-Mar-2024-14-16-09/nocc/R/BaseIm_occ_1_ori_m__CR_1_CRpw_4_CRpl_0_CRal_1_CRob_1_CRdm_3315_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_51_lum2_51_LWT_1_Lal_1_Ocon_1.png"));
          parzIn.noccImgL=imread(strcat(main_dir,"/images/test5_CROnly-11-Mar-2024-14-16-09/nocc/L/BaseIm_occ_1_ori_z__CR_1_CRpw_4_CRpl_0_CRal_1_CRob_1_CRdm_3315_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_51_lum2_51_LWT_1_Lal_1_Ocon_1.png"));
          parzIn.fix_im=imread(strcat(main_dir,"/images/test5_CROnly-11-Mar-2024-14-16-09/fixation.png"));
          parzIn.noiseDir=strcat(main_dir,"/noise/crOnlyBI_20000frms_krnlNz_imzPerKrnl_111111100");
          parzIn.nWt=0.5;
      
      end
      
      if whichVersion=="tex"
      
          % for tex
          parzIn=struct; % initialize
          parzIn.occImgR=imread(strcat(main_dir,"/images/test6_TexOnly-19-Apr-2024-11-43-14/occ/R/BaseIm_occ_0_ori_m__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3315_T_1checkrz_Tr1_2.5_Tr2_0.4_Tal_1_L_0_lum1_127_lum2_127_LWT_1_Lal_1_Ocon_1.png"));
          parzIn.occImgL=imread(strcat(main_dir,"/images/test6_TexOnly-19-Apr-2024-11-43-14/occ/L/BaseIm_occ_0_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3315_T_1checkrz_Tr1_2.5_Tr2_0.4_Tal_1_L_0_lum1_127_lum2_127_LWT_1_Lal_1_Ocon_1.png"));
          parzIn.noccImgR=imread(strcat(main_dir,"/images/test6_TexOnly-19-Apr-2024-11-43-14/nocc/R/BaseIm_occ_1_ori_m__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3315_T_1checkrz_Tr1_2.5_Tr2_0.4_Tal_1_L_0_lum1_127_lum2_127_LWT_1_Lal_1_Ocon_1.png"));
          parzIn.noccImgL=imread(strcat(main_dir,"/images/test6_TexOnly-19-Apr-2024-11-43-14/nocc/L/BaseIm_occ_1_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3315_T_1checkrz_Tr1_2.5_Tr2_0.4_Tal_1_L_0_lum1_127_lum2_127_LWT_1_Lal_1_Ocon_1.png"));
          parzIn.fix_im=imread(strcat(main_dir,"/images/test6_TexOnly-19-Apr-2024-11-43-14/fixation.png"));
          parzIn.noiseDir=strcat(main_dir,"/noise/texOnlyBI_v2_20000frms_krnlNz_imzPerKrnl_111111100");   
          parzIn.nWt=0.5;
      
      end
      end
      
      if noizType=="white"
      if whichVersion=="lum"
      
          % for lum
          parzIn=struct; % initialize
          parzIn.occImgR=imread(strcat(main_dir,"/images/test4_LumOnly-05-Mar-2024-10-12-48/occ/R/BaseIm_occ_0_ori_m__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3214_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_0_lum2_255_LWT_1_Lal_1_Ocon_1.png"));
          parzIn.occImgL=imread(strcat(main_dir,"/images/test4_LumOnly-05-Mar-2024-10-12-48/occ/L/BaseIm_occ_0_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3214_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_0_lum2_255_LWT_1_Lal_1_Ocon_1.png"));
          parzIn.noccImgR=imread(strcat(main_dir,"/images/test4_LumOnly-05-Mar-2024-10-12-48/nocc/R/BaseIm_occ_1_ori_m__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3214_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_0_lum2_255_LWT_1_Lal_1_Ocon_1.png"));
          parzIn.noccImgL=imread(strcat(main_dir,"/images/test4_LumOnly-05-Mar-2024-10-12-48/nocc/L/BaseIm_occ_1_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3214_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_0_lum2_255_LWT_1_Lal_1_Ocon_1.png"));
          parzIn.fix_im=imread(strcat(main_dir,"/images/test4_LumOnly-05-Mar-2024-10-12-48/fixation.png"));
          parzIn.noiseDir=strcat(main_dir,"/noise/512by512_whiteNoise_20000frms_smpl1");
          parzIn.nWt=0.5;
      
      end
      
      if whichVersion=="cr"
      
          % for cr
          parzIn=struct; % initialize
          parzIn.occImgR=imread(strcat(main_dir,"/images/test5_CROnly-11-Mar-2024-14-16-09/occ/R/BaseIm_occ_0_ori_m__CR_1_CRpw_4_CRpl_0_CRal_1_CRob_1_CRdm_3315_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_51_lum2_51_LWT_1_Lal_1_Ocon_1.png"));
          parzIn.occImgL=imread(strcat(main_dir,"/images/test5_CROnly-11-Mar-2024-14-16-09/occ/L/BaseIm_occ_0_ori_z__CR_1_CRpw_4_CRpl_0_CRal_1_CRob_1_CRdm_3315_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_51_lum2_51_LWT_1_Lal_1_Ocon_1.png"));
          parzIn.noccImgR=imread(strcat(main_dir,"/images/test5_CROnly-11-Mar-2024-14-16-09/nocc/R/BaseIm_occ_1_ori_m__CR_1_CRpw_4_CRpl_0_CRal_1_CRob_1_CRdm_3315_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_51_lum2_51_LWT_1_Lal_1_Ocon_1.png"));
          parzIn.noccImgL=imread(strcat(main_dir,"/images/test5_CROnly-11-Mar-2024-14-16-09/nocc/L/BaseIm_occ_1_ori_z__CR_1_CRpw_4_CRpl_0_CRal_1_CRob_1_CRdm_3315_T_0checkrz_Tr1_NA_Tr2_NA_Tal_1_L_1_lum1_51_lum2_51_LWT_1_Lal_1_Ocon_1.png"));
          parzIn.fix_im=imread(strcat(main_dir,"/images/test5_CROnly-11-Mar-2024-14-16-09/fixation.png"));
          parzIn.noiseDir=strcat(main_dir,"/noise/512by512_whiteNoise_20000frms_smpl2");
          parzIn.nWt=0.5;
      
      end
      
      if whichVersion=="tex"
          % for tex
          parzIn=struct; % initialize
          parzIn.occImgR=imread(strcat(main_dir,"/images/test6_TexOnly-19-Apr-2024-11-43-14/occ/R/BaseIm_occ_0_ori_m__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3315_T_1checkrz_Tr1_2.5_Tr2_0.4_Tal_1_L_0_lum1_127_lum2_127_LWT_1_Lal_1_Ocon_1.png"));
          parzIn.occImgL=imread(strcat(main_dir,"/images/test6_TexOnly-19-Apr-2024-11-43-14/occ/L/BaseIm_occ_0_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3315_T_1checkrz_Tr1_2.5_Tr2_0.4_Tal_1_L_0_lum1_127_lum2_127_LWT_1_Lal_1_Ocon_1.png"));
          parzIn.noccImgR=imread(strcat(main_dir,"/images/test6_TexOnly-19-Apr-2024-11-43-14/nocc/R/BaseIm_occ_1_ori_m__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3315_T_1checkrz_Tr1_2.5_Tr2_0.4_Tal_1_L_0_lum1_127_lum2_127_LWT_1_Lal_1_Ocon_1.png"));
          parzIn.noccImgL=imread(strcat(main_dir,"/images/test6_TexOnly-19-Apr-2024-11-43-14/nocc/L/BaseIm_occ_1_ori_z__CR_0_CRpw_4_CRpl_0_CRal_1_CRob_NA_CRdm_3315_T_1checkrz_Tr1_2.5_Tr2_0.4_Tal_1_L_0_lum1_127_lum2_127_LWT_1_Lal_1_Ocon_1.png"));
          parzIn.fix_im=imread(strcat(main_dir,"/images/test6_TexOnly-19-Apr-2024-11-43-14/fixation.png"));
          parzIn.noiseDir=strcat(main_dir,"/noise/512by512_whiteNoise_20000frms_smpl3");   
          parzIn.nWt=0.5;
      
      end
      end
     ```
      You will need to adjust the paths to the right and left occluded and unoccluded example images (occImgR,occImgL,noccImgR,noccImgL), the fixation image (fix_im), and the noise directory path (noiseDir), to the paths on your local machine.

      This will need to be done for the texture (tex), common region (cr), and luminance (lum) versions along with the kernel and white noise directories.

2) Similar hard-coded path adjustments need to be made in RevCorr_QSTmain7.m

   - Hard-coded string portions of "noise_dir" and "image_dir" parameters defined in the section below need to be adjusted to match the paths on your machine:
  
   ```
   if noizType=="white"
    if L==1 && T ==0 && C==0
        noise_dir = strcat(main_dir,"/noise/512by512_whiteNoise_20000frms_smpl1");
        nTag=".png";
    end
    if T==1 && L==0 && L==0
        noise_dir = strcat(main_dir,"/noise/512by512_whiteNoise_20000frms_smpl2");
        nTag=".png";
    end
    if C==1 && L==0 && T==0
        noise_dir = strcat(main_dir,"/noise/512by512_whiteNoise_20000frms_smpl3");
        nTag=".png";
    end
    end
    
    if noizType=="krnlNz"
        if L==1 && T ==0 && C==0
            noise_dir = strcat(main_dir,"/noise/lumOnlyBI_20000frms_krnlNz_imzPerKrnl_111111100");
            nTag=".png";
        end
        if T==1 && L==0 && L==0
            noise_dir = strcat(main_dir,"/noise/texOnlyBI_v2_20000frms_krnlNz_imzPerKrnl_111111100");
            nTag=".png";
        end
        if C==1 && L==0 && T==0
            noise_dir = strcat(main_dir,"/noise/crOnlyBI_20000frms_krnlNz_imzPerKrnl_111111100");
            nTag=".png";
        end
    end
    
    
    %                      QUEST OPTIMIZED PARAMETER:
    %##########################################################################
    
    %  ---------  Difference Parameter Between-Objects Scenarios:  ---------
    % -------------------------------------------------------------------------
    
    if L==1 && T ==0 && C==0
        stimVer="lum";
        % LUMINANCE ONLY
        image_dir = strcat(main_dir,"/images/test4_LumOnly-05-Mar-2024-10-12-48");
        qstParStr="lumDiff"; % diff between lum1/lum2
        parTagz = ["lum1","lum2"];
        obj1Par="lum1";
        obj2Par="lum2";
        nWt=0.75;
    end
    
    if T==1 && L==0 && C==0
        stimVer="tex";
        % TEXTURE ONLY
        image_dir = strcat(main_dir,"/images/test6_TexOnly-19-Apr-2024-11-43-14"); % generalized
        qstParStr="texDiff"; % diff between texture resolution1/texture resolution2
        parTagz = ["Tr1","Tr2"];
        obj1Par="Tr1";
        obj2Par="Tr2";
        nWt=0.75;
    end
    % -------------------------------------------------------------------------
    
    %   ---------  Individual Parameter Single-Object Scenarios:  ---------
    if C==1 && L==0 && T==0
        stimVer="cr";
        % COMMON REGION ONLY
        image_dir = strcat(main_dir,"/images/test5_CROnly-11-Mar-2024-14-16-09"); % generalized
        qstParStr="CRpl"; % common region boundary luminance
        parTagz = ["CRpl"];
        obj1Par="CRpl";
        obj2Par="CRpl";
        nWt=0.75;
    end
    % -------------------------------------------------------------------------
   ```

   - under the section that begins with:
   ```
   %%% get some info. on computer running this script
   ```

   You'll need to add a conditional statement for your local machine.. something like:
   ```
    elseif strcmp(comp.machineName,'YOUR_MACHINE_NAME')
          directory = figs_dir_base;
          kboard = GetKeyboardIndices();
          room = 'ROOMNAME'; % subject testing room 2026, machine closest to door..
          Screen('Preference', 'SkipSyncTests', 1); %EJD: added to skip sync tests
          Screen('Preference','VisualDebugLevel', 0);
          Screen('Preference','SuppressAllWarnings', 1);
   ```

## Stimulus Presentation

The Matlab code for running the stimulus is located in the '/revCorrStim/matlab/stimulus' subdirectory.

A printable copy of the instruction script for experimenters to read to subjects is located at: /revCorrStim/matlab/stimulus/instructions/RevCorr_OGT_Instructions_v1.docx

The "main script" for running the stimulus/experiments is: /revCorrStim/matlab/stimulus/RevCorr_main7.m

How to use it:
1) Type RevCorr_main7 in the command window and hit enter.
2) A welcome banner should appear that looks something like this:

   <p align="center">

  </p>
  
  <p align="center">
    <img src="https://github.com/ejduwell/revCorrStim/blob/main/images/githubExamples/welcomeBannerSnap.png" alt="stim1">
    <br>
  </p>
  <p align="center">

  Press Enter/Return to continue...

3) You will then be prompted to enter a subject id. Enter a 3-letter ID and hit enter:

   <p align="center">

  </p>
  
  <p align="center">
    <img src="https://github.com/ejduwell/revCorrStim/blob/main/images/githubExamples/subjInitPromptSnap.png" alt="stim2">
    <br>
  </p>
  <p align="center">

4) It will then ask whether your want to run the occluded base images, non-occluded base images or both. Specify as requested and hit enter:
   (I chose both)

   <p align="center">

  </p>
  
  <p align="center">
    <img src="https://github.com/ejduwell/revCorrStim/blob/main/images/githubExamples/occNoccBothSnap.png" alt="stim3">
    <br>
  </p>
  <p align="center">

5) 
## CI / SI Generation

The Matlab code for computing CIs/SIs is located in the '/revCorrStim/matlab/CIGen' subdirectory.
