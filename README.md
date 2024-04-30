# revCorrStim

<p align="center">

</p>

<p align="center">
  <img src="https://github.com/ejduwell/revCorrStim/blob/main/images/asciiSnapshot.png" alt="TextyBeast Logo">
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

   Note: This may take some time to run. It will also take over a portion of your screen as many Psychtoolbox functions used to generate the stimuli write the stimuli on screen windows.

   When its finished there should be a new subdirectory in the revCorrStim/images folder corresponding to the name assigned to the 'outdir' parameter earlier plus a unique time/date string.
   
   The output directory will contain two subdirectories: 'nocc' and 'occ' which contain the non-occluded and occluded versions of the base images respectively.

   These will be further subdivided into "L" and "R" subdirectories which contain the left and right angled base images respectively.
   
   It will also contain a .mat file which also contains copies of the images and workspace parameters from when the images were made such that all aspects of how the set was made can be referenced later.
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

## Stimulus Presentation

The Matlab code for running the stimulus is located in the '/revCorrStim/matlab/stimulus' subdirectory.

## CI / SI Generation

The Matlab code for computing CIs/SIs is located in the '/revCorrStim/matlab/CIGen' subdirectory.
