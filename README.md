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
- Generate Base Images ( [jump to section](https://github.com/ejduwell/revCorrStim/tree/main?tab=readme-ov-file#base-image-generation) )
- Generate Noise Images ( [jump to section](https://github.com/ejduwell/revCorrStim/tree/main?tab=readme-ov-file#noise-image-generation) )
- Run the Stimulus ( [jump to section](https://github.com/ejduwell/revCorrStim/tree/main?tab=readme-ov-file#running-an-experiment) )
- Generate CIs/SIs ( [jump to section](https://github.com/ejduwell/revCorrStim/tree/main?tab=readme-ov-file#ci--si-generation) )

In general, each of these aspects of the revCorrStim repository were intentionally created to stand alone as seperate "modules." This approach allows users to create multiple different versions of an experiment with various base image / noise combinations, but run them all with the same centralized stimulus presentation program and analyze the output data with the same analysis software.

<p align="center">
  <div style="text-align: center;">
  ** NOTE: before you can run the stimulus, you will first need to:
    
  1) generate base & noise images
  2) update a number of parameters to point the stimulus presentation program to these images. 
  
  Also, perhaps it goes without saying, but before you can generate CIs/SIs, you must first run an experiment and have data to analyze..
  </div>
</p>

## Base Image Generation

The Matlab code for generating base images is located in the '/revCorrStim/matlab/BaseImgGen' subdirectory.

The most recent copy of the 'main function' for generating base images is GenRevCorrBaseImsRF_v4.m
and
## Noise Image Generation

The Matlab code for generating noise images is located in the 'revCorrStim/matlab/NoiseGen/imgKernelNoise' subdirectory.

## Running an Experiment

The Matlab code for running the stimulus is located in the '/revCorrStim/matlab/stimulus' subdirectory.

## CI / SI Generation

The Matlab code for computing CIs/SIs is located in the '/revCorrStim/matlab/CIGen' subdirectory.
