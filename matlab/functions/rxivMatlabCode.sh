 #!/bin/bash
 
 # Get input parameters
 # --------------------------------------------------------------------------
 baseDir=$1      # input variable #1 should be the full path to the program base directory in which your matlab code resides ...
 outDirBase=$2   # input varuable #2 should be the full path to the output parent/base directory where you want to save the archived code ...
 # --------------------------------------------------------------------------
 
 echo " "
 echo "#######################################################################"
 echo "SAVING ARCHIVED COPIES OF ALL MATLAB CODE ('.M' FILES) IN CURRENT STATE"
 echo "#######################################################################"
 echo " "
 echo "input #1 to rxivMatlabCode.sh  (baseDir) is : $baseDir"
 echo "input #2 to rxivMatlabCode.sh (outDirBase) is : $outDirBase"
 echo " "
 # Go to base directory and copy all matlab '.m'  files in directory structure
 # as they are in their present state into the output directory for archive
 # --------------------------------------------------------------------------
 strtDir=$(pwd) # save starting position
 
 cd $outDirBase # go to output base directory
 mkdir matlabCodeRxiv # make rxiv directory
 cd matlabCodeRxiv # enter it
 outDir=$(pwd) # save full path to currentl location as the output dir
 now=$(date) # get current date and time..
 
 # echo some descriptive info into a README.txt
 echo "This directory was created by the bash function 'rxivMatlabCode.sh' (written by E.J. Duwell, PhD) and contains an archived copy of all files ending in '*.m' contained within the $baseDir directory on $now" > README.txt
 
 cd $baseDir # go to the base directory
 cp --parents `find -name \*.m` $outDir/ # copy all .m files into the output directory..
 
 # compress it..
 cd $outDirBase # go to output base directory
 dateStr=$(date +'%m-%d-%Y-%H%M%S')
 tar -czvf matlabCodeRxiv_$dateStr.tar.gz matlabCodeRxiv # compress to .tar.gz
 rm -rf matlabCodeRxiv
 
 echo " "
 echo "FINISHED ARCHIVING MATLAB CODE..."
 echo "ARCHIVE SAVED AT: $outDirBase/matlabCodeRxiv_$dateStr.tar.gz"
 echo " "
 
