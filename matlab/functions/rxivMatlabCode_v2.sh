#!/bin/bash

# Function to copy .m files preserving the directory structure
copy_m_files() {
  # Check if exactly two arguments are passed
  if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <source_directory> <output_directory>"
    exit 1
  fi

  local src_dir="$1"
  local out_dir="$2"
  
 echo " "
 echo "#######################################################################"
 echo "SAVING ARCHIVED COPIES OF ALL MATLAB CODE ('.M' FILES) IN CURRENT STATE"
 echo "#######################################################################"
 echo " "
 echo "input #1 to rxivMatlabCode.sh  (baseDir) is : $src_dir"
 echo "input #2 to rxivMatlabCode.sh (outDirBase) is : $out_dir"
 echo " "

  # Ensure the source directory exists
  if [ ! -d "$src_dir" ]; then
    echo "Source directory does not exist: $src_dir"
    exit 1
  fi

  # Create the output directory if it doesn't exist
  mkdir -p "$out_dir/matlabCodeRxiv"
  
  strtDir=$(pwd) # save starting position
  cd $out_dir/matlabCodeRxiv
  out_dir=$(pwd) # save full path to currentl location as the output dir
  now=$(date) # get current date and time..
 
  # echo some descriptive info into a README.txt
  echo "This directory was created by the bash function 'rxivMatlabCode.sh' (written by E.J. Duwell, PhD) and contains an archived copy of all files ending in '*.m' contained within the $baseDir directory on $now" > README.txt

  cd $strtDir #return to start position..

  # Find all '.m' files and copy them to the output directory
  # while preserving the directory structure
  find "$src_dir" -type f -name '*.m' -print0 | while IFS= read -r -d '' file; do
    # Construct the target directory path
    local target_dir="$out_dir/$(dirname "${file#$src_dir}")"

    # Create the target directory if it doesn't exist
    mkdir -p "$target_dir"

    # Copy the file
    cp "$file" "$target_dir/"
  done

  echo "Copy operation completed."
  echo "Compressing to .tar.gz"
  
  # compress it..
  cd $out_dir # go to output base directory
  cd ..
  dateStr=$(date +'%m-%d-%Y-%H%M%S')
  tar -czvf matlabCodeRxiv_$dateStr.tar.gz matlabCodeRxiv # compress to .tar.gz
  rm -rf matlabCodeRxiv
}

# Calling the function with command line arguments
copy_m_files "$1" "$2"

