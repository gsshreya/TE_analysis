#!/bin/bash

# Define paths
dir1="/gpfs/data/user/shweta_lab/data/TE/TE_discovery/joint_call_samples/insertions_55of57" # source
dir2="/gpfs/data/user/shweta_lab/data/TE/TE_discovery/joint_call_samples/insertions" # destination
listofsamples="/gpfs/data/user/ruchitha/MELT/CallingMELT/SampleListFreeze/2025_01_18/ToInclude/combinedrerun_minus2samples.txt" # samples list

# Define the full set of suffixes to match
suffixes=(
    ".aligned.final.sorted.bam"
    ".aligned.final.sorted.bam.bai"
    ".aligned.pulled.sorted.bam"
    ".aligned.pulled.sorted.bam.bai"
    ".hum_breaks.sorted.bam"
    ".hum_breaks.sorted.bam.bai"
    ".tmp.bed"
    "_md.bam"   # Soft link
    "_md.bam.bai" # Soft link
)

# Include directories ending with _mdtmp
dir_suffix="_mdtmp"

# Create destination directory if it does not exist
mkdir -p "$dir2"

# Read sample names into an array
mapfile -t samples < "$listofsamples"

# Loop through all files and directories in dir1
for item in "$dir1"/*; do
    # Extract the basename of the item
    name=$(basename "$item")
    
    # Extract the first 12 characters of the name
    prefix=${name:0:12}
    
    # Check if the prefix exists in the list of samples
    if printf "%s\n" "${samples[@]}" | grep -q "^$prefix$"; then
        # Check for directories ending with _mdtmp
        if [[ -d "$item" && "$name" == *"$dir_suffix" ]]; then
            # Move the directory
            mv "$item" "$dir2"
            echo "Moved directory: $item"
            continue
        fi

        # Check if the name ends with any of the specified suffixes
        for suffix in "${suffixes[@]}"; do
            if [[ "$name" == *"$suffix" ]]; then
                # Move the file or softlink
                mv "$item" "$dir2"
                echo "Moved: $item"
                break
            fi
        done
    fi
done

echo "Selected files and directories moved successfully!"