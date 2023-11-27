#!/bin/zsh
# Purpose: Merge multiple bam files into one bam file

# Default values
chunk_size=1000
threads=12
bam_directory="."

# Parse command line arguments using getopts
while getopts ":c:t:d:" opt; do
  case $opt in
    c)
      chunk_size=$OPTARG
      ;;
    t)
      threads=$OPTARG
      ;;
    d)
      bam_directory=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG"
      echo "Usage: $0 -c <chunk_size> -t <threads> [-d <bam_directory>]"
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument."
      echo "Usage: $0 -c <chunk_size> -t <threads> [-d <bam_directory>]"
      exit 1
      ;;
  esac
done

# Create an array to store bam file names
bam_array=()
bn=$(basename $(realpath $bam_directory))

# Change to the specified directory or use the current directory
cd $bam_directory || exit 1

# Loop through all files with a .bam extension in the current directory
for file in *.bam; do
    # Add each bam file to the array
    bam_array+=($file)
done

# Get the length of the array
arraylength=${#bam_array[@]}

# Print the length of the array
echo "Array length is $arraylength"

# Loop through the array in chunks
for ((i=0; i<${arraylength}; i+=${chunk_size})); do
    echo "Processing chunk starting from index $i"

    # Get a chunk of the array
    array_chunk=("${bam_array[@]:$i:$chunk_size}")

    # Print the length of the chunk
    echo "Chunk length is ${#array_chunk[@]}"

    # Print the names of files in the chunk
    echo "Files in the chunk: ${array_chunk[@]}"

    # Merge bam files in the chunk using samtools
    nice samtools merge -@ ${threads} -o ${bn}_merged_${i}.bam ${array_chunk[@]}

    # Wait for all jobs to finish
    wait
done

# Wait for all jobs to finish
wait

# Merge all previously created chunks into a single bam file
nice samtools merge -@ ${threads} -o ${bn}_merged_all.bam ${bn}_merged_*.bam

# Print a message indicating successful completion
echo "All bam files merged successfully."
