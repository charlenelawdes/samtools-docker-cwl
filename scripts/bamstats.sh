#!/bin/zsh

# Default values
file=""
chr=""
threads=10

# Parse command line arguments using getopts
while getopts ":f:c:t:" opt; do
  case $opt in
    f)
      file=$OPTARG
      ;;
    c)
      chr=$OPTARG
      ;;
    t)
      threads=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG"
      echo "Usage: $0 -f <file> [-c <chr>] -t <threads>"
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument."
      echo "Usage: $0 -f <file> [-c <chr>] -t <threads>"
      exit 1
      ;;
  esac
done

# Check if required input file is provided
if [ -z "$file" ]; then
    echo "Usage: $0 -f <file> [-c <chr>] -t <threads>"
    exit 1
fi

echo OS Info:
echo 
nice lsb_release -a
echo 
echo Samtools version:
echo 
nice samtools --version
echo 
echo Run Info:
echo 
echo "Starting bamstats.sh with $threads threads"
echo 
echo "Input file: $file"



# Check if $chr is empty
if [ -n "$chr" ]; then
    # Calculate average depth for the specified chromosome
    nice samtools depth -@ $threads -r $chr $file | awk '{sum += $3} END {print "Average depth: " sum/NR}'

    # Calculate average depth including deletions for the specified chromosome
    nice samtools depth -J -@ $threads -r $chr $file | awk '{sum += $3} END {print "Average depth including deletions: " sum/NR}'

    # Calculate average read length for the specified chromosome
    nice samtools view -@ $threads -h $file | \
    awk -v chrom="${chr}" -v file="${file}" \
      'BEGIN {
        total_length = 0;
        count = 0
      }
      {
        # Skip header lines and reads not mapping to the specified chromosome
        if ($1 ~ /^@/ || $3 != chrom) next;
        total_length += length($10);
        count++
      }
      END {
        # Print the average read length or a message if no reads were found
        if (count > 0)
          print "Average read length for " file, chrom ": " total_length/count;
        else
          print "No reads found for " chrom
      }'
else
    # Calculate average depth for all chromosomes if $chr is not provided
    nice samtools depth -@ $threads $file | awk '{sum += $3} END {print "Average depth: " sum/NR}'

    # Calculate average depth including deletions for all chromosomes if $chr is not provided
    nice samtools depth -J -@ $threads $file | awk '{sum += $3} END {print "Average depth including deletions: " sum/NR}'

    # Calculate average read length for all chromosomes if $chr is not provided
    nice samtools view -@ $threads -h $file | \
    awk -v file="${file}" \
      'BEGIN {
        total_length = 0;
        count = 0
      }
      {
        # Skip header lines
        if ($1 ~ /^@/) next;
        total_length += length($10);
        count++
      }
      END {
        # Print the average read length or a message if no reads were found
        if (count > 0)
          print "Average read length for " file ": " total_length/count;
        else
          print "No reads found"
      }'
fi
