#!/bin/bash

# change the parameters to path on your working directories
raw_bed_dir="/scratch/zh4nh/data/hmm_universe/training_3col" # folder of raw bed files
chr_bed_dir="/scratch/zh4nh/data/hmm_universe/training_by_chr_3cols" # folder where the transitional files are stored (sorted and unsorted bed files by chromosomes)
sorted_file_path="/scratch/zh4nh/data/hmm_universe/combined/combined_chrsort.bed.bed" # output file path

mkdir -p "$chr_bed_dir"

total_files=$(find "$raw_bed_dir" -type f \( -name "*.bed" -o -name "*.bed.gz" \) | wc -l)
processed_files=0
bar_length=40

echo "Gather regions from input files by chromosome:"

# Process each bed or bed.gz file
for file in "$raw_bed_dir"/*.{bed,bed.gz}; do
    if [[ "$file" == *.gz ]]; then
        # Decompress gzipped files on the fly
        zcat "$file" | awk -v out_dir="$chr_bed_dir" '{print >> out_dir "/" $1 ".chr"}'
    else
        # Process uncompressed BED files
        cat "$file" | awk -v out_dir="$chr_bed_dir" '{print >> out_dir "/" $1 ".chr"}'
    fi
    
    # Update progress
    processed_files=$((processed_files + 1))
    completed=$((processed_files * bar_length / total_files))
    bar=$(printf "%-${bar_length}s" "$(printf "%0.s=" $(seq 1 $completed))")
    printf "\rProcessing Files: [%s] %d/%d" "$bar" "$processed_files" "$total_files"
    
done

echo "Sorting each chromosome file:"

# Sorting and combining
chr_files=("$chr_bed_dir"/*.chr)
total_chr_files=${#chr_files[@]}
processed_chr_files=0

for chr_file in "$chr_bed_dir"/*.chr; do
    chr_name=$(basename "$chr_file")  # Get filename with .bed extension

    paste <(awk -F'\t' '{print $1}' "$chr_file") \
      <(awk -F'\t' '{print $2}' "$chr_file" | sort -n) \
      <(awk -F'\t' '{print $3}' "$chr_file" | sort -n) > "$chr_bed_dir/$chr_name.sorted"
    
    processed_chr_files=$((processed_chr_files + 1))
    completed=$((processed_chr_files * bar_length / total_chr_files))
    bar=$(printf "%-${bar_length}s" "$(printf "%0.s=" $(seq 1 $completed))")
    printf "\rSorting Files: [%s] %d/%d" "$bar" "$processed_chr_files" "$total_chr_files"
done


cd "$chr_bed_dir"

# write sorted chromosome bed files in order to the output file
cat $(ls *.sorted | sort -V) > "$sorted_file_path"
