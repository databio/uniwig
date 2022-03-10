#!/bin/sh
# source /project/shefflab/rivanna_config/env.sh
# bulker activate databio/peppro
export PATH=~/.local/bin:$PATH

cd "/project/shefflab/processed/bedhmm/ctcf/"
cat *.bed > /home/ys4aj/research/hmm/uniwig/data/bed/non_ctcf_combined.bed

# cd "/home/ys4aj/research/hmm/uniwig/data/bed/"
# grep -P "chr1\t" non_ctcf_combined.bed > non_ctcf_combined_chr1.bed

cd "/home/ys4aj/research/hmm/uniwig/"
time ./bin/uniwig -m 25 data/bed/non_ctcf_combined.bed `refgenie seek hg38/fasta.chrom_sizes`