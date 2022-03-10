#!/bin/sh
# source /project/shefflab/rivanna_config/env.sh
# bulker activate databio/peppro
export PATH=~/.local/bin:$PATH

cd "/project/shefflab/processed/bedhmm/ctcf/"
cat *.bed > /home/ys4aj/research/hmm/uniwig/data/bed/non_ctcf_combined.bed

cd "/home/ys4aj/research/hmm/uniwig/data/bed/"
sort -k1,1V non_ctcf_combined.bed > non_ctcf_combined_chrsort.bed

cd "/home/ys4aj/research/hmm/uniwig/"
time ./bin/uniwig -s -m 25 data/bed/non_ctcf_combined_chrsort.bed `refgenie seek hg38/fasta.chrom_sizes`