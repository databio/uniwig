#!/bin/sh
# source /project/shefflab/rivanna_config/env.sh
# bulker activate databio/peppro
export PATH=~/.local/bin:$PATH

# uniwig directory
UNIWIG_DIR="./"
# refgenie chrom_size file path
CHROMSIZE="./test/hg38.chrom.sizes"

# directory for the raw data (bed files)
RAWDATA_DIR="./data/raw/"
# directory for combined data
COMBDATA_DIR="./data/combined/"
# directory for bw files output
BW_DIR="./data/bw/"

# raw data filename
raw="*.bed"
# unsorted combined data filename
unsorted="combined_unsort.bed"
# file header for bw files output
header="test"

cat $RAWDATA_DIR$raw > $COMBDATA_DIR$unsorted

time ./bin/uniwig -m 25 $COMBDATA_DIR$unsorted $CHROMSIZE $BW_DIR$header