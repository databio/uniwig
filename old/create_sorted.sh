#!/bin/sh
# source /project/shefflab/rivanna_config/env.sh
# bulker activate databio/peppro
export PATH=~/.local/bin:$PATH

cd "/project/shefflab/processed/bedhmm/ctcf/"
for i in *.bed; do sort -k1,1V -k2,2n $i > /home/ys4aj/research/hmm/uniwig/data/tempsort/$i.sorted; done;

cd "/home/ys4aj/research/hmm/uniwig/data/tempsort/"
sort -k1,1V -k2,2n --merge *.sorted > /home/ys4aj/research/hmm/uniwig/data/bed/non_ctcf_combined.bed

cd "/home/ys4aj/research/hmm/uniwig/data/bed/"
sort -k1,1V -k3,3n non_ctcf_combined.bed > non_ctcf_combined_endsort.bed

grep -P "chr1\t" non_ctcf_combined.bed > non_ctcf_combined_chr1.bed
grep -P "chr1\t" non_ctcf_combined_endsort.bed > non_ctcf_combined_chr1_endsort.bed

cd "/home/ys4aj/research/hmm/uniwig/"
time ./bin/uniwig -m 25 -e data/bed/non_ctcf_combined_chr1_endsort.bed | wigToBigWig -clip stdin `refgenie seek hg38/fasta.chrom_sizes` data/bw/non-ctcf-combined-ends-chr1.bw
time ./bin/uniwig -m 25 data/bed/non_ctcf_combined_chr1.bed | wigToBigWig -clip stdin `refgenie seek hg38/fasta.chrom_sizes` data/bw/non-ctcf-combined-chr1.bw

cd "/home/ys4aj/research/hmm/uniwig/data/"
grep -P "chr1\t" `refgenie seek hg38/fasta.chrom_sizes` > temp/chr1.chrom_sizes
time bedtools genomecov -d -i bed/non_ctcf_combined_chr1.bed -g temp/chr1.chrom_sizes | cut -f 3 | sed "1i fixedStep chrom=chr1 start=1 step=1" | wigToBigWig -clip stdin temp/chr1.chrom_sizes bw/non-ctcf-combined-coverage-chr1.bw

# cd "/home/ys4aj/research/hmm/uniwig/"
# time ./bin/uniwig -m 25 /home/ys4aj/research/hmm/yc_uniwig/non_ctcf_combined_chr1.bed > /home/ys4aj/research/hmm/yc_uniwig/temp/non-ctcf-combined-chr1.wig
# time ./bin/uniwig -m 25 /home/ys4aj/research/hmm/yc_uniwig/non_ctcf_combined_chr1.bed > /home/ys4aj/research/hmm/yc_uniwig/temp/chr1show.txt