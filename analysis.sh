# This file contains some interactive exploratory
# analysis using uniwig.


# From within a folder with all your .bed files in it:
```
for i in *.bed; do sort -k1,1V -k3,3n $i > $i.sorted; done;
sort -k1,1V -k2,2n --merge *.sorted > all_sorted.bed
cut -f1,2 all_sorted.bed | uniq -c > starts2.txt
```

# In order for this to work, the file has to be sorted
# we do the sort in 2 steps: First sort each file,
# then use the sort --merge to combine them. This is much faster
# than merge them first, and then sorting the combined file.

# to sort by start coordinate, we first sort by chrom then start:

for i in `findname Ctcf *.bed`; do sort -k1,1V -k2,2n $i > $i.sorted; done;
sort -k1,1V -k2,2n --merge *.sorted > all_sorted.bed

# you could then use a command like this to count the frequency of each position:
cut -f1,2 all_sorted.bed | uniq -c > starts.txt


# to sort by end coordinate, we first sort by chrom then end:
for i in `findname Ctcf *.bed`; do sort -k1,1V -k3,3n $i > $i.sorted; done;
sort -k1,1V -k3,3n --merge *.sorted > all_sorted.bed


# Here's how you'd create it for the whole genome? not sure if this works.
#time ./uniwig -m 25 -e non_ctcf_combined_endsort.bed | wigToBigWig -clip stdin `refgenie seek hg38/fasta.chrom_sizes` non-ctcf-combined-ends.bw      


# Once we have it sorted by start, we could just resort that by end
sort -k1,1V -k3,3n non_ctcf_combined.bed > non_ctcf_combined_endsort.bed

	
# Here is how to take the combined file and filter it to just chr1
grep -P "chr1\t" non_ctcf_combined.bed > non_ctcf_combined_chr1.bed
grep -P "chr1\t" non_ctcf_combined_endsort.bed > non_ctcf_combined_chr1_endsort.bed

# confirm that it worked:
wc -l non_ctcf_combined_chr1.bed
wc -l non_ctcf_combined_chr1_endsort.bed
tail non_ctcf_combined_chr1_endsort.bed
tail non_ctcf_combined_chr1.bed

# Given the sorted files, we can pass to uniwig to smooth
time ./uniwig -m 25 -e non_ctcf_combined_chr1_endsort.bed | wigToBigWig -clip stdin `refgenie seek hg38/fasta.chrom_sizes` data/non-ctcf-combined-ends-chr1.bw      
time ./uniwig -m 25 non_ctcf_combined_chr1.bed | wigToBigWig -clip stdin `refgenie seek hg38/fasta.chrom_sizes` data/non-ctcf-combined-chr1.bw      

# Now copy this data over to rivanna
ls data/*chr1*
scp data/*chr1* rivi:/project/shefflab/processed/bedhmm/
scp data/*coverage*chr1*.bw rivi:/project/shefflab/processed/bedhmm/
scp data/

# Here are some other related commands to play around with
echo "fixedStep chrom=chr1 start=1 step=1" > data/non-ctcf-combined-coverage-chr1.wig
bedtools genomecov -d -i non_ctcf_combined_chr1_endsort.bed -g `refgenie seek hg38/fasta.chrom_sizes` | grep "chr1" | cut -f3 >> data/non-ctcf-combined-coverage-chr1.wig

wigToBigWig -clip data/non-ctcf-combined-coverage-chr1.wig `refgenie seek hg38/fasta.chrom_sizes` data/non-ctcf-combined-coverage-chr1.bw

bedtools genomecov -d -i non_ctcf_combined_chr1_endsort.bed -g `refgenie seek hg38/fasta.chrom_sizes` | grep "chr1" | cut -f3 | sed 1i"fixedStep chrom=chr1 start=1 step=1" | head

bedtools genomecov -d -i non_ctcf_combined_chr1.bed-g `refgenie seek hg38/fasta.chrom_sizes` | grep "chr1" | cut -f3 | sed 1i"fixedStep chrom=chr1 start=1 step=1" | wigToBigWig -clip stdin `refgenie seek hg38/fasta.chrom_sizes` data/non-ctcf-combined-coverage-chr1.bw


time bedtools genomecov -d -i non_ctcf_combined_chr1.bed -g `refgenie seek hg38/fasta.chrom_sizes` | grep "chr1" | cut -f3 | sed 1i"fixedStep chrom=chr1 start=1 step=1" > data/non-ctcf-combined-coverage-chr1.wig

wigToBigWig -clip data/non-ctcf-combined-coverage-chr1.wig `refgenie seek hg38/fasta.chrom_sizes` data/non-ctcf-combined-coverage-chr1.bw

