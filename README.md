# Uniwig

## usage:

```
uniwig bedfile stepsize smoothSize variableformat > file.wig 

smoothsize 0 - no smoothing
variableFormat 1 - for variable format, 0 - for fixed format

```


test command for chromosome print: 

```
./uniwig test.bed 1 5 0
```

## Makefile:

```
make uniwig
```

```
make test
```

```
make clean
```



```
for i in *.bed; do sort -k1,1V -k2,2n $i > $i.sorted; done;
sort -k1,1V -k2,2n --merge *.sorted > all_sorted.bed
cut -f1,2 all_sorted.bed | uniq -c > starts2.txt
```


for i in `findname Ctcf *.bed`; do sort -k1,1V -k2,2n $i > $i.sorted; done;
sort -k1,1V -k2,2n --merge *.sorted > all_sorted.bed
cut -f1,2 all_sorted.bed | uniq -c > starts2.txt



time ./uniwig -m 25 -e non_ctcf_combined_endsort.bed | wigToBigWig -clip stdin `refgenie seek hg38/fasta.chrom_sizes` non-ctcf-combined-ends.bw      


sort -k1,1V -k3,3n non_ctcf_combined.bed > non_ctcf_combined_endsort.bed

	
grep -P "chr1\t" non_ctcf_combined.bed > non_ctcf_combined_chr1.bed
grep -P "chr1\t" non_ctcf_combined_endsort.bed > non_ctcf_combined_chr1_endsort.bed
wc -l non_ctcf_combined_chr1.bed
wc -l non_ctcf_combined_chr1_endsort.bed
tail non_ctcf_combined_chr1_endsort.bed
tail non_ctcf_combined_chr1.bed

time ./uniwig -m 25 -e non_ctcf_combined_chr1_endsort.bed | wigToBigWig -clip stdin `refgenie seek hg38/fasta.chrom_sizes` data/non-ctcf-combined-ends-chr1.bw      
time ./uniwig -m 25 non_ctcf_combined_chr1.bed | wigToBigWig -clip stdin `refgenie seek hg38/fasta.chrom_sizes` data/non-ctcf-combined-chr1.bw      

ls data/*chr1*
scp data/*chr1* rivi:/project/shefflab/processed/bedhmm/
scp data/*coverage*chr1*.bw rivi:/project/shefflab/processed/bedhmm/
scp data/

echo "fixedStep chrom=chr1 start=1 step=1" > data/non-ctcf-combined-coverage-chr1.wig
bedtools genomecov -d -i non_ctcf_combined_chr1_endsort.bed -g `refgenie seek hg38/fasta.chrom_sizes` | grep "chr1" | cut -f3 >> data/non-ctcf-combined-coverage-chr1.wig

wigToBigWig -clip data/non-ctcf-combined-coverage-chr1.wig `refgenie seek hg38/fasta.chrom_sizes` data/non-ctcf-combined-coverage-chr1.bw

bedtools genomecov -d -i non_ctcf_combined_chr1_endsort.bed -g `refgenie seek hg38/fasta.chrom_sizes` | grep "chr1" | cut -f3 | sed 1i"fixedStep chrom=chr1 start=1 step=1" | head

bedtools genomecov -d -i non_ctcf_combined_chr1.bed-g `refgenie seek hg38/fasta.chrom_sizes` | grep "chr1" | cut -f3 | sed 1i"fixedStep chrom=chr1 start=1 step=1" | wigToBigWig -clip stdin `refgenie seek hg38/fasta.chrom_sizes` data/non-ctcf-combined-coverage-chr1.bw


time bedtools genomecov -d -i non_ctcf_combined_chr1.bed -g `refgenie seek hg38/fasta.chrom_sizes` | grep "chr1" | cut -f3 | sed 1i"fixedStep chrom=chr1 start=1 step=1" > data/non-ctcf-combined-coverage-chr1.wig

wigToBigWig -clip data/non-ctcf-combined-coverage-chr1.wig `refgenie seek hg38/fasta.chrom_sizes` data/non-ctcf-combined-coverage-chr1.bw


