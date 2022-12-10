# Uniwig

Given a set of bed files, we want to produce 2 [BigWig](http://genome.ucsc.edu/goldenPath/help/bigWig.html) files: one track of the start coordinates, one track of the end coordinates, and one track for core coordinates.

## Prerequisites
- [libBigWig](https://github.com/dpryan79/libBigWig) - a C library for processing BigWig files.

You will need to install the `libBigWig` library locally by following these steps:
1. git clone the `libBigWig` repository
```
git clone https://github.com/dpryan79/libBigWig
```
We recommand that you clone this repository outside the `uniwig` repository so that they are in the same level. Which means that if you are in the `uniwig` local repository, you can visit the `libBigWig` local repository by `cd ../libBigWig`.

2. install `libBigWig` with provided `Makefile`
```
cd libBigWig
make install prefix=lib
```
This should add two library files `libBigWig.a` and `libBigWig.so` in your `libBigWig` local repository.

## Compiling uniwig
You can use the included `Makefile` to compile, test, or remove `uniwig`. The compiled program will have the path of `./bin/uniwig` (assuming you are already in the `uniwig` local repository). Below are some commands with different `Makefile` usage.
- to complie `uniwig`
```
make uniwig
```
- to test `uniwig`
```
make tests
```
- to remove `uniwig` (remove `./bin/uniwig`)
```
make clean
```
- to recompile a new version of `uniwig` after change
```
make rebuild
```
Make sure to change the `LIB_DIR` if your `libBigWig` and `uniwig` local repositories are not located at the same level. The default path for `LIB_DIR` is a relative path from `uniwig` local repository to `libBigWig` local repository, which has the library files already installed (see previous steps for installing `libBigWig`)


## Usage:
To use the compiled `uniwig` program located at `./bin/uniwig`, use the command with the following format
```
./bin/uniwig (-s) -m (5) $combined_bed_file_path $chrom_size_file_path $bw_output_file_header
```
With the parameters
- `-s` for specifying whether the provided combined bed file is already sorted by the chromosome number. If `-s` flag is not given, `uniwig` will sort the combined bed file by chromosome number.
- `-m` for specifying the smooth size. The positive integer number given after the `-m` flag will be used to select the window size for smoothing the coordinates.
- `$combined_bed_file_path` for specifying the path of combined bed file. A relative path given would start changing directories from `uniwig` local repository (it is recommanded to put the combined bed file in `./data/combined/`)
- `$chrom_size_file_path` for specifying the path of chromosome size reference file. There is a reference file provided in this repository, located at `./test/hg38.chrom.sizes`, but you are welcome to use any other reference files (such as [refgenie](https://refgenie.databio.org/en/latest/))
- `$bw_output_file_header` for specifying the header of the BigWig output file. There will be three `.bw` files produced, ending with `_start.bw`, `_end.bw`, and `_core.bw`. The input for this parameter will be added to the output file paths as prefixes (it is recommanded to put the output BigWig file in `./data/bw/`)

You can also use the provided shell script for executing the entire workflow: from a set of raw data (bed files) to the final BigWig file output. There are two choices provided in this repository
- `create_unsorted.sh` is the shell script for executing `uniwig` without sorting the combined bed file by chromosome number before it's submitted to `uniwig`. In other words, `uniwig` will be responsible for sorting the chromosomes.
- `create_sorted.sh` is the shell script for executing `uniwig` with sorting the combined bed file by chromosome number before it's submitted to `uniwig`. The sorting is done by shell script commmand `sort`.

For both shell script, the default paths for the required files are provided. It is recommanded to organize your input files by the default choice, but make sure to change the paths to corresponding repositories if you are not using the default paths.
