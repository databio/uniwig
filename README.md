# Uniwig

Given a set of bed files, we want to produce 3 [BigWig](http://genome.ucsc.edu/goldenPath/help/bigWig.html) output files: one track of the start coordinates, one track of the end coordinates, and one track for core coordinates.

## Prerequisites

Uniwig will rely on [libBigWig](https://github.com/dpryan79/libBigWig), a C library for processing BigWig files. You will need to compile the `libBigWig` library locally first, then you can compile uniwig. The steps are described below.

## Compiling uniwig

### 1. Clone the `libBigWig` and `uniwig` repositories

We recommend you clone this repository so `libBigWig` and `uniwig` are parallel folders at the same level. 

```
git clone https://github.com/dpryan79/libBigWig
git clone https://github.com/databio/uniwig
```

### 2. Compile `libBigWig` with provided `Makefile`

This should add two library files `libBigWig.a` and `libBigWig.so` in your `libBigWig` local repository:

```
cd libBigWig
make install prefix=lib
```

### 3. Compile uniwig

With libBigWig compiled, you can use the `Makefile` to compile, test, or remove `uniwig`. The compiled program will have the path of `./bin/uniwig` (assuming you are already in the `uniwig` local repository). Below are some commands with different `Makefile` usage.

To compile `uniwig`:

```
cd uniwig
make uniwig
```

To test `uniwig`:

```
make tests
```

To remove `uniwig` (remove `./bin/uniwig`):
```
make clean
```

To recompile a new version of `uniwig` after change:
```
make rebuild
```

If your `libBigWig` and `uniwig` local repositories are *not* located at the same level, then change the `LIB_DIR` in the `Makefile`. The default path for `LIB_DIR` is a relative path from `uniwig` local repository to `libBigWig` local repository, which has the library files already compiled.

## Uniwig usage:

To use the compiled `uniwig` program located at `./bin/uniwig`, use the command with the following format

```
./bin/uniwig (-s) -m (5) -w (1) $combined_bed_file_path $chrom_size_file_path $bw_output_file_header
```

With these parameters:

- `-s` for specifying whether the provided combined bed file is already sorted by the chromosome number. If `-s` flag is not given, `uniwig` will sort the combined bed file by chromosome number.
- `-m` for specifying the smooth size. The positive integer number given after the `-m` flag will be used to select the window size for smoothing the coordinates.
- `-w` for specifying the write size. The data will be written into `.bw` files in chunks to reduce memory usage, and the chunk size (number of lines) will be determined by the `-w` flag.
- `$combined_bed_file_path` for specifying the path of combined bed file. A relative path given would start changing directories from `uniwig` local repository (it is recommended to put the combined bed file in `./data/combined/`)
- `$chrom_size_file_path` for specifying the path of chromosome size reference file. There is a reference file provided in this repository, located at `./test/hg38.chrom.sizes`, but you are welcome to use any other reference files (such as [refgenie](https://refgenie.databio.org/en/latest/))
- `$bw_output_file_header` for specifying the header of the BigWig output file. There will be three `.bw` files produced, ending with `_start.bw`, `_end.bw`, and `_core.bw`. The input for this parameter will be added to the output file paths as prefixes (it is recommended to put the output BigWig file in `./data/bw/`)

You can also use the provided shell script for executing the entire workflow: from a set of raw data (bed files) to the final BigWig file output. There are two choices provided in this repository:
- `create_unsorted.sh` is the shell script for executing `uniwig` without sorting the combined bed file by chromosome number before it's submitted to `uniwig`. In other words, `uniwig` will be responsible for sorting the chromosomes.
- `create_sorted.sh` is the shell script for executing `uniwig` with sorting the combined bed file by chromosome number before it's submitted to `uniwig`. The sorting is done by shell script command `sort`.

For both shell script, the default paths for the required files are provided. It is recommended to organize your input files by the default choice, but make sure to change the paths to corresponding repositories if you are not using the default paths.
