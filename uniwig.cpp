#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <bits/stdc++.h>
#include <vector>
#include <deque>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kseq.h"
#include "khash.h"
#include "bigWig.h"

KSTREAM_INIT(gzFile, gzread, 0x10000)

/*TODO list: 
* chromosome size - how? chromosome sizes file?  
* handle header present in bed file
* add header to wiggle format - template in smugwig Python code  
* make other smooth/exact functions work 
* output not as stdout - call bigwig etc. 
* call big wig library instead of calling bigwig using PIPE 
*uncompressed wig file as an out */

//one genomic region from bed file containing chromosome number, start and end of the region
typedef struct chromosome
{
    char chrom;
    std::vector<int> starts;
    std::vector<int> ends;
    //int size; //region count
};

void showChromosomes(std::vector<chromosome> chroms)
{
    std::cout << "\nRegions: ";

    for (int chr_nr = 0; chr_nr < chroms.size(); chr_nr++)
    {
        chromosome chromosome = chroms[chr_nr];
        int length = chromosome.starts.size(); //number of regions
        for (int reg = 0; reg < length; reg++)
        {
            std::cout << "\n"
                      << chromosome.chrom << "\t" << chromosome.starts[reg] << "\t" << chromosome.ends[reg];
        }
        std::cout << "\n";
    }
}

static bool exactFixedFormat(int chrSize, int stepSize, std::vector<int> input)
{
    int countIndex = 1;
    int previousCut = 0;
    int cutSite = 0;
    int iterator = 0;
    cutSite = input[iterator];
    iterator++;
    //std::cin >> cutSite; // Grab the first cut
    std::cout << "\nExact Fixed Format\n";

    // Use fixedStep wiggle format
    // Print out 0s until the first cut
    while (countIndex < cutSite)
    {
        std::cout << 0 << "\n";
        countIndex += stepSize;
    }
    int currentCount = 1;
    previousCut = cutSite;
    // std::cerr << "First read: " << cutSite << "\n";
    // std::cerr << "Chrom size: " << chrSize << "\n";

    // Loop through cuts, converting to wiggle format
    while (iterator < input.size())
    {
        cutSite = input[iterator];
        // if it's a duplicate read...
        if (cutSite == previousCut)
        { // sum up all reads for this spot.
            currentCount++;
            iterator++;
            continue;
        }

        // otherwise, it makes it past this loop;
        // output the sum of counts for the previous spot
        std::cout << currentCount << "\n";
        countIndex++;
        // reset for the current spot
        currentCount = 1;
        // and print out all 0s between them
        while (countIndex < cutSite)
        {
            if (countIndex % stepSize == 0)
            {
                std::cout << 0 << "\n";
            }
            countIndex++;
        }
        previousCut = cutSite;
        iterator++;
    } // end while

    // Finish the last one cut
    std::cout << currentCount << "\n";

    // Finish chromosome by printing 0s until we each the end.
    while (countIndex <= chrSize)
    {
        if (countIndex != chrSize)
        {
            if (countIndex % stepSize == 0)
            {
                std::cout << 0 << "\n";
            }
        }
        countIndex++;
    }
    return true;
}

//TODO - change to operate on vector of genomic regions

static bool exactVariableFormat(int chrSize, int stepSize, std::vector<int> input)
{
    // All the countIndex stuff is only required for fixedFormat
    int previousCut = 0;
    int cutSite = 0;
    int currentCount = 1;
    int iterator = 0;

    int nReads = 1;
    std::cout << "\nExact Variable Format\n";
    //std::cin >> previousCut; // Grab the first cut
    // Loop through cuts, converting to wiggle format
    previousCut = input[iterator];
    iterator++;
    while (iterator < input.size())
    {
        cutSite = input[iterator];
        nReads++;
        // if it's a duplicate read...
        if (cutSite == previousCut)
        { // sum up all reads for this spot.
            currentCount++;
            iterator++;
            continue;
        }

        // otherwise, it makes it past this loop;
        // output the sum of counts for the previous spot
        std::cout << previousCut << "\t" << currentCount << "\n";

        // reset for the current spot
        currentCount = 1;
        previousCut = cutSite;
        iterator++;
    } // end while

    // Finish it off with the last 'previousCut'
    std::cout << previousCut << "\t" << currentCount << "\n";

    // This is in here for debugging. We could remove this read counting
    // to speed things up further once things are stable.
    // std::cerr << "Reads processed: " << nReads << "\n";

    return true;
}

//TODO - change to operate on vector of genomic regions

static bool smoothVariableFormat(int chrSize, int stepSize, int smoothSize, std::vector<int> input)
{
    // All the countIndex stuff is only required for fixedFormat
    int previousCut = 0;
    int cutSite = 0;
    int currentCount = 0;
    int iterator = 0;
    std::deque<int> closers;

    cutSite = input[iterator];
    iterator++;
    //std::cin >> cutSite; // Grab the first cut
    cutSite -= smoothSize;
    int endSite = cutSite + 1 + smoothSize * 2;

    if (cutSite < 1)
    {
        cutSite = 1;
    }

    previousCut = cutSite;
    int countIndex = cutSite;
    std::cout << "\nSmooth Variable format\n";

    // Loop through cuts, converting to wiggle format
    //while (std::cin >> cutSite)
    while (iterator < input.size())
    {
        cutSite = input[iterator];
        cutSite -= smoothSize;
        ++currentCount;
        closers.push_back(cutSite + 1 + smoothSize * 2);
        if (cutSite < 1)
        {
            cutSite = 1;
        }

        // if it's a duplicate read...
        if (cutSite == previousCut)
        {
            iterator++;
            continue; // skip to next read
        }

        //int whileloop = 0;
        while (countIndex < cutSite)
        {
            while (endSite == countIndex)
            {
                --currentCount;
                if (closers.empty())
                {
                    endSite = 0; // Must reset endSite to break return loop
                }
                else
                {
                    endSite = closers.front(); // return reference to first element
                    closers.pop_front();       // removes the first element
                }
            }
            if (currentCount > 0)
            {
                std::cout << countIndex << "\t" << currentCount << "\n";
            }
            iterator++;
            ++countIndex;
        }
        previousCut = cutSite;
    } // end while

    // std::cout << countIndex << "\t" << currentCount << "\n";
    ++currentCount;
    //do it one last time
    // std::cout << "Chromsize:" << "\t" << chrSize << "\n";
    int end = std::min(cutSite + 1 + smoothSize * 2, chrSize);
    if (chrSize == 0)
    {
        end = cutSite + 1 + smoothSize * 2;
    }
    while (countIndex <= end)
    {
        while (endSite == countIndex)
        {
            --currentCount;
            if (closers.empty())
            {
                endSite = 0; // Must reset endSite to break return loop
            }
            else
            {
                endSite = closers.front(); // return reference to first element
                closers.pop_front();       // removes the first element
            }
        }
        if (currentCount > 0)
        {
            std::cout << countIndex << "\t" << currentCount << "\n";
        }

        ++countIndex;
    }

    return true;
}

//TODO - change to operate on vector of genomic regions

static bool smoothFixedFormat(int chrSize, int stepSize, int smoothSize, std::vector<int> input)
{
    int countIndex = 1;
    int currentCount = 0;
    int cutSite = 0, previousCut = 0, endSite = 0, iterator = 0;

    std::deque<int> closers;
    cutSite = input[iterator];
    iterator++;
    //std::cin >> cutSite; // Grab the first cut
    cutSite -= smoothSize;
    endSite = cutSite + 1 + smoothSize * 2;
    if (cutSite < 1)
    {
        cutSite = 1;
    }

    std::cout << "\nSmooth Fixed format\n";

    // Print out 0s until the first cut
    while (countIndex < cutSite)
    {
        std::cout << 0 << "\n";
        countIndex += stepSize; //step
    }
    previousCut = cutSite;

    // Loop through cuts, converting to wiggle format
    //while (std::cin >> cutSite)
    while (iterator < input.size())
    {
        cutSite = input[iterator];
        cutSite -= smoothSize;
        //std::cout << "Push: " << cutSite << "\n";
        ++currentCount;
        closers.push_back(cutSite + 1 + smoothSize * 2);
        if (cutSite < 1)
        {
            cutSite = 1;
        }

        // if it's a duplicate read...
        if (cutSite == previousCut)
        {
            iterator++;
            continue; // skip to next read
        }

        //int whileloop = 0;
        while (countIndex < cutSite)
        {
            while (endSite == countIndex)
            {
                --currentCount;
                if (closers.empty())
                {
                    endSite = 0; // Must reset endSite to break return loop
                }
                else
                {
                    endSite = closers.front(); // return reference to first element
                    closers.pop_front();       // removes the first element
                }
            }
            if (countIndex % stepSize == 0)
            {
                std::cout << currentCount << "\n";
            }
            ++countIndex;
        }
        iterator++;
        previousCut = cutSite;
    } // end while

    // In c we have to add one here for some reason.
    ++currentCount;
    // Finish chromosome by printing 0s until we each the end.
    while (countIndex <= chrSize)
    {
        while (endSite == countIndex)
        {
            --currentCount;
            if (closers.empty())
            {
                endSite = 0;
            }
            else
            {
                endSite = closers.front(); // return reference to first element
                closers.pop_front();       // removes the first element
            }
        }
        if (countIndex % stepSize == 0)
        {
            std::cout << currentCount << "\n";
        }
        ++countIndex;
    }
    //std::cout << "\nHere";
    return true;
}

// Parent functions that will be called from python. It will select either the
// fixed or variable function according to argument choice.

static bool sitesToExactWig(int chrSize, int stepSize, int smoothSize, bool variableStep, std::vector<int> input)
{
    if (variableStep)
    {
        exactVariableFormat(0, stepSize, input);
    }
    else
    {
        exactFixedFormat(chrSize, stepSize, input);
    }
    return true;
}
//TODO - change to operate on vector of genomic regions
static bool sitesToSmoothWig(int chrSize, int stepSize, int smoothSize, bool variableStep, std::vector<int> input)
{
    if (variableStep)
    {
        smoothVariableFormat(0, stepSize, smoothSize, input);
    }
    else
    {
        smoothFixedFormat(chrSize, stepSize, smoothSize, input);
    }
    return true;

    // The strategy here is to make a smoothed signal track (bigwig file) given
    // the exact base-pair locations of the signals. We want to extend those
    // signals out +/- some number. The problem is, this messes up sorting, so
    // you can't simply split every value into a range surrounding it, because
    // then you have to re-sort. This script uses an alternative algorithm that
    // avoids that resorting step, resulting in better performance.

    // We conceptualize a nucleotide signal (or 'cut site') as a start and an
    // end. We loop through each value and handle it as a start, while pushing
    // the corresponding end onto a deque. We will then pull out the oldest
    // 'end' items from the deque as we process through the values. Each new
    // value increments the emitted value, while each 'closing' value
    // decrements it.

    // We initiate an deque of 'closers', which are positions that will
    // decrement the signal output (end points of a smoothed cut).
}

char *parse_bed(char *s, int32_t *st_, int32_t *en_)
{
    char *p, *q, *ctg = 0;
    int32_t i, st = -1, en = -1;
    for (i = 0, p = q = s;; ++q)
    {
        if (*q == '\t' || *q == '\0')
        {
            int c = *q;
            *q = 0;
            if (i == 0)
                ctg = p;
            else if (i == 1)
                st = atol(p);
            else if (i == 2)
                en = atol(p);
            ++i, p = q + 1;
            if (c == '\0')
                break;
        }
    }
    *st_ = st, *en_ = en;
    return i >= 3 ? ctg : 0;
}

std::vector<chromosome> read_bed(const char *bedPath)
{
    //vector of vector of regions to store regions in one vector per chromosme

    std::cout << "\nInput file: " << bedPath << "\n";
    gzFile fp;
    kstream_t *ks;
    kstring_t str = {0, 0, 0};
    int32_t k = 0;
    fp = bedPath && strcmp(bedPath, "-") ? gzopen(bedPath, "r") : gzdopen(0, "r");
    if (fp == 0)
    {
        fprintf(stderr, "ERROR: failed to open the input file\n");
        exit(1);
    }
    ks = ks_init(fp);
    chromosome chr;
    char chrom = 0;
    std::vector<chromosome> chromosomes;
    while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0)
    {
        char *ctg;
        int32_t st, en;
        ctg = parse_bed(str.s, &st, &en);
        //std:: cout << "\n" << ctg << "\t" << st << "\t" << en;

        if (chrom == 0)
        {
            chrom = *ctg;
        }
        if (ctg)
        {
            //a vector of genomic regions to store all bed file regions - not the most efficient way - could be changed to AIList if needed as in IGD
            if (*ctg != chrom)
            {
                chromosomes.push_back(chr);
                chrom = *ctg;
                chr.chrom = chrom;
                std::vector<int> start;
                std::vector<int> end;
                chr.ends = end;
                chr.starts = start;

                //std:: cout << "\nRegions size should be 0: "<<regions.size() << "\nChromosome: " << chrom;
            }
            chr.chrom = chrom;
            chr.starts.push_back(st);
            chr.ends.push_back(en);
            //std:: cout << "\nPushing back regions ...";
        }
    }
    chromosomes.push_back(chr);
    //std::cout << "\nFinished reading";
    free(str.s);
    ks_destroy(ks);
    gzclose(fp);
    return chromosomes;
}

int main(int argc, char *argv[])
{ //uniwig bedfile stepsize smoothSize variableformat(pass only if variable formatg wanted)
    std::cout << "\nNumber of arguments: " << argc;
    if (argc < 2)
    {
        std::cout << "No file specified";
        return 1;
    }
    const char *bedPath = argv[1];
    int stepSize = atoi(argv[2]);
    int smoothSize = atoi(argv[3]);
    bool variableFormat = argv[4];
    std::cout << "\nFile: " << bedPath << " StepSize: " << stepSize << " SmoothSize: " << smoothSize << " VariableFormat: " << variableFormat << "\n";
    std::vector<chromosome> chromosomes;
    chromosomes = read_bed(bedPath);
    showChromosomes(chromosomes);
    //Python code for reference
    /*            if self.compress_output:
                # Output will be a pipe to wigToBigWig
                smugwig_shell = False
                smugwig_cmd = ["smugwigi",
                            "--mode", "smooth",
                            "--chromsize", "0" if self.variable_step else str(chrom_size),
                            "--step-type", "variable" if self.variable_step else "fixed",
                            "--step-size", str(self.step_size) or "1",
                            "--smooth-length", str(self.smooth_length) or "25"]
                smugwig_out = subprocess.PIPE           
            else:
                # Output will go direct to the file if not compressing
                smugwig_cmd = " ".join(["smugwigi",
                            "--mode", "exact",
                            "--chromsize", "0" if self.variable_step else str(chrom_size),
                            "--step-type", "variable" if self.variable_step else "fixed",
                            "--step-size", str(self.step_size) or "1",
                            "--smooth-length", str(self.smooth_length) or "25",
                            ">",  chromOutFileBwSm])
                smugwig_out = None
                smugwig_shell = True
    */
    //bigwig
    /*
        if self.variable_step:
            header_line = "variableStep chrom={chrom}\n".format(chrom=chrom)
        else: 
            header_line = "fixedStep chrom={chrom} start={begin} step={step}\n".format(
                chrom=chrom, begin=begin, step=self.step_size) */

    bigWigFile_t *fp = NULL;

    // Easy to paralallize with OpenMPI
    for (int chrom; chrom < chromosomes.size(); chrom++)
    {
        int chrsize = chromosomes[chrom].ends.back();
        std::cout << "\n************** STARTS ***********\n";
        bool result_ef = sitesToExactWig(chrsize, stepSize, 3, false, chromosomes[chrom].starts);
        bool result_vf = sitesToExactWig(chrsize, stepSize, 3, true, chromosomes[chrom].starts);
        bool result = sitesToSmoothWig(chrsize, stepSize, 3, true, chromosomes[chrom].starts);
        bool result_fixed = sitesToSmoothWig(chrsize, stepSize, 3, false, chromosomes[chrom].starts);

        std::cout << "\n************** ENDS ***********\n";
        bool result_efen = sitesToExactWig(chrsize, stepSize, 3, false, chromosomes[chrom].ends);
        bool result_vfen = sitesToExactWig(chrsize, stepSize, 3, true, chromosomes[chrom].ends);
        bool resulten = sitesToSmoothWig(chrsize, stepSize, 3, true, chromosomes[chrom].ends);
        bool result_fixeden = sitesToSmoothWig(chrsize, stepSize, 3, false, chromosomes[chrom].ends);
    }

    return 0;
}