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

#include <stdbool.h>
#include <getopt.h>

KSTREAM_INIT(gzFile, gzread, 0x10000)


//one genomic region from bed file containing chromosome number, start and end of the region
typedef struct chromosome
{
    std::string chrom;
    std::vector<int> starts;
    std::vector<int> ends;
};

//function to print the regions from bed file
void showChromosomes(std::vector<chromosome> chroms)
{
    std::cout << "\nRegions: ";
    for (int chr_nr = 0; chr_nr < chroms.size(); chr_nr++)
    {
        chromosome chromosome = chroms[chr_nr];
        int length = chromosome.starts.size(); //number of regions
        for (int reg = 0; reg < length; reg++)
        {
            std::string c = chromosome.chrom;
            std::cout <<"\n" << c << "\t" << chromosome.starts[reg] << "\t" << chromosome.ends[reg];
        }
    }
    std::cout << "\n";
}

static bool exactFixedFormat(int chrSize, int stepSize, std::vector<int> input, std::string chrom)
{

    int countIndex = 1;
    int previousCut = 0;
    int cutSite = 0;
    int iterator = 0;
    cutSite = input[iterator];
    iterator++;
    std::cout<<"fixedStep chrom=" <<chrom<< " start=" << input[0] << " step=" << stepSize << "\n";

    // Use fixedStep wiggle format
    // Print out 0s until the first cut
    while (countIndex < cutSite)
    {
        // std::cout << 0 << "\n";
        countIndex += stepSize;
    }
    int currentCount = 1;
    previousCut = cutSite;

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


static bool exactVariableFormat(int chrSize, int stepSize, std::vector<int> input, std::string chrom)
{
    std::cout << "variableStep chrom=" << chrom << "\n";

    // All the countIndex stuff is only required for fixedFormat
    int previousCut = 0;
    int cutSite = 0;
    int currentCount = 1;
    int iterator = 0;

    int nReads = 1;

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

    return true;
}

//TODO - change to operate on vector of genomic regions

static bool smoothVariableFormat(int chrSize, int stepSize, int smoothSize, std::vector<int> input, std::string chrom)
{
    std::cout << "variableStep chrom=" << chrom << "\n";

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
    //std::cout << "\nSmooth Variable format\n";

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

static bool smoothFixedFormat(int chrSize, int stepSize, int smoothSize, std::vector<int> input, std::string chrom)
{
    std::cout<<"fixedStep chrom=" <<chrom<< " start=" << input[0] << " step=" << stepSize << "\n";
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

    //std::cout << "\nSmooth Fixed format\n";

    // Print out 0s until the first cut
    while (countIndex < cutSite)
    {
        //std::cout << 0 << "\n";
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

static bool sitesToExactWig(int chrSize, int stepSize, int smoothSize, bool variableStep, std::vector<int> input,  std::string chrom)
{

    //std::cout << "\nchromosome in sites to exact wig: \n" << chrom;
    if (variableStep)
    {
        exactVariableFormat(0, stepSize, input, chrom);
    }
    else
    {
        exactFixedFormat(chrSize, stepSize, input, chrom);
    }
    return true;
}
//TODO - change to operate on vector of genomic regions
static bool sitesToSmoothWig(int chrSize, int stepSize, int smoothSize, bool variableStep, std::vector<int> input,  std::string chrom)
{
    if (variableStep)
    {
        smoothVariableFormat(0, stepSize, smoothSize, input, chrom);
    }
    else
    {
        smoothFixedFormat(chrSize, stepSize, smoothSize, input, chrom);
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

char *parse_bed(char *s, int32_t *st_, int32_t *en_, char **r)
{
    char *p, *q, *ctg = 0;
    int32_t i, st = -1, en = -1;
    if (r)
        *r = 0;
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
            {
                en = atol(p);
                if (r && c != 0)
                    *r = q, *q = c;
            }
            ++i, p = q + 1;
            if (i == 3 || c == '\0')
                break;
        }
    }
    *st_ = st, *en_ = en;
    return i >= 3 ? ctg : 0;
}

std::vector<chromosome> read_bed(const char *bedPath)
{
    //vector of vector of regions to store regions in one vector per chromosome

    //std::cout << "\nInput file: " << bedPath << "\n";
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
    char chrom[100] = "";
    std::vector<chromosome> chromosomes;

    while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0)
    {
        char *ctg, *rest;
        int32_t st, en;
        ctg = parse_bed(str.s, &st, &en, &rest);

        //std:: cout << "\n" << ctg << "\t" << st << "\t" << en;

        if (strcmp(chrom, "") == 0)
        {
            strcpy(chrom, ctg);
            chr.chrom = std::string(chrom);
            chr.starts.push_back(st);
            chr.ends.push_back(en);
            continue;
        }
        if (ctg)
        {
            if(strcmp(chrom, ctg) < 0)
            {
                //std::cout << "\nI'm here, chrom = " << chrom << ", ctg = " << ctg <<  ", compare result: "<< strcmp(chrom, ctg) <<"\n";
                chromosomes.push_back(chr);
                strcpy(chrom, ctg);
                chr.chrom = std::string(chrom);
                std::vector<int> start;
                std::vector<int> end;
                chr.ends = end;
                chr.starts = start;
            }

            chr.starts.push_back(st);
            chr.ends.push_back(en);
        }
    }

    // sort the starts and ends respectively
    std::sort(chr.starts.begin(),chr.starts.end());
    std::sort(chr.ends.begin(),chr.ends.end());

    chromosomes.push_back(chr);
    //std::cout << "\nFinished reading";
    free(str.s);
    ks_destroy(ks);
    gzclose(fp);
    return chromosomes;
}

int main(int argc, char *argv[])
{ //uniwig bedfile stepsize smoothSize(0 for no smoothing) variableformat(pass 1 if variable format wanted, 0 for fixed) 
   
    bool variableFormat = false;
    int stepSize = 1;
    int smoothSize = 1;
    // bool startMode = true;

  static struct option long_options[] =
    {
      {"variableFormat",    required_argument, 0, 'v'},
      {"stepSize",    required_argument, 0, 't'},
      {"smoothSize",  required_argument, 0, 'm'},
    //   {"startMode",  required_argument, 0, 'e'},
      {"bedPath",  required_argument, 0, 'b'},
      {0, 0, 0, 0}
    };

    int option_index = 0;
    int opt; 
    while ((opt = getopt_long(argc, argv, "vt:m:e", long_options, &option_index)) != -1) {
        switch (opt) {
        case 0:
            fprintf (stderr, "positional argument 1?\n");
        case 'v':
            fprintf (stderr, "option -v\n");
            variableFormat = true; break;
        case 't': 
            fprintf (stderr, "option -t with value '%s'\n", optarg);
            stepSize = atoi(optarg); break;
        case 'm':
            fprintf (stderr, "option -m with value '%s'\n", optarg);
            smoothSize = atoi(optarg); 
            break;
        // case 'e':
        //     fprintf (stderr, "option -v\n");
        //     startMode = false; break;
        default:
            fprintf(stderr, "Usage: %s [-vtme] [file...]\n", argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    // std::cerr << "option_index: " << optind << std::endl;
    // char** positionals;
    // positionals = &argv[optind];
    // for (; *positionals; positionals++)
    //     fprintf(stderr, "Positional: %s\n", *positionals);

    const char *bedPath = argv[optind];
    std::cerr << "Variable format: " << variableFormat << std::endl;
    std::cerr << "Step size: " << stepSize << std::endl;
    std::cerr << "Smooth size: " << smoothSize << std::endl;
    // std::cerr << "Start mode: " << startMode << std::endl;
    std::cerr << "bedPath: " << bedPath << std::endl;
    if (bedPath == 0) {
        fprintf(stderr, "ERROR: failed to open the input file\n");
    }
    // return 0;


    // if (argc < 2)
    // {
    //     std::cout << "No file specified";
    //     return 1;
    // }

    // const char *bedPath = argv[1];
    // int stepSize = atoi(argv[2]); //default 1
    // int smoothSize = atoi(argv[3]); //default 25 - if 0 - no smoothing
    // bool variableFormat = argv[4]; //default - fixed step 
    // bool starts = argv[5]; 

    /* previous CLi: 
        "--mode", "smooth",
                            "--chromsize", "0" if self.variable_step else str(chrom_size),
                            "--step-type", "variable" if self.variable_step else "fixed",
                            "--step-size", str(self.step_size) or "1",
                            "--smooth-length", str(self.smooth_length) or "25"] 
                            */





    //std::cout << "\nFile: " << bedPath << " StepSize: " << stepSize << " SmoothSize: " << smoothSize << " VariableFormat: " << variableFormat << "\n";
    std::vector<chromosome> chromosomes;
    chromosomes = read_bed(bedPath);
    // showChromosomes(chromosomes);



    // Easy to paralallize with OpenMPI
    //std::cout << "Number of chromosomes: " << chromosomes.size();

    for (int chrom; chrom<chromosomes.size(); chrom++)
    {
        chromosome chromosome = chromosomes[chrom];
        std::string c = chromosome.chrom;
        int last_end_id = chromosome.ends.size() - 1;
        int chrsize = chromosome.ends[last_end_id];
        if (smoothSize==0) 
        {
            bool result_st = sitesToExactWig(chrsize, stepSize, smoothSize, variableFormat, chromosome.starts, chromosome.chrom);
            bool result_en = sitesToExactWig(chrsize, stepSize, smoothSize, variableFormat, chromosome.ends, chromosome.chrom);
        }
        else
        {
            bool result_st = sitesToSmoothWig(chrsize, stepSize, smoothSize, variableFormat, chromosome.starts, chromosome.chrom);
            bool result_en = sitesToSmoothWig(chrsize, stepSize, smoothSize, variableFormat, chromosome.ends, chromosome.chrom);
        }
        
    }



    // if(smoothSize == 0)
    // {
    //         for (int chrom; chrom < chromosomes.size(); chrom++)
    //         {
    //             chromosome chromosome = chromosomes[chrom];
    //             std::string c = chromosome.chrom;
    //             int last_end_id = chromosome.ends.size() - 1;
    //             int chrsize = chromosome.ends[last_end_id];

    //             if (startMode) {
    //                 bool result_st = sitesToExactWig(chrsize, stepSize, smoothSize, variableFormat, chromosome.starts, chromosome.chrom);
    //             } else {
    //                 bool result_en = sitesToExactWig(chrsize, stepSize, smoothSize, variableFormat, chromosome.ends, chromosome.chrom);
    //             }
    //         }
        
    // }
    // else if(smoothSize > 0)
    // {
    //      for (int chrom; chrom < chromosomes.size(); chrom++)
    //         {
    //             chromosome chromosome = chromosomes[chrom];
    //             std::string c = chromosome.chrom;
    //             int last_end_id = chromosome.ends.size() - 1;
    //             int chrsize = chromosome.ends[last_end_id];
    //             if (startMode) {
    //                 bool result_starts = sitesToSmoothWig(chrsize, stepSize, smoothSize, variableFormat, chromosome.starts, chromosome.chrom);
    //             } else {
    //                 bool result_ends = sitesToSmoothWig(chrsize, stepSize, smoothSize, variableFormat, chromosome.ends, chromosome.chrom);
    //             }
    //         }
    // }

    return 0;
}