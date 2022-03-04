#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <bits/stdc++.h>
#include <vector>
#include <deque>
#include <map>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kseq.h"
#include "khash.h"
#include "bigWig.h"
#include "kxsort.h"

#include <stdbool.h>
#include <getopt.h>

KSTREAM_INIT(gzFile, gzread, 0x10000)


//one genomic region from bed file containing chromosome number, start and end of the region
struct chromosome
{
    std::string chrom;
    std::vector<int> starts;
    std::vector<int> ends;
};

//function to print the regions from bed file
void showChromosomes(std::map<std::string, chromosome> chroms)
{
    std::cout << "\nRegions: ";
    for (std::map<std::string, chromosome>::iterator it=chroms.begin(); it!=chroms.end(); it++)
    {
        chromosome chromosome = it->second;
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

static int smoothFixedFormatBW(int chrSize, int stepSize, int smoothSize, std::vector<int> input, std::string chrom, std::string order)
{
    std::vector<float> temp_values; // for later converting into array values to add to bw

    int countIndex = 1;
    int currentCount = 0;
    int cutSite = 0, previousCut = 0, endSite = 0, iterator = 0;

    std::deque<int> closers;
    cutSite = input[iterator];
    iterator++;
    cutSite -= smoothSize;
    endSite = cutSite + 1 + smoothSize * 2;
    if (cutSite < 1)
    {
        cutSite = 1;
    }

    // Skip until the first cut
    while (countIndex < cutSite)
    {
        countIndex += stepSize; //step
    }
    previousCut = cutSite;

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
            if (countIndex % stepSize == 0)
            {
                // std::cout << currentCount << "\n";
                temp_values.push_back((float) currentCount);
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
            // std::cout << currentCount << "\n";
            temp_values.push_back((float) currentCount);
        }
        ++countIndex;
    }

    char *tempchrom = (char *) malloc(chrom.length()+1);
    strcpy(tempchrom,chrom.c_str());
    // std::cout << "till here segfault not occured" << "\n";
    bigWigFile_t *fp = NULL;
    std::string tempfname;
    if (order.compare("start")) {
        tempfname = "/home/ys4aj/research/hmm/uniwig/data/bw/non-ctcf-combined-ends-chr"+chrom+".bw";
    } else {
        tempfname = "/home/ys4aj/research/hmm/uniwig/data/bw/non-ctcf-combined-chr"+chrom+".bw";
    }
    char *fname = (char *) malloc(tempfname.length()+1);
    strcpy(fname,tempfname.c_str());

    char *chroms[] = {tempchrom};
    uint32_t chrLens[] = {chrSize};

    int n = temp_values.size();
    float values[n];
    for (int i=0; i<n; i++) {
        values[i] = temp_values[i];
    }
    // float values[] = &int_values[0];
    if (bwInit(1<<17) != 0) {
        fprintf(stderr, "Error in bwInit\n");
        return false;
    }

    fp = bwOpen(fname, NULL, "w");
    if (!fp) {
        fprintf(stderr, "Error while opening file\n");
        return false;
    }

    if (bwCreateHdr(fp, 10)) goto createHdrError;
    fp->cl = bwCreateChromList(chroms, chrLens, 1);
    if (!fp->cl) goto createChromListError;
    if (bwWriteHdr(fp)) goto writeHdrError;
    if (bwAddIntervalSpanSteps(fp,tempchrom,input[0],1,1,values,n)) goto addIntervalSpanStepsError;

    bwClose(fp);
    bwCleanup();

    return true;

    createHdrError:
        fprintf(stderr, "Received createHdrError for chr%s\n", chrom.c_str());
        goto error;

    createChromListError:
        fprintf(stderr, "Received createChromListError for chr%s\n", chrom.c_str());
        goto error;

    writeHdrError:
        fprintf(stderr, "Received writeHdrError for chr%s\n", chrom.c_str());
        goto error;

    addIntervalSpanStepsError:
        fprintf(stderr, "Received addIntervalSpanStepsError for chr%s\n", chrom.c_str());
        goto error;

    error:
        bwClose(fp);
        bwCleanup();
        return false;
}

// Parent functions that will be called from python. It will select either the
// fixed or variable function according to argument choice.

static bool sitesToExactWig(int chrSize, int stepSize, int smoothSize, bool variableStep, std::vector<int> input,  std::string chrom)
{
    //std::cout << "\nchromosome in sites to exact wig: \n" << chrom;
    bool status = false;
    if (variableStep)
    {
        status = exactVariableFormat(0, stepSize, input, chrom);
    }
    else
    {
        status = exactFixedFormat(chrSize, stepSize, input, chrom);
    }
    return status;
}
//TODO - change to operate on vector of genomic regions
static bool sitesToSmoothWig(int chrSize, int stepSize, int smoothSize, bool variableStep, std::vector<int> input,  std::string chrom)
{
    bool status = false;
    if (variableStep)
    {
        status = smoothVariableFormat(0, stepSize, smoothSize, input, chrom);
    }
    else
    {
        status = smoothFixedFormat(chrSize, stepSize, smoothSize, input, chrom);
    }
    return status;

    /*
        The strategy here is to make a smoothed signal track (bigwig file) given
        the exact base-pair locations of the signals. We want to extend those
        signals out +/- some number. The problem is, this messes up sorting, so
        you can't simply split every value into a range surrounding it, because
        then you have to re-sort. This script uses an alternative algorithm that
        avoids that resorting step, resulting in better performance.

        We conceptualize a nucleotide signal (or 'cut site') as a start and an
        end. We loop through each value and handle it as a start, while pushing
        the corresponding end onto a deque. We will then pull out the oldest
        'end' items from the deque as we process through the values. Each new
        value increments the emitted value, while each 'closing' value
        decrements it.

        We initiate an deque of 'closers', which are positions that will
        decrement the signal output (end points of a smoothed cut).
    */
}

static bool sitesToSmoothBigWig(int chrSize, int stepSize, int smoothSize, bool variableStep, std::vector<int> input, std::string chrom, std::string order)
{
    bool status = false;
    if (variableStep)
    {
        status = smoothVariableFormat(0, stepSize, smoothSize, input, chrom);
        // smoothVariableFormatBW(0, stepSize, smoothSize, input, chrom, order);
    }
    else
    {
        status = smoothFixedFormatBW(chrSize, stepSize, smoothSize, input, chrom, order);
    }
    return status;
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

std::map<std::string, chromosome> read_bed(const char *bedPath)
{
    //vector of vector of regions to store regions in one vector per chromosome
    std::cout << "\nReading chromosomes" << std::endl;
    //std::cout << "\nInput file: " << bedPath << "\n";
    gzFile fp;
    kstream_t *ks;
    kstring_t str = {0, 0, 0};
    fp = bedPath && strcmp(bedPath, "-") ? gzopen(bedPath, "r") : gzdopen(0, "r");
    if (fp == 0)
    {
        fprintf(stderr, "ERROR: failed to open the input file\n");
        exit(1);
    }
    ks = ks_init(fp);
    // chromosome chr;
    char chrom[100] = "";
    std::map<std::string, chromosome> chromosomes;
    // std::vector<chromosome> chromosomes;

    while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0)
    {
        char *ctg, *rest;
        int32_t st, en;
        ctg = parse_bed(str.s, &st, &en, &rest);

        // std:: cout << "\n" << ctg << "\t" << st << "\t" << en;

        if (ctg) {
            if (strcmp(chrom, ctg) != 0) {
                strcpy(chrom, ctg);
                if (chromosomes.find(chrom) == chromosomes.end()) {
                    chromosome chr;
                    chr.chrom = std::string(chrom);
                    // fprintf(stderr, "Creating a new chromosome: %s\n", chr.chrom.c_str());
                    chromosomes.insert(std::pair<std::string, chromosome>(chrom, chr));
                }
            }
            chromosomes[chrom].starts.push_back(st);
            chromosomes[chrom].ends.push_back(en);
        }
    }

    // sort the starts and ends respectively
    for (std::map<std::string, chromosome>::iterator it=chromosomes.begin(); it!=chromosomes.end(); it++) {
        kx::radix_sort(it->second.starts.begin(),it->second.starts.end());
        kx::radix_sort(it->second.ends.begin(),it->second.ends.end());
    }


    std::cout << "Reading finished" << std::endl;
    free(str.s);
    ks_destroy(ks);
    gzclose(fp);
    return chromosomes;
}

int main(int argc, char *argv[])
{
    bool variableFormat = false;
    int stepSize = 1;
    int smoothSize = 1;

  static struct option long_options[] =
    {
      {"variableFormat",    required_argument, 0, 'v'},
      {"stepSize",    required_argument, 0, 't'},
      {"smoothSize",  required_argument, 0, 'm'},
      {"bedPath",  required_argument, 0, 'b'},
      {"chromSizePath", required_argument, 0, 's'},
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
        default:
            fprintf(stderr, "Usage: %s [-vtme] [file...]\n", argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    const char *bedPath = argv[optind];
    const char *chromSizePath = argv[optind+1];
    std::cerr << "Variable format: " << variableFormat << std::endl;
    std::cerr << "Step size: " << stepSize << std::endl;
    std::cerr << "Smooth size: " << smoothSize << std::endl;
    std::cerr << "bedPath: " << bedPath << std::endl;
    std::cerr << "chromSizePath: " << chromSizePath << std::endl;
    if (bedPath == 0) {
        fprintf(stderr, "ERROR: failed to open the input file\n");
        return 1;
    }
    if (chromSizePath==0) {
        fprintf(stderr, "ERROR: failed to open the chrom size file\n");
        return 1;
    }

    std::map<std::string, chromosome> chromosomes;
    chromosomes = read_bed(bedPath);
    // showChromosomes(chromosomes);

    std::map<std::string, int> chromSizes;
    std::ifstream ReadChromSize(chromSizePath);
    std::string eachSize;
    std::string delim = "\t";
    while (getline(ReadChromSize, eachSize)) {
        // std::cout << "each line is: " << eachSize << std::endl;
        std::string chromname = eachSize.substr(0,eachSize.find(delim));
        int size = stoi(eachSize.substr(eachSize.find(delim),-1));
        // std::cout << chromname << " " << size << std::endl;
        chromSizes.insert(std::pair<std::string, int>(chromname, size));
    }

    // Easy to paralallize with OpenMPI
    // std::cout << "Number of chromosomes: " << chromosomes.size() << std::endl;

    // std::map<std::string, int>::iterator iter;
    // for (iter=chromSizes.begin(); iter!=chromSizes.end(); iter++) {
    //     std::cout << iter->first << " " << iter->second << std::endl;
    // }
    int success = 0, failure = 0;
    std::cout << "\nStart processing each chromosome" << std::endl;

    for (std::map<std::string, chromosome>::iterator it=chromosomes.begin(); it!=chromosomes.end(); it++)
    {
        chromosome chromosome = it->second;
        std::string c = chromosome.chrom;
        /* checking if the chr starts and ends are sorted */
        // fprintf(stderr, "%s\n", c.c_str());
        // std::cout << "start sorted\t\t" << std::is_sorted(chromosome.starts.begin(),chromosome.starts.end()) << std::endl;
        // std::cout << "end sorted\t\t" << std::is_sorted(chromosome.ends.begin(),chromosome.ends.end()) << std::endl;
        int chrsize = chromSizes[c];
        if (chrsize == 0) {
            fprintf(stderr, "%s - not matched in the chrom_size file\n", c.c_str());
            failure++;
            continue;
        }
        std::cout << c << " - uniwig with size " << chrsize << std::endl;
        std::string c_num = c.substr(3,-1); // this is used as chrom name in bigWig.h
        int result_st;
        int result_en;
        if (smoothSize==0) 
        {
            result_st = sitesToExactWig(chrsize, stepSize, smoothSize, variableFormat, chromosome.starts, chromosome.chrom);
            result_en = sitesToExactWig(chrsize, stepSize, smoothSize, variableFormat, chromosome.ends, chromosome.chrom);
        }
        else
        {
            result_st = sitesToSmoothBigWig(chrsize, stepSize, smoothSize, variableFormat, chromosome.starts, c_num, "start");
            result_en = sitesToSmoothBigWig(chrsize, stepSize, smoothSize, variableFormat, chromosome.ends, c_num, "end");
        }
        if (result_st && result_en) {
            success++;
        }
        else {
            failure++;
        }
    }

    std::cout << "Finished with " << success << " success and " << failure << " failure" << std::endl;

    return 0;
}