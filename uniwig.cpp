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

KSTREAM_INIT(gzFile, gzread, 0x10000)

//one genomic region from bed file containing chromosome number, start and end of the region
typedef struct chromosome
{
    char chrom;
    std::vector<int> start;
    std::vector<int> end;
    //int size; //region count
};

void showChromosomes(std::vector<chromosome> chroms)
{
    for (int chr_nr = 0; chr_nr < chroms.size(); chr_nr++)
    {
        chromosome chromosome = chroms[chr_nr];
        int length = chromosome.start.size(); //number of regions
        for (int reg = 0; reg < length; reg++)
        {
            std::cout << "\n"<< chromosome.chrom << "\t" << chromosome.start[reg] << "\t" << chromosome.end[reg];
        }
    }
}

static bool exactFixedFormat(int chrSize, int stepSize)
{
    int countIndex = 1;
    int previousCut = 0;
    int cutSite = 0;
    std::cin >> cutSite; // Grab the first cut

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
    while (std::cin >> cutSite)
    {
        // if it's a duplicate read...
        if (cutSite == previousCut)
        { // sum up all reads for this spot.
            currentCount++;
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

static bool exactVariableFormat(int chrSize, int stepSize)
{
    // All the countIndex stuff is only required for fixedFormat
    int previousCut = 0;
    int cutSite = 0;
    int currentCount = 1;

    int nReads = 1;
    std::cin >> previousCut; // Grab the first cut
    // Loop through cuts, converting to wiggle format
    while (std::cin >> cutSite)
    {
        nReads++;
        // if it's a duplicate read...
        if (cutSite == previousCut)
        { // sum up all reads for this spot.
            currentCount++;
            continue;
        }

        // otherwise, it makes it past this loop;
        // output the sum of counts for the previous spot
        std::cout << previousCut << "\t" << currentCount << "\n";

        // reset for the current spot
        currentCount = 1;
        previousCut = cutSite;
    } // end while

    // Finish it off with the last 'previousCut'
    std::cout << previousCut << "\t" << currentCount << "\n";

    // This is in here for debugging. We could remove this read counting
    // to speed things up further once things are stable.
    // std::cerr << "Reads processed: " << nReads << "\n";

    return true;
}

//TODO - change to operate on vector of genomic regions

static bool smoothVariableFormat(int chrSize, int stepSize, int smoothSize)
{
    // All the countIndex stuff is only required for fixedFormat
    int previousCut = 0;
    int cutSite = 0;
    int currentCount = 0;

    std::deque<int> closers;

    std::cin >> cutSite; // Grab the first cut
    cutSite -= smoothSize;
    int endSite = cutSite + 1 + smoothSize * 2;

    if (cutSite < 1)
    {
        cutSite = 1;
    }

    previousCut = cutSite;
    int countIndex = cutSite;

    // Loop through cuts, converting to wiggle format
    while (std::cin >> cutSite)
    {
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

static bool smoothFixedFormat(int chrSize, int stepSize, int smoothSize, vector<int> input)
{
    int countIndex = 1;
    int currentCount = 0;
    int cutSite = 0, previousCut = 0, endSite = 0;

    std::deque<int> closers;
    std::cin >> cutSite; // Grab the first cut
    cutSite -= smoothSize;
    endSite = cutSite + 1 + smoothSize * 2;
    if (cutSite < 1)
    {
        cutSite = 1;
    }

    // std::cout << "Smooth Fixed format\n";

    // Print out 0s until the first cut
    while (countIndex < cutSite)
    {
        std::cout << 0 << "\n";
        countIndex += stepSize; //step
    }
    previousCut = cutSite;

    // Loop through cuts, converting to wiggle format
    while (std::cin >> cutSite)
    {
        cutSite -= smoothSize;
        // std::cout << "Push: " << cutSite << "\n";
        ++currentCount;
        closers.push_back(cutSite + 1 + smoothSize * 2);
        if (cutSite < 1)
        {
            cutSite = 1;
        }

        // if it's a duplicate read...
        if (cutSite == previousCut)
        {
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

    return true;
}

// Parent functions that will be called from python. It will select either the
// fixed or variable function according to argument choice.
// TODO Can chromosome size potenatialy be a last record in Bed file?
// TODO - change to operate on vector of genomic regions

static bool sitesToExactWig(int chrSize, int stepSize, int smoothSize, bool variableStep)
{
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);

    // get first argument as an int
    // get second argument as a bool
    // get third argument as an int
    /*if(!PyArg_ParseTuple(args, "pii", &variableStep, &chrSize, &stepSize))
        return Py_False; */

    // We expect the first line given to be a wig header, which we just echo
    std::string header;
    getline(std::cin >> std::ws, header); // Grab the first line (header)
    std::cout << header << "\n";

    // std::cerr << "Chrom size: " << chrSize << "\n";
    // std::cerr << "Variable step: " << variableStep << "\n";
    // std::cerr << "Step size: " << stepSize << "\n";

    // Loop through cuts, converting to wiggle format
    if (variableStep)
    {
        exactVariableFormat(chrSize, stepSize);
    }
    else
    {
        exactFixedFormat(chrSize, stepSize);
    }
    return true;
}
//TODO - change to operate on vector of genomic regions
static bool sitesToSmoothWig(int chrSize, int stepSize, int smoothSize, bool variableStep, vector<int> input)
{
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);
    //bool variableStep;
    //int chrSize, stepSize, smoothSize;
    // We expect the first line given to be a wig header, which we just echo
    //std::string header;
    //getline(std::cin >> std::ws, header);  // Grab the first line (header)
    //std::cout << header << "\n";

    if (variableStep)
    {
        smoothVariableFormat(chrSize, stepSize, smoothSize);
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

    std::cout << "Bed file: " << bedPath;
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
                chr.end = start;
                chr.start = end;

                //std:: cout << "\nRegions size should be 0: "<<regions.size() << "\nChromosome: " << chrom;
            }
            chr.chrom = chrom;
            chr.start.push_back(st);
            chr.end.push_back(en);
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
{
    const char *bedPath = argv[1];
    std::vector<chromosome> chromosomes;
    chromosomes = read_bed(bedPath);
    showChromosomes(chromosomes);
    for(int chrom; chrom <= chromosomes.size(); chrom ++)
    {
        sitesToSmoothWig(6000, 1, 5, false, chromosomes[chrom].start);
    }
    //Needs checking if the file exists && better way to specify the file - CLI?
    //FIle can't contain header
}