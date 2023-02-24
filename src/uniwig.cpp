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
void showChromosomes_map(std::map<std::string, chromosome> chroms)
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

void showChromosomes_vec(std::vector<chromosome> chroms)
{
    std::cout << "\nRegions: ";
    for (int chr_nr=0; chr_nr<chroms.size(); chr_nr++)
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

static bool smoothFixedStartEndBW(bigWigFile_t *fp, int chrSize, int stepSize, int smoothSize, std::vector<int> input, std::string chrom, std::string order, int writeSize)
{
    std::vector<uint> temp_values; // for later converting into array values to add to bw

    char* tempchrom = new char[chrom.length()+1];
    strcpy(tempchrom,chrom.c_str());

    int countIndex = 1;
    int currentCount = 0;
    int cutSite = 0, previousCut = 0, endSite = 0, iterator = 0;
    int err = 0;
    float val = 0;
    float* valp = NULL;
    int n = 0;

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
                temp_values.push_back(currentCount);
                val = (float) currentCount;
                valp = &val;
                if ((input[iterator-1]+countIndex-previousCut)%writeSize==0) {
                    n = temp_values.size();
                    valp = new float[n];
                   // std::cout << "\n" << countIndex-n+1 << " - " << countIndex << std::endl;
                    for (int i=0; i<n; i++) {
                   //     std::cout << temp_values[i];
                        valp[i] = (float) temp_values[i];
                    }
                    err = bwAddIntervalSpanSteps(fp,tempchrom,countIndex-n+1,1,1,valp,n);
                    temp_values.clear();
                    delete[] valp;
                    if (err) goto addIntervalSpanStepsError;
                }
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
            temp_values.push_back(currentCount);
            val = (float) currentCount;
            valp = &val;
            if ((input[iterator-1]+countIndex-previousCut)%writeSize==0) {
                n = temp_values.size();
                valp = new float[n];
              //  std::cout << "\n" << input[iterator-1]+countIndex-previousCut-n+1 << " - " << input[iterator-1]+countIndex-previousCut << std::endl;
                for (int i=0; i<n; i++) {
               //     std::cout << temp_values[i];
                    valp[i] = (float) temp_values[i];
                }
                err = bwAddIntervalSpanSteps(fp,tempchrom,countIndex-n+1,1,1,valp,n);
                temp_values.clear();
                delete[] valp;
                if (err) goto addIntervalSpanStepsError;
            }
        }
        ++countIndex;
    }

    n = temp_values.size();
    valp = new float[n];
    // std::cout << "\n" << input[iterator-1]+countIndex-previousCut-n+1 << " - " << input[iterator-1]+countIndex-previousCut << std::endl;
    for (int i=0; i<n; i++) {
        // std::cout << temp_values[i];
        valp[i] = (float) temp_values[i];
    }
    err = bwAddIntervalSpanSteps(fp,tempchrom,countIndex-n+1,1,1,valp,n);
    temp_values.clear();
    delete[] valp;
    if (err) goto addIntervalSpanStepsError;
    delete[] tempchrom;

    return true;

    addIntervalSpanStepsError:
        std::cout << "\n\t\tFailed addIntervalSpanStepsError for " << chrom << " - Error code " << err << std::endl;
        delete[] tempchrom;
        return false;
}


static bool fixedCoreBW(bigWigFile_t *fp, int chrSize, int stepSize, std::vector<int> start, std::vector<int> end, std::string chrom, int writeSize)
{
    std::vector<uint> temp_values; // for later converting into array values to add to bw

    char* tempchrom = new char[chrom.length()+1];
    strcpy(tempchrom,chrom.c_str());

    int countIndex = 1;
    int currentCount = 0;
    int cutSite = 0, previousCut = 0, endSite = 0, iterator = 0;
    int err = 0;
    float val = 0;
    float* valp = NULL;
    int n = 0;

    std::deque<int> closers;
    cutSite = start[iterator];
    endSite = end[iterator];
    iterator++;
    // cutSite -= smoothSize;
    // endSite = cutSite + 1 + smoothSize * 2;
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
    while (iterator < start.size())
    {
        // cutSite = input[iterator];
        cutSite = start[iterator];
        // endSite = end[iterator];
        // cutSite -= smoothSize;
        ++currentCount;
        closers.push_back(end[iterator]);
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
                temp_values.push_back(currentCount);
                val = (float) currentCount;
                valp = &val;
                if ((start[iterator-1]+countIndex-previousCut)%writeSize==0) {
                    n = temp_values.size();
                    valp = new float[n];
                    // std::cout << "\n" << input[iterator-1]+countIndex-previousCut-n+1 << " - " << input[iterator-1]+countIndex-previousCut << std::endl;
                    for (int i=0; i<n; i++) {
                        // std::cout << temp_values[i];
                        valp[i] = (float) temp_values[i];
                    }
                    err = bwAddIntervalSpanSteps(fp,tempchrom,countIndex-n+1,1,1,valp,n);
                    temp_values.clear();
                    delete[] valp;
                    if (err) goto addIntervalSpanStepsError;
                }
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
            temp_values.push_back(currentCount);
            val = (float) currentCount;
            valp = &val;
            if ((start[iterator-1]+countIndex-previousCut)%writeSize==0) {
                n = temp_values.size();
                valp = new float[n];
                // std::cout << "\n" << input[iterator-1]+countIndex-previousCut-n+1 << " - " << input[iterator-1]+countIndex-previousCut << std::endl;
                for (int i=0; i<n; i++) {
                    // std::cout << temp_values[i];
                    valp[i] = (float) temp_values[i];
                }
                err = bwAddIntervalSpanSteps(fp,tempchrom,countIndex-n+1,1,1,valp,n);
                temp_values.clear();
                delete[] valp;
                if (err) goto addIntervalSpanStepsError;
            }
        }
        ++countIndex;
    }

    n = temp_values.size();
    valp = new float[n];
    // std::cout << "\n" << input[iterator-1]+countIndex-previousCut-n+1 << " - " << input[iterator-1]+countIndex-previousCut << std::endl;
    for (int i=0; i<n; i++) {
        // std::cout << temp_values[i];
        valp[i] = (float) temp_values[i];
    }
    err = bwAddIntervalSpanSteps(fp,tempchrom,countIndex-n+1,1,1,valp,n);
    temp_values.clear();
    delete[] valp;
    if (err) goto addIntervalSpanStepsError;
    delete[] tempchrom;

    return true;

    addIntervalSpanStepsError:
        std::cout << "\n\t\tFailed addIntervalSpanStepsError for " << chrom << " - Error code " << err << std::endl;
        delete[] tempchrom;
        return false;
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

std::map<std::string, chromosome> read_bed_map(const char *bedPath)
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


    std::cout << "Reading finished\n" << std::endl;
    free(str.s);
    ks_destroy(ks);
    gzclose(fp);
    return chromosomes;
}

std::vector<chromosome> read_bed_vec(const char *bedPath)
{
    //vector of vector of regions to store regions in one vector per chromosome
    //std::cout << "\nInput file: " << bedPath << "\n";
    std::cout << "\nReading chromosomes" << std::endl;
    gzFile fp;
    kstream_t *ks;
    kstring_t str = {0, 0, 0};
    // int32_t k = 0;
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
            if(strcmp(chrom, ctg) != 0)
            {
                // std::cout << "\nI'm here, chrom = " << chrom << ", ctg = " << ctg <<  ", compare result: "<< strcmp(chrom, ctg) <<"\n";
                kx::radix_sort(chr.starts.begin(),chr.starts.end());
                kx::radix_sort(chr.ends.begin(),chr.ends.end());
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
    kx::radix_sort(chr.starts.begin(),chr.starts.end());
    kx::radix_sort(chr.ends.begin(),chr.ends.end());
    chromosomes.push_back(chr);

    std::cout << "Reading finished\n" << std::endl;
    free(str.s);
    ks_destroy(ks);
    gzclose(fp);
    return chromosomes;
}

void print_help(char* argv) {
    fprintf(stderr, "Usage: %s [-h] [-v] [-s] [-tmw] [file...]\n\n", argv);
    fprintf(stderr, "uniwig -- produce wig/bigwig files from bed files\n\n");
    fprintf(stderr, "required arguments:\n");
    fprintf(stderr, "  -t                   step size\n");
    fprintf(stderr, "  -m                   smooth size\n");
    fprintf(stderr, "  -w                   write size\n");
    fprintf(stderr, "\noptional arguments:\n");
    fprintf(stderr, "  -h                   show help commands\n");
    fprintf(stderr, "  -v                   format variables\n");
    fprintf(stderr, "  -s                   bed files alreaday sorted\n");
}

int main(int argc, char *argv[])
{
    bool variableFormat = false;
    bool sorted = false;
    int stepSize = 1;
    int smoothSize = 1;
    int writeSize = 1;

    static struct option long_options[] =
    {
      {"variableFormat",    required_argument, 0, 'v'},
      {"sorted",    required_argument, 0, 's'},
      {"stepSize",    required_argument, 0, 't'},
      {"smoothSize",  required_argument, 0, 'm'},
      {"writeSize",     required_argument, 0, 'w'},
      {"bedPath",  required_argument, 0, 'b'},
      {"chromSizePath", required_argument, 0, 'c'},
      {"fileHeader", required_argument, 0, 'f'},
      {0, 0, 0, 0}
    };

    int option_index = 0;
    int opt; 
    while ((opt = getopt_long(argc, argv, "hvst:m:w:", long_options, &option_index)) != -1) {
        switch (opt) {
        case 0:
            fprintf (stderr, "positional argument 1?\n");
        case 'h':
            print_help(argv[0]); 
            exit(0);
        case 'v':
            fprintf (stderr, "option -v\n");
            variableFormat = true; break;
        case 's':
            fprintf (stderr, "option -s\n");
            sorted = true; break;
        case 't': 
            fprintf (stderr, "option -t with value '%s'\n", optarg);
            stepSize = atoi(optarg); break;
        case 'm':
            fprintf (stderr, "option -m with value '%s'\n", optarg);
            smoothSize = atoi(optarg); 
            break;
        case 'w':
            fprintf (stderr, "option -w with value '%s'\n", optarg);
            writeSize = atoi(optarg);
            break;
        default:
            print_help(argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    const char *bedPath = argv[optind];
    const char *chromSizePath = argv[optind+1];
    const char *fileHeader = argv[optind+2];
    std::cerr << "Variable format: " << variableFormat << std::endl;
    std::cerr << "Sorted bed file: " << sorted << std::endl;
    std::cerr << "Step size: " << stepSize << std::endl;
    std::cerr << "Smooth size: " << smoothSize << std::endl;
    std::cerr << "Write size: " << writeSize << std::endl;
    std::cerr << "bedPath: " << bedPath << std::endl;
    std::cerr << "chromSizePath: " << chromSizePath << std::endl;
    std::cerr << "FileHeader: " << fileHeader << std::endl;
    if (bedPath == 0) {
        fprintf(stderr, "ERROR: failed to open the input file\n");
        return 1;
    }
    if (chromSizePath==0) {
        fprintf(stderr, "ERROR: failed to open the chrom size file\n");
        return 1;
    }

    std::map<std::string, int> chromSizes;
    std::ifstream ReadChromSize(chromSizePath);
    std::string eachSize;
    std::string delim = "\t";
    // std::vector<char*> temp_chroms;
    // std::vector<uint32_t> temp_chrLens;
    // int x = 0;
    while (getline(ReadChromSize, eachSize)) {
        // std::cout << "each line is: " << eachSize << std::endl;
        std::string chromname = eachSize.substr(0,eachSize.find(delim));
        int size = stoi(eachSize.substr(eachSize.find(delim),-1));
        // std::cout << chromname << " " << size << std::endl;
        chromSizes.insert(std::pair<std::string, int>(chromname, size));
    }


    std::string fnames[3] = {
        std::string(fileHeader)+"_start.bw",
        std::string(fileHeader)+"_end.bw",
        std::string(fileHeader)+"_core.bw"
    };


    if (sorted) {
        std::vector<chromosome> chromosomes;
        chromosomes = read_bed_vec(bedPath);
        // showChromosomes_vec(chromosomes);

        int x = chromosomes.size();
        char *chroms[x];
        uint32_t chrLens[x];
        for (int i=0; i<x; i++) {
            chromosome chromosome = chromosomes[i];
            std::string c = chromosome.chrom;
            char* tempc = new char[c.length()+1];
            strcpy(tempc, c.c_str());
            chroms[i] = tempc;
            chrLens[i] = chromSizes[c];
        }

        for (int j=0; j<3; j++) { // for bw file
            char* fname = new char[fnames[j].length()+1];
            strcpy(fname,fnames[j].c_str());

            bigWigFile_t *fp = NULL;
            fp = bwOpen(fname, NULL, "w");
            if (!fp) {
                fprintf(stderr, "Error while opening file\n");
                return 1;
            }

            if (!bwCreateHdr(fp, 10)) {
                fp->cl = bwCreateChromList(chroms, chrLens, x);
                if (fp->cl) {
                    if (!bwWriteHdr(fp)) {
                        fprintf(stderr, "Successfully wrote header to %s\n", fname);
                    }
                }
            }

            // count for success
            int success = 0, failure = 0;
            std::cout << "Processing each chromosome" << std::endl;

            
            // for chrom, write
            for (int chrom=0; chrom<chromosomes.size(); chrom++) {
                chromosome chromosome = chromosomes[chrom];
                std::string c = chromosome.chrom;
                /* checking if the chr starts and ends are sorted */
                // fprintf(stderr, "%s\n", c.c_str());
                // std::cout << "start sorted\t\t" << std::is_sorted(chromosome.starts.begin(),chromosome.starts.end()) << std::endl;
                // std::cout << "end sorted\t\t" << std::is_sorted(chromosome.ends.begin(),chromosome.ends.end()) << std::endl;
                int chrSize = chromSizes[c];
                if (chrSize == 0) {
                    fprintf(stderr, "\t%s\t- not matched in the chrom_size file\n", c.c_str());
                    failure++;
                    continue;
                }
                else {
                    std::cout << "\t" << c << "\t- uniwig with size " << chrSize << "\t- ";
                    std::string c_num = c.substr(3,-1); // this is used as chrom name in bigWig.h

                    bool result = false;
                    // ignoring smoothSize = 0 and variableFormat = true
                    if (smoothSize!=0 && !variableFormat) {
                        switch (j) {
                            case 0: {
                                std::cout << "start";
                                result = smoothFixedStartEndBW(fp, chrSize, stepSize, smoothSize, chromosome.starts, c, "start", writeSize);
                                break;
                            }
                            case 1: {
                                std::cout << "end";
                                result = smoothFixedStartEndBW(fp, chrSize, stepSize, smoothSize, chromosome.ends, c, "end", writeSize);
                                break;
                            }
                            case 2: {
                                std::cout << "core";
                                result = fixedCoreBW(fp, chrSize, stepSize, chromosome.starts, chromosome.ends, c, writeSize);
                                break;
                            }                            
                        }
                    }
                    if (result) {
                        std::cout << "\t complete" << std::endl;
                        success++;
                    }
                    else {
                        std::cout << "\t failed" << std::endl;
                        failure++;
                    }
                }
            }
            std::cout << "Finished with " << success << " success and " << failure << " failure. Cleaning up buffer..." << std::endl;
            
            delete[] fname;
            bwClose(fp);

            std::cout << "Buffer cleaned\n" << std::endl;
        }
        for (int k=0; k<x; k++) {
            delete[] chroms[k];
        }
    }
    else {
        std::map<std::string, chromosome> chromosomes;
        chromosomes = read_bed_map(bedPath);
        // showChromosomes_map(chromosomes);

        int x = chromosomes.size();
        char *chroms[x];
        uint32_t chrLens[x];
        int i = 0;
        for (std::map<std::string, chromosome>::iterator it=chromosomes.begin(); it!=chromosomes.end(); it++) {
            std::string c = it->first;
            char* tempc = new char[c.length()+1];
            strcpy(tempc, c.c_str());
            chroms[i] = tempc;
            chrLens[i] = chromSizes[c];
            i++;
        }
        
        for (int j=0; j<3; j++) { // for bw file
            char* fname = new char[fnames[j].length()+1];
            strcpy(fname,fnames[j].c_str());

            bigWigFile_t *fp = NULL;
            fp = bwOpen(fname, NULL, "w");
            if (!fp) {
                fprintf(stderr, "Error while opening file\n");
                return 1;
            }

            if (!bwCreateHdr(fp, 10)) {
                fp->cl = bwCreateChromList(chroms, chrLens, x);
                if (fp->cl) {
                    if (!bwWriteHdr(fp)) {
                        fprintf(stderr, "Successfully wrote header to %s\n", fname);
                    }
                }
            }

            // count for success
            int success = 0, failure = 0;
            std::cout << "Processing each chromosome" << std::endl;

            // for chrom, write
            for (std::map<std::string, chromosome>::iterator it=chromosomes.begin(); it!=chromosomes.end(); it++) {
                chromosome chromosome = it->second;
                std::string c = chromosome.chrom;
                /* checking if the chr starts and ends are sorted */
                // fprintf(stderr, "%s\n", c.c_str());
                // std::cout << "start sorted\t\t" << std::is_sorted(chromosome.starts.begin(),chromosome.starts.end()) << std::endl;
                // std::cout << "end sorted\t\t" << std::is_sorted(chromosome.ends.begin(),chromosome.ends.end()) << std::endl;
                int chrSize = chromSizes[c];
                if (chrSize == 0) {
                    fprintf(stderr, "\t%s - not matched in the chrom_size file\n", c.c_str());
                    failure++;
                    continue;
                }
                else {
                    std::cout << "\t" << c << "\t- uniwig with size " << chrSize << "\t- ";
                    std::string c_num = c.substr(3,-1); // this is used as chrom name in bigWig.h

                    bool result = false;
                    // ignoring smoothSize = 0 and variableFormat = true
                    if (smoothSize!=0 && !variableFormat) {
                        switch (j) {
                            case 0: {
                                std::cout << "start";
                                result = smoothFixedStartEndBW(fp, chrSize, stepSize, smoothSize, chromosome.starts, c, "start", writeSize);
                                break;
                            }
                            case 1: {
                                std::cout << "end";
                                result = smoothFixedStartEndBW(fp, chrSize, stepSize, smoothSize, chromosome.ends, c, "end", writeSize);
                                break;
                            }
                            case 2: {
                                std::cout << "core";
                                result = fixedCoreBW(fp, chrSize, stepSize, chromosome.starts, chromosome.ends, c, writeSize);
                                break;
                            }                            
                        }
                    }
                    if (result) {
                        std::cout << "\t complete" << std::endl;
                        success++;
                    }
                    else {
                        std::cout << "\t failed" << std::endl;
                        failure++;
                    }
                }
            }
            std::cout << "Finished with " << success << " success and " << failure << " failure. Cleaning up buffer..." << std::endl;

            delete[] fname;
            bwClose(fp);

            std::cout << "Buffer cleaned\n" << std::endl;
        }
        for (int k=0; k<x; k++) {
            delete[] chroms[k];
        }
    }

    return 0;
}
