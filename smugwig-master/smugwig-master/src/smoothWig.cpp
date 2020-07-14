/*	
* cutsToWig	
* Name: cutsToWig.cpp	
* Author: Nathan Sheffield, University of Virginia, 2018	
* c++ conversion: Jason Smith University of Virginia, 2019	
* Date: 04.18.2019	
*	
* This is an fast utility (originally in Perl) that converts cut sites	
* (coordinates) into a wiggle-like output.	
* Input #1: a sorted list of integers, corresponding to bases of interest	
* Input #2: an integer (length of chromosome) to fill up with 0s in at the end	
*	
* Output: newline-separated integers counting how many of each site	
* 		  was present in the input.	
*	
* Perl is the right language for this utility; because this is a high	
* IO task (spitting out hundreds of millions of lines), and Perl is highly	
* optimized for this type of IO process. A corresponding python	
* program will take 10-100 fold longer to do the same thing.	
*	
* Run it like this.	
* 1. Pipe your cuts in via stdin:	
* cat cuts.txt | smoothWig CHROMSIZE SMOOTH_LENGTH > out.wig	
* or you can pass your cuts file on the command line:	
* cutsToWig CHROMSIZE SMOOTH_LENGTH cuts.txt > out.wig	
*/	

//#include <Python.h>	
#include <iostream>	
#include <string>	
#include <boost/container/deque.hpp>	
#include <bits/stdc++.h>	

int main(int argc, char *argv[])	
{	
	std::ios_base::sync_with_stdio(false);	
    std::cin.tie(NULL);	

	int countIndex = 1;	
	int currentCount = 0, chrSize = 0, smoothSize = 0, stepSize = 0;	
	int cutSite = 0, previousCut = 0, endSite = 0;	

	std::string header;	

	if (argc >= 3) {	
		chrSize = std::stoi(argv[1], nullptr);	
		smoothSize = std::stoi(argv[2], nullptr);	
		stepSize = std::stoi(argv[3], nullptr);	
	} else {	
		fprintf(stderr, "Missing required input\n");	
		return EXIT_FAILURE;	
	}	

	getline(std::cin >> std::ws, header);  // Get the header	
	std::cout << header << "\n";	

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

	// We initiate an array of 'closers', which are positions that will	
	// decrement the signal output (end points of a smoothed cut).	
	boost::container::deque<int> closers;	

	std::cin >> cutSite;        // Grab the first cut	
	cutSite -= smoothSize;	
	endSite = cutSite + smoothSize*2;	

	// Print out 0s until the first cut	
	while (countIndex < cutSite) {	
		std::cout << 0 << "\n";	
		countIndex += stepSize;	
	}	
	previousCut = cutSite;	

	// Loop through cuts, converting to wiggle format	
	while(std::cin >> cutSite) {	
		cutSite -= smoothSize;	
		++currentCount;	
		closers.push_back(cutSite + smoothSize*2);	

		// if it's a duplicate read...	
		if (cutSite == previousCut) {	
			continue;  // skip to next read	
		}	

		//int whileloop = 0;	
		while (countIndex < cutSite) {	
			while (endSite == countIndex) {	
				--currentCount;	
				if (closers.empty()) {	
					endSite = 0;  // Must reset endSite to break return loop	
				} else {	
					endSite = closers.front(); // return reference to first element	
					closers.pop_front();	   // removes the first element	
				}	
			}	
			if (countIndex % stepSize == 0) {	
				std::cout << currentCount << "\n";	
			}	
			++countIndex;		
		}	
		previousCut = cutSite;	
	} // end while	

	// Finish chromosome by printing 0s until we each the end.	
	while(countIndex <= chrSize) {	
		while (endSite == countIndex) {	
			--currentCount;	
			if (closers.empty()) {	
				endSite = 0;	
			} else {	
				endSite = closers.front();  // return reference to first element	
				closers.pop_front();		// removes the first element	
			}	
		}	
		if (countIndex % stepSize == 0) {	
			if (countIndex != chrSize) {	
				std::cout << currentCount << "\n";	
			}	
		}	
		++countIndex;	
	}	
	return EXIT_SUCCESS;	
}