
#include <Python.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <bits/stdc++.h>

#include <deque>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// static PyObject *echo(PyObject *self, PyObject *args) { 
//     char input[500];
//     while(fgets(input, 500, stdin)){  //read from STDIN (aka command-line)
//         printf("%s", input);  //print out what user typed in
//         memset(input, 0, strlen(input));  //reset string to all 0's
//     }
//     return Py_True;
// }

// do we need this?
bool StringToBool(std::string s)
{
   if (s == "true" || s == "TRUE" || s == "True" || s == "T" || s == "t") {
       return true;
   }
   else {
       return false;
   }
}


static PyObject *echo(PyObject *self, PyObject *args) {
    // std::ios_base::sync_with_stdio(false);
    // std::cin.tie(NULL);

    std::string header;
    getline(std::cin >> std::ws, header);  // Grab the first line (header)
    std::cout << header << "\n";
    
    int cutSite = 0;
    int nReads = 0;
    while(std::cin >> cutSite) {
        std::cout << cutSite << "\n";
    }
    std::cerr << "nreads: " << nReads << ":echo:";
    return Py_True;
}



static PyObject *exactFixedFormat(int chrSize, int stepSize) {
    int countIndex = 1;
    int previousCut = 0, cutSite = 0;
    std::cin >> cutSite;  // Grab the first cut

    // Use fixedStep wiggle format
    // Print out 0s until the first cut
    while (countIndex < cutSite) {
        std::cout << 0 << "\n";
        countIndex += stepSize;   
    }
    int currentCount = 1;
    previousCut = cutSite;
    // std::cerr << "First read: " << cutSite << "\n";
    // std::cerr << "Chrom size: " << chrSize << "\n";
    
    // Loop through cuts, converting to wiggle format
    while(std::cin >> cutSite) {
        // if it's a duplicate read...
        if (cutSite == previousCut) { // sum up all reads for this spot.
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
        while (countIndex < cutSite) {
            if (countIndex % stepSize == 0) {
                std::cout << 0 << "\n";
            }
            countIndex++;   
        }
        previousCut = cutSite;
    } // end while

    // Finish the last one cut
    std::cout << currentCount << "\n";

    // Finish chromosome by printing 0s until we each the end.
    while(countIndex <= chrSize) {
        if (countIndex != chrSize) {
            if (countIndex % stepSize == 0) {  
                std::cout << 0 << "\n";
            }
        }
        countIndex++;
    }
    return Py_True;
}


static PyObject *exactVariableFormat(int chrSize, int stepSize) {
    // All the countIndex stuff is only required for fixedFormat
    int previousCut = 0;
    int cutSite = 0;
    int currentCount = 1;

    int nReads = 1;
    std::cin >> previousCut;  // Grab the first cut
    // Loop through cuts, converting to wiggle format
    while(std::cin >> cutSite) {
        nReads++;
        // if it's a duplicate read...
        if (cutSite == previousCut) { // sum up all reads for this spot.
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

    return Py_True;
}


static PyObject *smoothVariableFormat(int chrSize, int stepSize, int smoothSize) {
    // All the countIndex stuff is only required for fixedFormat
    int previousCut = 0;
    int cutSite = 0;
    int currentCount = 0;

    std::deque <int> closers;

    std::cin >> cutSite;        // Grab the first cut   
    cutSite -= smoothSize;
    int endSite = cutSite + 1 + smoothSize*2;   

    if (cutSite < 1){
        cutSite = 1;
    }

    previousCut = cutSite;
    int countIndex = cutSite;

    // Loop through cuts, converting to wiggle format   
    while(std::cin >> cutSite) {   
        cutSite -= smoothSize;  
        ++currentCount; 
        closers.push_back(cutSite + 1 + smoothSize*2);  
        if (cutSite < 1){
            cutSite = 1;
        }

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
                    closers.pop_front();       // removes the first element 
                }   
            } 
            if (currentCount > 0) {
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
    int end = std::min(cutSite + 1 + smoothSize*2, chrSize);
    if (chrSize == 0) {
        end = cutSite + 1 + smoothSize*2;
    }
    while (countIndex <= end) {  
        while (endSite == countIndex) { 
            --currentCount; 
            if (closers.empty()) {  
                endSite = 0;  // Must reset endSite to break return loop    
            } else {    
                endSite = closers.front(); // return reference to first element 
                closers.pop_front();       // removes the first element 
            }   
        } 
        if (currentCount > 0) {
            std::cout << countIndex << "\t" << currentCount << "\n";
        }

        ++countIndex;       
    }   

    return Py_True;
}


static PyObject *smoothFixedFormat(int chrSize, int stepSize, int smoothSize) {
    int countIndex = 1;
    int currentCount = 0;
    int cutSite = 0, previousCut = 0, endSite = 0;  

    std::deque <int> closers;

    std::cin >> cutSite;        // Grab the first cut   
    cutSite -= smoothSize;
    endSite = cutSite + 1 + smoothSize*2;
    if (cutSite < 1){
        cutSite = 1;
    }

    // std::cout << "Smooth Fixed format\n";

    // Print out 0s until the first cut 
    while (countIndex < cutSite) {  
        std::cout << 0 << "\n"; 
        countIndex += stepSize;  //step
    }   
    previousCut = cutSite;  

    // Loop through cuts, converting to wiggle format   
    while(std::cin >> cutSite) {
        cutSite -= smoothSize;  
        // std::cout << "Push: " << cutSite << "\n";
        ++currentCount; 
        closers.push_back(cutSite + 1 +smoothSize*2);  
        if (cutSite < 1){
            cutSite = 1;
        }

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
                    closers.pop_front();       // removes the first element 
                }   
            }   
            if (countIndex % stepSize == 0) {   
                std::cout << currentCount <<"\n";  
            }
            ++countIndex;       
        }
        previousCut = cutSite;  
    } // end while  

    // In c we have to add one here for some reason.
    ++currentCount;

    // Finish chromosome by printing 0s until we each the end.  
    while(countIndex <= chrSize) {  
        while (endSite == countIndex) { 
            --currentCount; 
            if (closers.empty()) {  
                endSite = 0;    
            } else {    
                endSite = closers.front();  // return reference to first element    
                closers.pop_front();        // removes the first element    
            }   
        }   
        if (countIndex % stepSize == 0) {   
                std::cout << currentCount << "\n";  
        }   
        ++countIndex;   
    }   


    return Py_True;
}


// Parent functions that will be called from python. It will select either the
// fixed or variable function according to argument choice.
static PyObject *sitesToExactWig(PyObject *self, PyObject *args) {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL); 
    bool variableStep;
    int chrSize, stepSize;


    // get first argument as an int
    // get second argument as a bool
    // get third argument as an int    
    if(!PyArg_ParseTuple(args, "pii", &variableStep, &chrSize, &stepSize))
        return Py_False;

    // We expect the first line given to be a wig header, which we just echo
    std::string header;
    getline(std::cin >> std::ws, header);  // Grab the first line (header)
    std::cout << header << "\n";

    // std::cerr << "Chrom size: " << chrSize << "\n";
    // std::cerr << "Variable step: " << variableStep << "\n";
    // std::cerr << "Step size: " << stepSize << "\n";

    // Loop through cuts, converting to wiggle format
    if (variableStep) {
        exactVariableFormat(chrSize, stepSize);
    } else {
        exactFixedFormat(chrSize, stepSize);
    }
    return Py_True;
}

static PyObject *sitesToSmoothWig(PyObject *self, PyObject *args) {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL); 
    bool variableStep;
    int chrSize, stepSize, smoothSize;

    // get first argument as an int
    // get second argument as a bool
    // get third argument as an int    
    if(!PyArg_ParseTuple(args, "piii", &variableStep, &chrSize, &stepSize, &smoothSize))
        return Py_False;

    // We expect the first line given to be a wig header, which we just echo
    std::string header;
    getline(std::cin >> std::ws, header);  // Grab the first line (header)
    std::cout << header << "\n";

    if (variableStep) {
        smoothVariableFormat(chrSize, stepSize, smoothSize);
    } else {
        smoothFixedFormat(chrSize, stepSize, smoothSize);
    }
    return Py_True;

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


static PyMethodDef smugwigMethods[] = {
    {"sitesToExactWig", sitesToExactWig, METH_VARARGS, "Convert sites to exact wiggle."},
    {"sitesToSmoothWig", sitesToSmoothWig, METH_VARARGS, "Convert sites to smooth wiggle."},
    {"echo", echo, METH_VARARGS, "Echo back input."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef smugwig_module = {
   PyModuleDef_HEAD_INIT,
   "smugwigc",       /* name of module */
   NULL,            /* module documentation, may be NULL */
   -1,              /* size of per-interpreter state of the module,
                or -1 if the module keeps state in global variables. */
   smugwigMethods
};

PyMODINIT_FUNC PyInit_smugwigc(void) {
    return PyModule_Create(&smugwig_module);
}
