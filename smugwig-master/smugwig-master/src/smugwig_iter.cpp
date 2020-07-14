
#include <Python.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <bits/stdc++.h>

bool StringToBool(std::string s)
{
   if (s == "true" || s == "TRUE" || s == "True" || s == "T" || s == "t") {
       return true;
   }
   else {
       return false;
   }
}


static PyObject *cutsToWigStdin(PyObject *self, PyObject *args)
{
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);

    int countIndex = 1, currentCount = 1;
    int previousCut = 0, cutSite = 0, chrSize = 0;
    bool variableStep = true;

    if(!PyArg_ParseTuple(args, "i", &chrSize))
        return Py_False;

    // if (argc >= 2) {
    //     // chrSize = std::stoi(argv[1], 0);
    //     chrSize = atoi(argv[1]);
    //     // variableStep = StringToBool(argv[2]);   // convert string to bool
    // } else {
    //     fprintf(stderr, "Missing required input\n");
    //     return Py_False;
    // }

    std::string header;
    getline(std::cin >> std::ws, header);  // Grab the first line (header)
    std::cout << header << "\n";

    std::cin >> cutSite;  // Grab the first cut


    if (variableStep) {  // Use variableStep wiggle format
        // Increment until the first cut
        while (countIndex < cutSite) {
            countIndex++;   
        }
        previousCut = cutSite;

        // Loop through cuts, converting to wiggle format
        while(std::cin >> cutSite) {
            // if it's a duplicate read...
            if (cutSite == previousCut) { // sum up all reads for this spot.
                currentCount++;
                continue;
            }

            // otherwise, it makes it past this loop;
            // output the sum of counts for the previous spot
            std::cout << previousCut << "\t" << currentCount << "\n";
            countIndex++;
            // reset for the current spot
            currentCount = 1;
            // increment until cutSite
            while (countIndex < cutSite) {
                countIndex++;   
            }
            previousCut = cutSite;
        } // end while

        // Increment until we each the end.
        while(countIndex <= chrSize) {
            countIndex++;
        }
    } else {  // Use fixedStep wiggle format
              // Print out 0s until the first cut
        while (countIndex < cutSite) {
            std::cout << 0 << "\n";
            countIndex++;   
        }
        previousCut = cutSite;

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
                std::cout << 0 << "\n";
                countIndex++;   
            }
            previousCut = cutSite;
        } // end while

        // Finish chromosome by printing 0s until we each the end.
        while(countIndex <= chrSize) {
            if (countIndex != chrSize) {
                std::cout << 0 << "\n";
            }
            countIndex++;
        }
    }
    return Py_True;
}



static PyObject *cEcho(PyObject *self, PyObject *args) {
    int cutSite = 0;
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);
    while(std::cin >> cutSite) {
        std::cout << cutSite << "\n";
    }

    return Py_True;
}

static PyObject *sumStdin(PyObject *self, PyObject *args)
{
    int result = 0, fitem=0;

    while (std::cin >> fitem) {
        result += fitem;
    }    
    std::cout << result << "\n";
    return Py_BuildValue("d", result);
}


static PyObject *sumIter(PyObject *self, PyObject *args)
{
    // Credit: Luther Blissett for early example
    // Example adapted by Nathan Sheffield for modern python
    PyObject* seq;
    PyObject* item;
    double result;

    /* get one argument as an iterator */
    if(!PyArg_ParseTuple(args, "O", &seq))
        return 0;
    seq = PyObject_GetIter(seq);
    if(!seq)
        return 0;


    /* process data sequentially */
    result = 0.0;
    while((item=PyIter_Next(seq))) {
        PyObject *fitem;
        fitem = PyNumber_Float(item);
        if(!fitem) {
            Py_DECREF(seq);
            Py_DECREF(item);
            PyErr_SetString(PyExc_TypeError, "all items must be numbers");
            return 0;
        }
        result += PyFloat_AS_DOUBLE(fitem);
        Py_DECREF(fitem);
        Py_DECREF(item);
    }    

    /* clean up and return result */
    Py_DECREF(seq);
    return Py_BuildValue("d", result);
}

// This is just a demo function that shows how to do some stuff simply. It will
// accept an iterator and an output stream, and it will just "echo", write items
// from the iterator right into the output stream. It expects the iterator to
// yield ints, and it writes them as ints (so your output stream has to know
// what to do with an int).

static PyObject *echo(PyObject *self, PyObject *args) {
    int cutSite = 0;
    PyObject* seq;
    PyObject* item;
    PyObject* outputStream;
      /* get first argument as an iterator */
    // get second argument as an int
    if(!PyArg_ParseTuple(args, "OO", &seq, &outputStream))
        return 0;
    seq = PyObject_GetIter(seq);
    if(!seq)
        return 0;

   while((item = PyIter_Next(seq))) {
        // while(std::cin >> cutSite) {
        cutSite = PyLong_AsLong(item);
        std::string out_string;
        std::stringstream ss;
        ss << cutSite;
        out_string = ss.str();
        std::cout << cutSite << "\n";
        PyObject_CallMethod(outputStream, "write", "i", cutSite);
    }

    return Py_True;
}



static PyObject *fixedFormat(PyObject* seq, PyObject* stream, int chrSize) {
    int countIndex = 1, currentCount = 0;
    int previousCut = 0, cutSite = 0;
    PyObject* item;

    item = PyIter_Next(seq);
    cutSite = PyLong_AsLong(item);  // Grab the first cut

    while (countIndex < cutSite) {
        // std::cout << 0 << "\n";
        PyObject_CallMethod(stream, "write", "i", 0);
        countIndex++;   
    }
    previousCut = cutSite;
    std::cout << "First cut site: " << countIndex << "\n";

    while((item = PyIter_Next(seq))) {
    // while(std::cin >> cutSite) {
        
        cutSite = PyLong_AsLong(item);
        // if it's a duplicate read...
        if (cutSite == previousCut) { // sum up all reads for this spot.
            currentCount++;
            continue;
        }

        // otherwise, it makes it past this loop;
        // output the sum of counts for the previous spot
        if (currentCount > 10) {
            std::cout << countIndex << ": " << currentCount << "\n";

        }
        PyObject_CallMethod(stream, "write", "i", currentCount);
        countIndex++;
        // reset for the current spot
        currentCount = 1;
        // and print out all 0s between them
        while (countIndex < cutSite) {
            // std::cout << 0 << "\n";
            PyObject_CallMethod(stream, "write", "i", 0);
            countIndex++;   
        }
        previousCut = cutSite;
    } // end while

    // Finish chromosome by printing 0s until we each the end.
    while(countIndex <= chrSize) {
        if (countIndex != chrSize) {
            // std::cout << 0 << "\n";
            PyObject_CallMethod(stream, "write", "i", 0);
        }
        countIndex++;
    }
    return Py_True;
}


static PyObject *variableFormat(PyObject* seq, PyObject* stream, int chrSize) {
    int countIndex = 1, currentCount = 0;
    int previousCut = 0, cutSite = 0;
    PyObject* item;

    // char * myfifo = "/tmp/myfifo"; 
    // mkfifo(myfifo, 0666);
    // fd = open(myfifo, O_WRONLY | O_NONBLOCK);
    // int b = 3;
    // write(fd, &b, sizeof(b));
    // close(fd)

    item = PyIter_Next(seq);
    cutSite = PyLong_AsLong(item);  // Grab the first cut

    while (countIndex < cutSite) {
        // std::cout << 0 << "\n";
        PyObject_CallMethod(stream, "write", "i", 0);
        countIndex++;   
    }
    previousCut = cutSite;
    std::cout << "First cut site: " << countIndex << "\n";
    // Loop through cuts, converting to wiggle format
    while((item = PyIter_Next(seq))) {
         cutSite = PyLong_AsLong(item);
        // if it's a duplicate read...
        if (cutSite == previousCut) { // sum up all reads for this spot.
            currentCount++;
            continue;
        }

        // otherwise, it makes it past this loop;
        // output the sum of counts for the previous spot
        // std::cout << previousCut << "\t" << currentCount << "\n";
        char buffer [50];
        sprintf(buffer, "%i\t%i", previousCut, currentCount);
        std::cout << buffer << "\n";
        // PyObject_CallMethod(stream, "write", "s", buffer);
        countIndex++;
        // reset for the current spot
        currentCount = 1;
        // increment until cutSite
        while (countIndex < cutSite) {
            countIndex++;   
        }
        previousCut = cutSite;
    } // end while

    // Increment until we each the end.
    while(countIndex <= chrSize) {
        countIndex++;
    }
    return Py_True;
}


static PyObject *cutsToWig(PyObject *self, PyObject *args) {
    bool variableStep;
    int chrSize = 0;
    PyObject* seq;
    PyObject* stream;
    PyObject* cutsToWigProcessSm;
    PyObject* bedOut;
    double result;

    // get first argument as an iterator
    // get second argument as an int
    // third, fourth and fifth arguments are python objects (streams to write to)
    if(!PyArg_ParseTuple(args, "OpiOOO", &seq, &variableStep, &chrSize, &stream, &cutsToWigProcessSm, &bedOut))
        return 0;
    seq = PyObject_GetIter(seq);
    if(!seq)
        return 0;

    // Loop through cuts, converting to wiggle format
    if (variableStep) {
        variableFormat(seq, stream, chrSize);
    } else {
        fixedFormat(seq, stream, chrSize);
    }

    return Py_True;
}







static PyMethodDef smugwigMethods[] = {
    {"cutsToWig", cutsToWig, METH_VARARGS, "Convert cuts to wiggle."},
    {"sumIter", sumIter, METH_VARARGS, "Convert cuts to wiggle."},
    {"echo", echo, METH_VARARGS, "Convert cuts to wiggle."},
    {"cEcho", cEcho, METH_VARARGS, "Convert cuts to wiggle."},
    {"sumStdin", sumStdin, METH_VARARGS, "Convert cuts to wiggle."},
    {"cutsToWigStdin", cutsToWigStdin, METH_VARARGS, "Convert cuts to wiggle."},
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
