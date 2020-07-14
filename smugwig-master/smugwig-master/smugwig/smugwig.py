#!/usr/bin/env python3

__author__ = ["Nathan C. Sheffield", "Jason Smith"]
__credits__ = []
__license__ = "BSD2"
__version__ = "0.3.1"
__email__ = "nathan@code.databio.org"

import itertools # Used for nested region looping across reads
import logmuse
import numpy
import os
import subprocess
import smugwig
import sys
import tempfile
# import pararead
import pysam

from operator import methodcaller
from argparse import ArgumentParser
from pararead import ParaReadProcessor

import multiprocessing as mp

_LOGGER = None

MODES = ["dnase", "atac"]

# A function object like this will be pickled by the parallel call to map,
# So it cannot contain huge files or the pickling will limit everything.
# For this reason I must rely on global vars for the big stuff.
class CutTracer(ParaReadProcessor):
    """
    A function object that holds parameters and can be called
    on different chromosomes. This extends the ParaReadProcessor object so
    that it can be run in parallel.
    """
    def __init__(self, reads_filename, chrom_sizes_file, temp_parent, nProc,
        limit, verbosity, shift_factor={"+":0, "-":-0}, variable_step=False,
        exactbw=False, summary_filename=None, bedout=False, smoothbw=False,
        smooth_length=25, step_size=5, retain_temp=False, tail_edge=False,
        suppress_compress_output=False):
        """
        param smerge: Should we use the sort merge function that takes two
        separate sorted streams, instead of re-sorting them as one stream?
        """

        # The resultAcronym should be set for each class
        self.resultAcronym="cuttrace"
        self.chrom_sizes_file = chrom_sizes_file
        self.shift_factor = shift_factor

        # Saving a smooth bigwig doubles the processor use for each chrom, so we
        # need to run half as many chroms at a time
        if smoothbw and exactbw:
            _LOGGER.info("Cutting parallel chroms in half to accommodate two tracks.")
            nProc = max(int(nProc / 2), 1)

        super(CutTracer,self).__init__(reads_filename, nProc, outfile=None,
            action=self.resultAcronym, temp_folder_parent_path=temp_parent, limit=limit, 
            allow_unaligned=False, retain_temp=retain_temp)
        self.exactbw = exactbw
        self.summary_filename = summary_filename
        self.verbosity=verbosity
        self.variable_step = variable_step
        self.bedout = bedout
        self.smoothbw = smoothbw
        self.smooth_length = smooth_length
        self.step_size = step_size
        self.tail_edge = tail_edge
        self.compress_output = not suppress_compress_output
        # Confirm that all the commands we will need are callable

        try:
            self.check_command("wigToBigWig")
            self.check_command("perl")
            self.check_command("bigWigCat")
            if self.exactbw:
                self.check_command("smugwigi")
        except AttributeError:
            # backwards compatibility with earlier versions of pararead that did
            # not have a check_command function
            pass

    def register_files(self):
        super(CutTracer, self).register_files()

    def unbuffered_write(self, txt):
        """ 
        Writes unbuffered output by flushing after each stdout.write call
        """
        sys.stdout.write(txt)
        sys.stdout.flush()

    def __call__(self, chrom):
        """
        Workhorse function of the method.
        This function takes a chrom, and processes the reads in that chromosome
        in the input bamfile
        @param chrom: a string with a chromosome, used by bamFile.fetch to
        grab a subset of reads from the bamfile
        """
        chrom_size = self.get_chrom_size(chrom)

        _LOGGER.debug("[Name: " + chrom + "; Size: " + str(chrom_size) + "]")
        reads = self.fetch_chunk(chrom)
        chromOutFile = self._tempf(chrom)
        # os.mkdir(os.path.dirname(chromOutFile))
        chromOutFileBw = os.path.abspath(chromOutFile + ".bw")

        if self.exactbw:

            # Normal sort
            # Universal newlines here allows us to write in text mode.
            sort_cmd = "sort -n"
            sort_cmd = ['sort', '-n']
            sort_process = subprocess.Popen(sort_cmd, shell=False,
                stdout=subprocess.PIPE, stdin=subprocess.PIPE, universal_newlines=True)

            # This is where we'll write our read coordinates
            exact_writer = sort_process.stdin

            if self.compress_output:
                # Output will be a pipe to wigToBigWig
                smugwig_cmd = ["smugwigi",
                            "--mode", "exact",
                            "--chromsize", "0" if self.variable_step else str(chrom_size),
                            "--step-size", str(self.step_size) or "1",
                            "--step-type", "variable" if self.variable_step else "fixed"]
                smugwig_out = subprocess.PIPE
                smugwig_shell = False                
            else:
                # Output will go direct to the file if not compressing
                smugwig_cmd = " ".join(["smugwigi",
                            "--mode", "exact",
                            "--chromsize", "0" if self.variable_step else str(chrom_size),
                            "--step-size", str(self.step_size) or "1",
                            "--step-type", "variable" if self.variable_step else "fixed",
                            ">",  chromOutFileBw])
                smugwig_out = None
                smugwig_shell = True

            _LOGGER.debug(sort_cmd)
            _LOGGER.debug(smugwig_cmd)
            smugwig_process = subprocess.Popen(smugwig_cmd,
                stdin=sort_process.stdout, stdout=smugwig_out, shell=smugwig_shell)

            if self.compress_output:
                wigToBigWig_cmd = ['wigToBigWig', '-clip', '-fixedSummaries',
                     '-keepAllChromosomes', 'stdin',
                     self.chrom_sizes_file, chromOutFileBw]
                #wigToBigWig_cmd = ['head']  # for testing purposes
                wigToBigWig_process = subprocess.Popen(wigToBigWig_cmd, 
                     stdin=smugwig_process.stdout, shell=False)
                _LOGGER.debug(wigToBigWig_cmd)      

        smugwig_process_smooth = None

        if self.smoothbw:

            sort_cmd2 = ['sort', '-n']
            sort_process2 = subprocess.Popen(sort_cmd2, shell=False,
                stdout=subprocess.PIPE, stdin=subprocess.PIPE, universal_newlines=True)

            # This is where we'll write our read coordinates
            smooth_writer = sort_process2.stdin

            cutsToWigSm = os.path.join(os.path.dirname(__file__),
                                       "smoothWig.pl")
            chromOutFileBwSm = chromOutFile + "_smooth.bw"
            tmpFile = chromOutFile + "_cuts.txt"
            cmd = ("sort -n | tee " + tmpFile + " | perl " + cutsToWigSm +
                   " " + str(chrom_size) + " " +  str(self.smooth_length) +
                   " " + str(self.step_size) + " " + str(self.variable_step))


            if self.compress_output:
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


            _LOGGER.debug(smugwig_cmd)
            smugwig_process_smooth = subprocess.Popen(smugwig_cmd, shell=smugwig_shell,
                stdin=sort_process2.stdout, stdout=subprocess.PIPE)

            if self.compress_output:
                wigToBigWig_cmd_smooth = ['wigToBigWig', '-clip', '-fixedSummaries',
                     '-keepAllChromosomes', 'stdin', self.chrom_sizes_file, chromOutFileBwSm]
                # wigToBigWig_cmd_smooth = ['head']
                wigToBigWig_processSm = subprocess.Popen(wigToBigWig_cmd_smooth,
                     stdin=smugwig_process_smooth.stdout)

        bedOut = None
        if self.bedout:
            chromOutFileBed = chromOutFile + ".bed"
            bedOut = open(chromOutFileBed, "w")


        def get_shifted_pos(read, shift_factor):
            """
            Shifts a read according to a shift factor to account for either
            transposition insertion site shifting or DNAse read shifting,
            depending on the strand of the read. Returns the shifted position of
            interest.
            :param read: A pysam read object
            :param shift_factor: A dict with positive or negative integer values
                for keys ["+", "-"], indicating how much (and which direction)
                to shift reads on the + or - strand
            """
            # default
            shifted_pos = None
            if not self.tail_edge:
                if read.flag & 1:  # paired
                    if read.flag == 99:  # paired, mapped in pair, mate reversed, first in pair
                        shifted_pos = read.reference_start + shift_factor["+"]
                        #r.reference_length  # col 8
                    elif read.flag == 147:  # mate of 99
                        shifted_pos = read.reference_end + shift_factor["-"]
                    elif read.flag == 163:  # paired, mapped in pair, mate reversed, second in pair
                        shifted_pos = read.reference_start + shift_factor["+"]
                    elif read.flag == 83:   # mate of 163
                        shifted_pos = read.reference_end + shift_factor["-"]
                else:  # unpaired
                    if read.flag & 16:  # read reverse strand
                        shifted_pos = read.reference_end + shift_factor["-"]
                    else:
                        shifted_pos = read.reference_start + shift_factor["+"]
            else: # Take 3' end of read
                if read.flag & 1:  # paired
                    if read.flag == 99:  # paired, mapped in pair, mate reversed, first in pair
                        shifted_pos = read.reference_end + shift_factor["-"]
                        #r.reference_length  # col 8
                    elif read.flag == 147:  # mate of 99
                        shifted_pos = read.reference_start + shift_factor["+"]
                    elif read.flag == 163:  # paired, mapped in pair, mate reversed, second in pair
                        shifted_pos = read.reference_end + shift_factor["-"]
                    elif read.flag == 83:   # mate of 163
                        shifted_pos = read.reference_start + shift_factor["+"]
                else:  # unpaired
                    if read.flag & 16:  # read reverse strand
                        shifted_pos = read.reference_start + shift_factor["+"]
                    else:
                        shifted_pos = read.reference_end + shift_factor["-"]

            return shifted_pos


        begin = 1

        if self.variable_step:
            header_line = "variableStep chrom={chrom}\n".format(chrom=chrom)
        else: 
            header_line = "fixedStep chrom={chrom} start={begin} step={step}\n".format(
                chrom=chrom, begin=begin, step=self.step_size)

        if self.exactbw:
            # sys.stdout.write(line_to_write)
            exact_writer.write(header_line)
            exact_writer.flush()

        if self.smoothbw:
            # smooth_writer.write(header_line.encode('utf-8'))
            smooth_writer.write(header_line)
            smooth_writer.flush()

        try:
            nreads = 0
            for read in reads:
                nreads+=1
                shifted_pos = get_shifted_pos(read, self.shift_factor)
                line_to_write = str(shifted_pos) + "\n"
                if self.exactbw:
                    # sys.stdout.write(line_to_write)
                    exact_writer.write(line_to_write)
                    exact_writer.flush()

                if self.smoothbw:
                    # smooth_writer.write((str(shifted_pos) + "\n").encode('utf-8'))
                    smooth_writer.write(line_to_write)
                    smooth_writer.flush()

                if self.bedout:
                    strand = "-" if read.is_reverse else "+"
                    # The bed file needs 6 columns (even though some are dummy) 
                    # because MACS says so.
                    bedOut.write("\t".join([
                        chrom,
                        str(shifted_pos - self.smooth_length),
                        str(shifted_pos + self.smooth_length), 
                        "N", 
                        "0",
                        strand]) + "\n")

            # For chroms with no reads, the 'read' variable will not be bound.

            _LOGGER.debug("Total number of reads processed for '{}': {}".format(chrom, nreads))
            import time
            if self.exactbw:
                # Now we have to close the stream to tell the sort process that
                # no more data is incoming. After this, we cannot use
                # `communicate()` because it returns 'ValueError: I/O operation
                # on closed file', because communicate() tries to send
                # information via stdin, which we just closed. So we instead
                # rely on wait, despite the warnings listed it
                # but actually the smugwig process is the one 
                exact_writer.close()

                # I believe we need to call communicate() only on the *last*
                # listening process, (which will wait for all preceding pipes)
                if self.compress_output:
                    wigToBigWig_process.communicate()
                else:
                    smugwig_process.communicate()
                # sort_process.wait()  # Here is where it closes out the file incomplete.
                # sort_process.communicate()

            if self.bedout:
                bedOut.close()

            if self.smoothbw:
                smooth_writer.close()
                # _LOGGER.debug("Encoding smooth bigwig for " + chrom +
                              # " (last read position:" + str(read.pos) + ")...")
                if self.compress_output:                              
                    wigToBigWig_processSm.communicate()
                else:
                    smugwig_process_smooth.communicate()

        except StopIteration as e:
            print("StopIteration error for chrom ", chrom, ": ", e)
            raise e

        return chrom


    def combine(self, good_chromosomes):
        """
        After running the process in parallel, this 'reduce' step will simply
        merge all the temporary files into one, and rename it to the output
        file name.
        """
        if not good_chromosomes:
            _LOGGER.info("No successful chromosomes, so no combining.")
            return
        elif len(good_chromosomes) == 1:
            subprocess.call(["mv", self._tempf(good_chromosomes[0]) +
                             ".bw", self.exactbw])
            if self.smoothbw:
                subprocess.call(["mv", self._tempf(good_chromosomes[0]) +
                                 "_smooth.bw", self.smoothbw])

        else:
            if self.exactbw:
                _LOGGER.info("Merging {} files into output file: '{}'".
                      format(len(good_chromosomes), self.exactbw))
                temp_files = [self._tempf(chrom) + ".bw" for chrom in good_chromosomes]
                if self.compress_output:
                    cmd = "bigWigCat {} {}".format(self.exactbw, " ".join(temp_files))
                else:
                    cmd = "cat {} > {}".format(" ".join(temp_files), self.exactbw)
                _LOGGER.debug(cmd)
                p = subprocess.call(cmd, shell=True)

            if self.smoothbw:
                _LOGGER.info("Merging {} files into output file: '{}'".
                      format(len(good_chromosomes), self.smoothbw))
                temp_files = [self._tempf(chrom) + "_smooth.bw" for chrom in good_chromosomes]
                cmd = "bigWigCat " + self.smoothbw + " " + " ".join(temp_files)
                _LOGGER.debug(cmd)
                p = subprocess.call(['bigWigCat', self.smoothbw] + temp_files)

            if self.bedout:
                # root, ext = os.path.splitext(self.exactbw)
                temp_files = [self._tempf(chrom) + ".bed" for chrom in good_chromosomes]
                cmd = "cat " + " ".join(temp_files) + " > " + self.bedout
                _LOGGER.debug(cmd)
                p = subprocess.call(cmd, shell=True)


def write_queue(queue, writer, id=""):
    """
    A function to process a particular queue. It does nothing more than read
    items from the queue, and write them to some output writer (and handles
    flushing, closing, etc).

    param queue: Queue with lines to grab and print.
    param writer: a file (or fifo) to write the lines to
    param id: some identifier for the queue
    """
    while True:
        line_to_write = queue.get()
        # print(id, line_to_write)
        if line_to_write is None:
            print("Queue {} finished. Closing fifo...".format(id))
            # Now we have to tell the sort process that nothing else is coming
            # from this fifo.
            writer.flush()
            writer.close()
            break
        writer.write(line_to_write)
        writer.flush()



def build_parser(cmdl):
    parser = ArgumentParser(description='Smugwig sam/bam reads processor.')
    parser.add_argument('-i', '--infile', dest='infile',
        help="Input file (in bam or sam format)", required=True)
    parser.add_argument('-c', '--chrom-sizes-file',
        help="Chromosome sizes file", required=True)
    parser.add_argument('-s', '--summary-file',
        help="Summary file")
    parser.add_argument("--step-type", default="variable", choices=["fixed", "variable"],
        help="Choose fixed or variable (default) step format.")
    parser.add_argument('-o', '--exactbw', dest='exactbw', default=None,
        help="Output filename for exact bigwig. Default: None")
    parser.add_argument('-w', '--smoothbw', dest='smoothbw', default=None,
        help="Output filename for smooth bigwig. Default: None")
    parser.add_argument('-r', '--step-size', default=5,
        help="Step size for smooth tracks. Default: 5")
    parser.add_argument('-b', '--bedout', default=None,
        help="Output filename for bed file. Default: None")
    parser.add_argument('-l', '--smooth-length',
        help="Smooth length for bed file", default=25, type=int)
    parser.add_argument('-d', '--tail-edge', action='store_true', default=False,
        help="Output the 3' end of the sequence read. Default: False")
    parser.add_argument('-m', '--mode', dest='mode', default=None, choices=MODES,
        help="Choose DNase or ATAC mode (adjusts shift parameters)")
    parser.add_argument('-t', '--limit', dest='limit',
        help="Limit to these chromosomes", nargs = "+", default=None)
    parser.add_argument('-p', '--cores', dest='cores',
        help="Number of cores to use", default=2, type=int)
    parser.add_argument('-x', '--suppress-compress-output', default=False, action="store_true",
        help="Summary file")    
    parser.add_argument('-e', '--temp-parent',
        default="",#os.getcwd(),
        help="Temporary file parent location. By default it will use the working"
        " directory, but you can place this elsewhere if you'd like."
        " The actual folder will be based on the exactbw filename.")
    parser.add_argument('--retain-temp', action='store_true', default=False,
        help="Retain temporary files? Default: False")

    parser = logmuse.add_logging_options(parser)
    return parser


def main():
    parser = build_parser(sys.argv[1:])
    args = parser.parse_args()
    if not (args.exactbw or args.smoothbw):
        parser.error('No output requested, use --exactbw and/or --smoothbw')
    global _LOGGER
    _LOGGER = logmuse.logger_via_cli(args)
    if args.mode == "dnase":
        shift_factor = {"+":1, "-":0}  # DNase
    elif args.mode == "atac":
        shift_factor = {"+":4, "-":-5}  # ATAC
    else:
        shift_factor = {"+":0, "-":0}

    _LOGGER.info("Using shift factor: {}".format(shift_factor))
    ct = CutTracer( reads_filename=args.infile,
                    chrom_sizes_file=args.chrom_sizes_file,
                    summary_filename=args.summary_file,
                    variable_step=args.step_type == "variable",
                    exactbw=args.exactbw,
                    smoothbw=args.smoothbw,
                    step_size=args.step_size,
                    bedout=args.bedout,
                    shift_factor=shift_factor,
                    smooth_length=args.smooth_length,
                    tail_edge=args.tail_edge,
                    suppress_compress_output=args.suppress_compress_output,
                    limit=args.limit,
                    nProc=args.cores,
                    temp_parent=args.temp_parent,
                    retain_temp=args.retain_temp,
                    verbosity=args.verbosity)

    ct.register_files()
    good_chromosomes = ct.run()
    
    _LOGGER.info("Reduce step (merge files)...")
    ct.combine(good_chromosomes)

if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        _LOGGER.error("Program canceled by user!")
        sys.exit(1)
