from __future__ import absolute_import
from __future__ import print_function
import sys, os, yaml, glob
import subprocess
import pandas as pd
import gzip
import re
import string
import shutil
import numpy as np
from matplotlib import pyplot as plt
from nougat import common, align
from nougat.pdf.peakdetect import peakdet 

def run(global_config, sample_config):
    sorted_libraries_by_insert = common._sort_libraries_by_insert(
            sample_config)
    sample_config["commands"] = ""
    if "tools" in sample_config:
        """If so, execute them one after the other in the specified order
        (might not work)"""
        for command in sample_config["tools"]:
            """with this I pick up at run time the correct function in the
            current module"""
            command_fn = getattr(sys.modules[__name__],
                    "_run_{}".format(command))
            """Update sample config, each command return sample_config and if
            necessary it modifies it"""
            sample_config = command_fn(global_config, sample_config,
                    sorted_libraries_by_insert)
    else:
        #run default pipeline for QC
        sample_config = _run_trimmomatic(global_config, sample_config,
                sorted_libraries_by_insert)
        sample_config = _run_fastqc(global_config, sample_config,
                sorted_libraries_by_insert)
        sample_config = _run_abyss(global_config, sample_config,
                sorted_libraries_by_insert)
    with open("{}.nougat".format(sample_config.get("output", "sample")), "w") as f:
        yaml.dump(sample_config, f)


def _run_align(global_config, sample_config,sorted_libraries_by_insert):

    if "reference" not in sample_config:
        print("reference sequence not provided, skypping alignment step.",
        "Please provide a reference if you are intrested in aligning the reads",
        "against a reference")
        return sample_config
    if not os.path.exists("alignments"):
        os.makedirs("alignments")
    os.chdir("alignments")

    sorted_libraries_by_insert =  align._align_reads(global_config,
            sample_config,  sorted_libraries_by_insert) # align reads
    sorted_alignments_by_insert = align._merge_bam_files(global_config,
            sample_config, sorted_libraries_by_insert) # merge alignments
    sorted_alignments_by_insert = align.picard_CGbias(global_config,
            sample_config,sorted_alignments_by_insert) # compute picard stats
    sorted_alignments_by_insert = align.picard_collectInsertSizeMetrics(
            global_config, sample_config,sorted_alignments_by_insert)
    sorted_alignments_by_insert = align.picard_markDuplicates(global_config,
            sample_config,sorted_alignments_by_insert)

    os.chdir("..")
    sample_config["alignments"] = sorted_alignments_by_insert

    return sample_config


def _run_fastqc(global_config, sample_config, sorted_libraries_by_insert):
    mainDir = os.getcwd()
    FastqcFolder = os.path.join(os.getcwd(), "fastqc")
    if not os.path.exists(FastqcFolder):
        os.makedirs(FastqcFolder)

    program=global_config["Tools"]["fastqc"]["bin"]
    program_options=global_config["Tools"]["fastqc"]["options"]
    for library, libraryInfo in sorted_libraries_by_insert:
        command = [program]
        for option in program_options:
            command.append(option)
        read1=libraryInfo["pair1"]
        read2=libraryInfo["pair2"]
        command.append(read1)
        if read2 is not None:
            command.append(read2)
        common.print_command(command)
        sample_config["commands"] += "\n" + common.get_command_str(command)
        folder_output_name = os.path.join(FastqcFolder,
                os.path.basename(read1).split(".fastq.gz")[0])
        if not common.check_dryrun(sample_config) and not \
                os.path.exists("{}_fastqc.zip".format(folder_output_name)):
            fastq_stdOut = open(os.path.join(FastqcFolder,
                    "{}_fastqc.stdout".format(library)), "a")
            fastq_stdErr = open(os.path.join(FastqcFolder,
                    "{}_fastqc.stderr".format(library)), "a")
            subprocess.call(command, stdout=fastq_stdOut, stderr=fastq_stdErr)
    sample_config["fastqc"] = FastqcFolder
    return sample_config


def _run_abyss(global_config, sample_config, sorted_libraries_by_insert):
    mainDir = os.getcwd()
    ABySS_Kmer_Folder = os.path.join(os.getcwd(), "abyss_kmer")
    if "kmer" not in sample_config:
        sys.exit("error in _run_abyss QCcontrol: kmer must be present in \
                sample_config.yaml")

    kmer = sample_config["kmer"]
    if not os.path.exists(ABySS_Kmer_Folder):
        os.makedirs(ABySS_Kmer_Folder)

    os.chdir(ABySS_Kmer_Folder)

    program = global_config["Tools"]["abyss"]["bin"]
    program = os.path.join(os.path.dirname(program), "ABYSS-P")
    program_options=global_config["Tools"]["abyss"]["options"]
    if "abyss" in sample_config:
        program_options=sample_config["abyss"]

    threads = 16 # default for UPPMAX
    if "threads" in sample_config :
        threads = sample_config["threads"]

    command = "mpirun -np {} {} ".format(threads, program)
    command += "-k {} ".format(kmer)
    command += "--coverage-hist=histogram.hist -o preUnitgs.fa"
    for library, libraryInfo in sorted_libraries_by_insert:
        read1=libraryInfo["pair1"]
        read2=libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        if orientation=="innie" or orientation=="outtie":
            command += " {} ".format(read1)
            if read2 is not None:
                command += " {} ".format(read2)
        if orientation == "none":
            command += " {} ".format(read1)

    common.print_command(command)
    sample_config["commands"] += "\n" + common.get_command_str(command)

    if not common.check_dryrun(sample_config) and not \
            os.path.exists("histogram.hist"):
        ABySS_Kmer_stdOut = open("ABySS_Kmer_Folder.stdOut", "a")
        ABySS_Kmer_stdErr = open("ABySS_Kmer_Folder.stdErr", "a")
        returnValue = subprocess.call(command, shell=True, \
                stdout=ABySS_Kmer_stdOut, stderr=ABySS_Kmer_stdErr)
        if returnValue > 0:
            print("ABySS kmer plotting failed: unkwnown reason")
        else :
            subprocess.call(("rm", "preUnitgs.fa"))
            _plotKmerFixed(1,200, kmer, "kmer_coverage_1_200.png")
            _plotKmerFixed(1,500, kmer, "kmer_coverage_1_500.png")
            _plotKmerFixed(15,200, kmer, "kmer_coverage_15_200.png")
            _plotKmerFixed(15,500, kmer, "kmer_coverage_15_500.png")
            _plotKmer(kmer, "kmer_coverage.png")

    os.chdir("..")
    sample_config["abyss"] = ABySS_Kmer_Folder
    return sample_config


def _plotKmer(kmer, output_name):
    """Kmer abundance as a single plot, suitable for the report
    """
    Kmer_histogram = pd.io.parsers.read_csv("histogram.hist", sep='\t',
            header=None)
    Kmer_coverage = Kmer_histogram[Kmer_histogram.columns[0]].tolist()
    Kmer_count = Kmer_histogram[Kmer_histogram.columns[1]].tolist()

    # Not interested in coverage > 5000
    kcov = [c for c in Kmer_coverage if c <= 5000]
    # Multiply the abundance by a gradient
    kcount_gradient = [kcov[i] * Kmer_count[i] for i in range(len(kcov))]

    # Lazily drift towards the most spacious area under the curve
    # using divide and conquer.
    def get_bisect(chunk):
        left = chunk[:int(len(chunk)/2)]
        right = chunk[int(len(chunk)/2):]
        lweight = sum(map(lambda x: x[0] * x[1], left)) / len(left)
        rweight = sum(map(lambda x: x[0] * x[1], right)) / len(right)
        if lweight > rweight:
            return left
        else:
            return right

    # Perform six bisections
    cov_count = list(zip(kcov, kcount_gradient))
    for i in range(0,6):
        try:
            cov_count = get_bisect(cov_count)
        except ZeroDivisionError: # Already at the leftmost position
            pass
    xmax = cov_count[-1][0]

    # We could always use more space
    xmax = int(xmax * 1.3)
    ymax = max(kcount_gradient)

    # Find the highest peak x > 1. Works 70% of the time.
    maxtab, mintab = peakdet(kcount_gradient, 100000.0)
    first_peak = list(np.array(maxtab)[:,0])[0]
    # Discard x = 0 peak
    if first_peak == 0 and maxtab.size > 2:
        maxtab = np.delete(maxtab, 0, 0)
    peak = np.argmax(maxtab, axis=0)[1]
    peak = maxtab[peak][0]

    plt.xlim((0, xmax))
    plt.ylim((0, ymax))
    plt.plot(kcov, kcount_gradient)

    plt.vlines(peak, 1, kcount_gradient[peak], colors='r',
        linestyles='--')
    plt.text(peak, kcount_gradient[peak], str(peak))

    plotname = "{}".format(output_name)
    plt.savefig(plotname)
    plt.clf()
    return 0


def _plotKmerFixed(min_limit, max_limit, kmer, output_name):
    """Old kmerplot, kept just in case...
    """
    Kmer_histogram = pd.io.parsers.read_csv("histogram.hist", sep='\t',
            header=None)
    Kmer_coverage = Kmer_histogram[Kmer_histogram.columns[0]].tolist()
    Kmer_count = Kmer_histogram[Kmer_histogram.columns[1]].tolist()
    Kmer_freq = [Kmer_coverage[i]*Kmer_count[i] for i in \
            range(len(Kmer_coverage))]
    #coverage peak, disregarding initial peak
    kmer_freq_peak = Kmer_freq.index(max(Kmer_freq[min_limit:max_limit]))
    kmer_freq_peak_value=max(Kmer_freq[min_limit:max_limit])

    xmax = max_limit
    ymax = kmer_freq_peak_value + (kmer_freq_peak_value*0.30)

    plt.plot(Kmer_coverage, Kmer_freq)
    plt.title("K-mer length = {}".format(kmer))
    plt.xlim((0,xmax))
    plt.ylim((0,ymax))
    plt.vlines(kmer_freq_peak, 0, kmer_freq_peak_value, colors='r',
            linestyles='--')
    plt.text(kmer_freq_peak, kmer_freq_peak_value+2000, str(kmer_freq_peak))
    plotname = "{}".format(output_name)
    plt.savefig(plotname)
    plt.clf()
    return 0


def _run_trimmomatic(global_config, sample_config, sorted_libraries_by_insert):
    program        = global_config["Tools"]["trimmomatic"]["bin"]
    program_folder = os.path.dirname(program)
    if "adapters" not in sample_config:
        sys.exit("running MP pipeline, adapters file to be used in trimming"
                "are needed for Trimmomatic. Please specify them"
                "in the sample configuration file and rerun")
    adapterFile    = sample_config["adapters"]
    if not os.path.exists(adapterFile):
        sys.exit("Trimmomatic cannot be run as adapter file is not specified"
                "or points to unknown position: {}".format(adapterFile))

    mainDirectory = os.getcwd()
    trimmomaticDir = os.path.join(mainDirectory, "Trimmomatic")
    if not os.path.exists(trimmomaticDir):
        os.makedirs(trimmomaticDir)
    os.chdir(trimmomaticDir)
    #now I am in running dir, I need to process one by one the libraries
    threads = 8
    if "threads" in sample_config:
        threads = sample_config["threads"]

    for library, libraryInfo in sorted_libraries_by_insert:
        read1=libraryInfo["pair1"]
        read2=libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        if read2 is not None:
            read1_baseName = os.path.split(read1)[1].split(".")[0]
            read2_baseName = os.path.split(read2)[1].split(".")[0]
            output_read1_pair = os.path.join(trimmomaticDir,
                    "{}.fastq.gz".format(read1_baseName))
            output_read1_sing = os.path.join(trimmomaticDir,
                    "{}_u.fastq.gz".format(read1_baseName))
            output_read2_pair = os.path.join(trimmomaticDir,
                    "{}.fastq.gz".format(read2_baseName))
            output_read2_sing = os.path.join(trimmomaticDir,
                    "{}_u.fastq.gz".format(read2_baseName))
            command = ["java",  "-jar", program, "PE", "-threads",
                    "{}".format(threads),  "-phred33",  read1, read2,
                    output_read1_pair, output_read1_sing, output_read2_pair,
                    output_read2_sing,
                    "ILLUMINACLIP:{}:2:30:10".format(adapterFile),
                    "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15",
                    "MINLEN:30"]
            common.print_command(command)
            sample_config["commands"] += "\n" + common.get_command_str(command)

            # do not execute is files have been already gennerated
            if not common.check_dryrun(sample_config) and not \
                    os.path.exists(output_read1_pair):
                stdOut = open("{}_trimmomatic.stdOut".format(read1_baseName),
                        "w")
                stdErr = open("{}_trimmomatic.stdErr".format(read1_baseName),
                        "w")
                returnValue = subprocess.call(command, stdout=stdOut,
                        stderr=stdErr) # run the program
                if returnValue != 0:
                    print("error while running command: {}".format(command))
            libraryInfo["pair1"] = output_read1_pair
            libraryInfo["pair2"] = output_read2_pair
            libraryInfo["trimmomatic"] = os.path.join(trimmomaticDir,
                    "{}_trimmomatic.stdErr".format(read1_baseName))
    os.chdir(mainDirectory)
    return sample_config


def _kmergenie_plot(hist_file):
    """Kmergenie outputs pdf plots. We want pngs without resorting to \
    imagemagick
    TODO: Abstract this to a common plotting function if possible"""
    kgenie_hist = pd.io.parsers.read_csv(hist_file, sep=" ", header=0)
    kmer_lengths = kgenie_hist[kgenie_hist.columns[0]].tolist()
    genomic_kmers = kgenie_hist[kgenie_hist.columns[1]].tolist()
    peak_value = max(genomic_kmers)
    peak_idx = genomic_kmers.index(peak_value)
    best_k = kmer_lengths[peak_idx]

    plt.plot(kmer_lengths, genomic_kmers)
    plt.title("Best K-mer length: {}".format(best_k))
    plt.xlabel("K-mer size")
    plt.ylabel("Number of genomic k-mers")
    y_margin = (min(genomic_kmers) + peak_value) / 2 * 0.01
    y_min = min(genomic_kmers) - y_margin
    y_max = peak_value + y_margin
    plt.ylim(y_min, y_max)
    plt.xlim(min(kmer_lengths) - 5, max(kmer_lengths) + 5)
    plt.vlines(best_k, 0, peak_value, colors = "r", linestyles='--')

    plt.savefig(hist_file + ".png")
    plt.clf()


def _run_kmergenie(global_config, sample_config, sorted_libraries_by_insert):
    """Runs kmergenie to establish a recommended kmer size for assembly"""

    maindir = os.getcwd()
    kmerdir = os.path.join(maindir, "kmergenie")
    if not os.path.exists(kmerdir):
        os.makedirs(kmerdir)
    os.chdir(kmerdir)

    #Write a list of input fastq files for kmergenie
    kmer_input = os.path.join(kmerdir,
            "{}kmerinput.txt".format(sample_config.get("output","")))

    program = global_config["Tools"]["kmergenie"]["bin"]
    program_options=global_config["Tools"]["kmergenie"]["options"]
    # Could be useful to add --diploid if sample is highly heterozygous
    if "kmergenie" in sample_config:
        program_options=sample_config["kmergenie"]

    threads = "" # Kmergenie will spawn number_of_cores - 1 threads by default
    if "threads" in sample_config :
        threads = sample_config["threads"]

    cmd_list = [program, kmer_input]
    for option in filter(None, program_options):
        cmd_list.append(option)
    if threads:
        cmd_list.append("-t {}".format(threads))
    command = " ".join(cmd_list)
    common.print_command(command)
    sample_config["commands"] += "\n" + common.get_command_str(command)


    if not common.check_dryrun(sample_config):
        with open(kmer_input, "w") as f:
            for lib, lib_info in sorted_libraries_by_insert:
                f.write(lib_info["pair1"] + "\n")
                f.write(lib_info["pair2"] + "\n")

        stdOut = open("kmergenie.stdOut", "w")
        stdErr = open("kmergenie.stdErr", "w")
        returnValue = subprocess.call(cmd_list, stdout=stdOut, stderr=stdErr)
        if returnValue != 0:
            print("error while running command: {}".format(command))
        else:
            _kmergenie_plot("histograms.dat")
    sample_config["kmergenie"] = kmerdir
    os.chdir(maindir)
    return sample_config


