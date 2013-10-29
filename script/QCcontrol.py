import sys, os, yaml, glob
import subprocess
import pandas as pd
from matplotlib import pyplot as plt



def run(global_config, sample_config, sorted_libraries_by_insert):
    _run_fastqc(global_config, sample_config, sorted_libraries_by_insert)
    _run_abyss_kmerPlot(global_config, sample_config, sorted_libraries_by_insert)


def _run_fastqc(global_config, sample_config, sorted_libraries_by_insert):
    print "running fastqc ..."
    mainDir = os.getcwd()
    FastqcFolder = os.path.join(os.getcwd(), "fastqc")
    if not os.path.exists(FastqcFolder):
        os.makedirs(FastqcFolder)
    else:
        print "done (fastqc folder already present, assumed already run)"
        return
    fastq_stdOut = open("fastqc.stdOut", "a")
    fastq_stdErr = open("fastqc.stdErr", "a")
    program=global_config["QCcontrol"]["fastqc"]["bin"]
    program_options=global_config["QCcontrol"]["fastqc"]["options"]
    for library, libraryInfo in sorted_libraries_by_insert:
        command = [program]
        for option in program_options:
            command.append(option)
        read1=libraryInfo["pair1"]
        read2=libraryInfo["pair2"]
        command.append(read1)
        if read2 is not None:
            command.append(read2)
        print command
        subprocess.call(command, stdout=fastq_stdOut, stderr=fastq_stdErr)
    print "fastqc succesfully exectued"
    return


def _run_abyss_kmerPlot(global_config, sample_config, sorted_libraries_by_insert):
    print "I am running abyss to check kmer content"
    mainDir = os.getcwd()
    ABySS_Kmer_Folder = os.path.join(os.getcwd(), "abyss_kmer")
    if not os.path.exists(ABySS_Kmer_Folder):
        os.makedirs(ABySS_Kmer_Folder)
    else:
        print "done (ABySS_Kmer_Folder folder already present, assumed already run)"
        #return
   
    os.chdir(ABySS_Kmer_Folder)
    ABySS_Kmer_stdOut = open("ABySS_Kmer_Folder.stdOut", "a")
    ABySS_Kmer_stdErr = open("ABySS_Kmer_Folder.stdErr", "a")

    program=global_config["QCcontrol"]["abyss"]["bin"]
    program_options=global_config["QCcontrol"]["abyss"]["options"]
    command = ['mpirun', '-np', '8', program]
    for option in program_options:
            command.append(option)
    for library, libraryInfo in sorted_libraries_by_insert:
        read1=libraryInfo["pair1"]
        read2=libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        if orientation=="innie":
            command.append(read1)
            if read2 is not None:
                command.append(read2)
        if orientation == "none":
                command.append(read1)
    print command
    subprocess.call(command, stdout=ABySS_Kmer_stdOut, stderr=ABySS_Kmer_stdErr)
    subprocess.call(("rm", "preUnitgs.fa"))
    _plotKmerPlot()

    os.chdir("..")
    return

def _plotKmerPlot():
    Kmer_histogram = pd.io.parsers.read_csv("histogram.hist", sep='\t', header=None)
    Kmer_coverage  = Kmer_histogram[Kmer_histogram.columns[0]].tolist()
    Kmer_count     = Kmer_histogram[Kmer_histogram.columns[1]].tolist()
    Kmer_freq      = [Kmer_coverage[i]*Kmer_count[i] for i in range(len(Kmer_coverage))]
    kmer_freq_peak = Kmer_freq.index(max(Kmer_freq[7:]))	#coverage peak, disregarding initial peak
    kmer_freq_peak_value=max(Kmer_freq[7:])
    
    xmax = 200
    ymax = kmer_freq_peak_value+(kmer_freq_peak_value*0.30)
    
    plt.plot(Kmer_coverage, Kmer_freq)
    plt.title('K-mer length = %s' % 54)
    plt.xlim((0,xmax))
    plt.ylim((0,ymax))
    plt.vlines(kmer_freq_peak, 0, kmer_freq_peak_value, colors='r', linestyles='--')
    plt.text(kmer_freq_peak, kmer_freq_peak_value+20000, str(kmer_freq_peak))
    plotname = "kmer_coverage.png"
    plt.savefig(plotname)
    plt.clf()
    return 0


def _run_jellyfish(global_config, sample_config):
    print "Jellyfish still to be implemented"




