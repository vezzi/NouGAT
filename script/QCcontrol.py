import sys, os, yaml, glob
import subprocess
import pandas as pd
from matplotlib import pyplot as plt



def run(global_config, sample_config):
    sorted_libraries_by_insert = sorted(sample_config["libraries"].iteritems(), key=lambda (k,v): v["insert"])
    if os.path.exists(os.path.join(os.getcwd(), "DATA")):
        sorted_libraries_by_insert = update_sample_config(sorted_libraries_by_insert)
    else:
        sorted_libraries_by_insert = prepare_folder_structure(sorted_libraries_by_insert)

    if "tools" in sample_config:
        #need to follow the commands listed here
        for command in sample_config["tools"]:
            command_fn = getattr( sys.modules[__name__] , "_run_{}".format(command))
            sample_config = command_fn(global_config, sample_config, sorted_libraries_by_insert)
    else:
        #run default pipeline for QC
        sample_config = _run_fastqc(global_config, sample_config, sorted_libraries_by_insert)
        sample_config = _run_abyss(global_config, sample_config, sorted_libraries_by_insert)



def prepare_folder_structure(sorted_libraries_by_insert):
    mainDir = os.getcwd()
    DataFolder = os.path.join(os.getcwd(), "DATA")
    if os.path.exists(DataFolder):
        sys.exit("DATA dir already exists: danger to over-write data: terminate execution")
    os.makedirs(DataFolder)
    os.chdir(DataFolder)
    CurrentDir = os.getcwd()
    #now prepare softlinks to data and give to libraries human readable names
    currentLibraryNumber = 1;
    type = ["SE", "PE", "MP"]
    for library, libraryInfo in sorted_libraries_by_insert:
        pair1 = libraryInfo["pair1"]
        pair2 = libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        pair1, pair2 = createSoftLinks(pair1, pair2, orientation, type, currentLibraryNumber)
        libraryInfo["pair1"] = pair1
        libraryInfo["pair2"] = pair2
        currentLibraryNumber += 1
    os.chdir("..")
    return sorted_libraries_by_insert

def update_sample_config(sorted_libraries_by_insert):
    mainDir = os.getcwd()
    DataFolder = os.path.join(os.getcwd(), "DATA")
    if not os.path.exists(DataFolder):
        sys.exit("DATA dir does not exists: we should not be here!!!!")
    os.chdir(DataFolder)
    currentLibraryNumber = 1
    type = ["SE", "PE", "MP"]
    for library, libraryInfo in sorted_libraries_by_insert:
        pair1 = libraryInfo["pair1"]
        pair2 = libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        libraryInfo["pair1"] = _new_name(pair1, orientation, type, currentLibraryNumber, 1)
        libraryInfo["pair2"] = _new_name(pair2, orientation, type, currentLibraryNumber, 2)
        currentLibraryNumber += 1
    os.chdir("..")
    return sorted_libraries_by_insert

def createSoftLinks(pair1, pair2, orientation, type, currentLibraryNumber):
    pair1NewName = _new_name(pair1, orientation, type, currentLibraryNumber, 1)
    pair2NewName = _new_name(pair2, orientation, type, currentLibraryNumber, 2)
    os.symlink(pair1, pair1NewName)
    if pair2NewName is not None:
         os.symlink(pair2, pair2NewName)
    return pair1NewName, pair2NewName

def _new_name(oldPathName, orientation, type, currentLibraryNumber, pairNumber):
    if oldPathName is None:
        return oldPathName;
    oldName = os.path.split(oldPathName)[1]
    oldNameHead , oldNameTail = oldName.split(".",1)
    newName = "lib{}_".format(currentLibraryNumber)
    if orientation == "none":
        newName += "SE."
    elif orientation == "innie":
        newName += "PE_{}.".format(pairNumber)
    elif orientation == "outtie":
        newName += "MP_{}.".format(pairNumber)
    newName += oldNameTail
    newName = os.path.join(os.getcwd(), newName)
    return newName



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
        print command
        subprocess.call(command, stdout=fastq_stdOut, stderr=fastq_stdErr)
    print "fastqc succesfully exectued"
    return sample_config


def _run_abyss(global_config, sample_config, sorted_libraries_by_insert):
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

    program = global_config["Tools"]["abyss"]["bin"]
    program = os.path.join(os.path.dirname(program), "ABYSS-P")
    program_options=global_config["Tools"]["abyss"]["options"]
    if "abyss" in sample_config:
        program_options=sample_config["abyss"]
    
    threads = 8 # default for UPPMAX
    if "threads" in sample_config :
        threads = sample_config["threads"]

    command = "mpirun -np {} {} ".format(threads, program)
    kmer = sample_config["kmer"]
    command += "-k {} ".format(kmer)
    command += "--coverage-hist=histogram.hist -o preUnitgs.fa"
    for library, libraryInfo in sorted_libraries_by_insert:
        read1=libraryInfo["pair1"]
        read2=libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        if orientation=="innie":
            command += " {} ".format(read1)
            if read2 is not None:
                command += " {} ".format(read2)
        if orientation == "none":
            command += " {} ".format(read1)
    print command
    subprocess.call(command, shell=True, stdout=ABySS_Kmer_stdOut, stderr=ABySS_Kmer_stdErr)
    subprocess.call(("rm", "preUnitgs.fa"))
    _plotKmerPlot()

    os.chdir("..")
    return

def _plotKmerPlot():
    Kmer_histogram = pd.io.parsers.read_csv("histogram.hist", sep='\t', header=None)
    Kmer_coverage  = Kmer_histogram[Kmer_histogram.columns[0]].tolist()
    Kmer_count     = Kmer_histogram[Kmer_histogram.columns[1]].tolist()
    Kmer_freq      = [Kmer_coverage[i]*Kmer_count[i] for i in range(len(Kmer_coverage))]
    kmer_freq_peak = Kmer_freq.index(max(Kmer_freq[15:400]))	#coverage peak, disregarding initial peak
    kmer_freq_peak_value=max(Kmer_freq[15:400])
    
    xmax = 200
    ymax = kmer_freq_peak_value + (kmer_freq_peak_value*0.30)
    
    plt.plot(Kmer_coverage, Kmer_freq)
    plt.title('K-mer length = %s' % 54)
    plt.xlim((0,xmax))
    plt.ylim((0,ymax))
    plt.vlines(kmer_freq_peak, 0, kmer_freq_peak_value, colors='r', linestyles='--')
    plt.text(kmer_freq_peak, kmer_freq_peak_value+2000, str(kmer_freq_peak))
    plotname = "kmer_coverage.png"
    plt.savefig(plotname)
    plt.clf()
    return 0


def _run_jellyfish(global_config, sample_config, sorted_libraries_by_insert):
    print "Jellyfish still to be implemented"




