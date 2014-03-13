import sys, os, yaml, glob
import subprocess
import pandas as pd
from matplotlib import pyplot as plt
from de_novo_scilife import common


def run(global_config, sample_config):
    sorted_libraries_by_insert = common._sort_libraries_by_insert(sample_config)
    if "tools" in sample_config:
        """If so, execute them one after the other in the specified order (might not work)"""
        for command in sample_config["tools"]:
            """with this I pick up at run time the correct function in the current module"""
            command_fn    = getattr(sys.modules[__name__] , "_run_{}".format(command))
            """Update sample config, each command return sample_config and if necessary it modifies it"""
            sample_config = command_fn(global_config, sample_config, sorted_libraries_by_insert)
    else:
        #run default pipeline for QC
        sample_config = _run_fastqc(global_config, sample_config, sorted_libraries_by_insert)
        sample_config = _run_abyss(global_config, sample_config, sorted_libraries_by_insert)
        sample_config = _run_trimmomatic(global_config, sample_config, sorted_libraries_by_insert)



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
        folder_output_name = os.path.join(FastqcFolder, os.path.basename(read1).split(".fastq.gz")[0])
        if not common.check_dryrun(sample_config) and not os.path.exists("{}_fastqc.zip".format(folder_output_name)):
            fastq_stdOut = open(os.path.join(FastqcFolder , "{}_fastqc.stdout".format(library)), "a")
            fastq_stdErr = open(os.path.join(FastqcFolder , "{}_fastqc.stderr".format(library)), "a")
            subprocess.call(command, stdout=fastq_stdOut, stderr=fastq_stdErr)
    sample_config["fastqc"] = FastqcFolder
    return sample_config


def _run_abyss(global_config, sample_config, sorted_libraries_by_insert):
    mainDir = os.getcwd()
    ABySS_Kmer_Folder = os.path.join(os.getcwd(), "abyss_kmer")
    if "kmer" not in sample_config:
        sys.exit("error in _run_abyss QCcontrol: kmer must be present in sample_config.yaml")
    
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
    if not common.check_dryrun(sample_config) and not os.path.exists("histogram.hist"):
        ABySS_Kmer_stdOut = open("ABySS_Kmer_Folder.stdOut", "a")
        ABySS_Kmer_stdErr = open("ABySS_Kmer_Folder.stdErr", "a")
        returnValue = subprocess.call(command, shell=True, stdout=ABySS_Kmer_stdOut, stderr=ABySS_Kmer_stdErr)
        if returnValue > 0:
            print "ABySS kmer plotting failed: unkwnown reason"
        else :
            subprocess.call(("rm", "preUnitgs.fa"))
            _plotKmerPlot(1,200, kmer, "kmer_coverage_1_200.png")
            _plotKmerPlot(1,500, kmer, "kmer_coverage_1_500.png")
            _plotKmerPlot(15,200, kmer, "kmer_coverage_15_200.png")
            _plotKmerPlot(15,500, kmer, "kmer_coverage_15_500.png")

    os.chdir("..")
    sample_config["abyss"] = ABySS_Kmer_Folder
    return sample_config

def _plotKmerPlot(min_limit, max_limit,kmer, output_name):
    Kmer_histogram = pd.io.parsers.read_csv("histogram.hist", sep='\t', header=None)
    Kmer_coverage  = Kmer_histogram[Kmer_histogram.columns[0]].tolist()
    Kmer_count     = Kmer_histogram[Kmer_histogram.columns[1]].tolist()
    Kmer_freq      = [Kmer_coverage[i]*Kmer_count[i] for i in range(len(Kmer_coverage))]
    kmer_freq_peak = Kmer_freq.index(max(Kmer_freq[min_limit:max_limit]))	#coverage peak, disregarding initial peak
    kmer_freq_peak_value=max(Kmer_freq[min_limit:max_limit])
    
    xmax = max_limit
    ymax = kmer_freq_peak_value + (kmer_freq_peak_value*0.30)
    
    plt.plot(Kmer_coverage, Kmer_freq)
    plt.title("K-mer length = {}".format(kmer))
    plt.xlim((0,xmax))
    plt.ylim((0,ymax))
    plt.vlines(kmer_freq_peak, 0, kmer_freq_peak_value, colors='r', linestyles='--')
    plt.text(kmer_freq_peak, kmer_freq_peak_value+2000, str(kmer_freq_peak))
    plotname = "{}".format(output_name)
    plt.savefig(plotname)
    plt.clf()
    return 0



def _run_trimmomatic(global_config, sample_config, sorted_libraries_by_insert):
    print "running trimmomatic ..."
    
    program        = global_config["Tools"]["trimmomatic"]["bin"]
    program_folder = os.path.dirname(program)
    adapterFile    = sample_config["adapters"]
    
    if not os.path.exists(adapterFile):
        print "Trimmomatic cannot be run as adapter file is not specified or points to unknown position: {}".format(adapterFile)
        return sample_config

    rootDir       = os.getcwd()
    trimmomaticDir = os.path.join(rootDir, "Trimmomatic")
    if not os.path.exists(trimmomaticDir):
        os.makedirs(trimmomaticDir)
    else:
        print "trimmomatic assumed already run"
        return sample_config #

    os.chdir(trimmomaticDir)
    #now I am in running dir, I need to process one by one the libraries
    threads = 8
    if "threads" in sample_config:
        threads = sample_config["threads"]

    for library, libraryInfo in sorted_libraries_by_insert:
        read1=libraryInfo["pair1"]
        read2=libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        if orientation=="innie":
            if read2 is not None:
                read1_baseName = os.path.split(read1)[1].split(".")[0]
                read2_baseName = os.path.split(read2)[1].split(".")[0]
                output_read1_pair = os.path.join(trimmomaticDir,  "{}.fastq.gz".format(read1_baseName))
                output_read1_sing = os.path.join(trimmomaticDir, "{}_u.fastq.gz".format(read1_baseName))
                output_read2_pair = os.path.join(trimmomaticDir, "{}.fastq.gz".format(read2_baseName))
                output_read2_sing = os.path.join(trimmomaticDir,  "{}_u.fastq.gz".format(read2_baseName))
                command = ["java",  "-jar", program, "PE", "-threads", "{}".format(threads),  "-phred33",  read1, read2,  output_read1_pair ,output_read1_sing , output_read2_pair, output_read2_sing ,"ILLUMINACLIP:{}:2:30:10".format(adapterFile), "MINLEN:30" ]
                #print ' '.join(map(str,command))
                stdOut = open("{}_trimmomatic.stdOut".format(read1_baseName), "w")
                stdErr = open("{}_trimmomatic.stdErr".format(read1_baseName), "w")
                
                returnValue = subprocess.call(command, stdout=stdOut, stderr=stdErr) # run the program
                if returnValue != 0:
                    print "error while running command: {}".format(command)
                else:
                    libraryInfo["pair1"] = output_read1_pair
                    libraryInfo["pair2"] = output_read2_pair

    os.chdir(rootDir)
    print "trimmomatic done"
    print sample_config
    return sample_config


