import sys, os, yaml, glob
import subprocess
import pandas as pd
from matplotlib import pyplot as plt
import align
import common



def run(global_config, sample_config):
    sorted_libraries_by_insert = sorted(sample_config["libraries"].iteritems(), key=lambda (k,v): v["insert"])
    if os.path.exists(os.path.join(os.getcwd(), "DATA")):
        sorted_libraries_by_insert = common.update_sample_config(sorted_libraries_by_insert)
    else:
        sorted_libraries_by_insert = common.prepare_folder_structure(sorted_libraries_by_insert)

    print "remove adator"
    sample_config = _run_trimmomatic(global_config, sample_config,sorted_libraries_by_insert)
    print "align sequences"
    sorted_libraries_by_insert = sorted(sample_config["libraries"].iteritems(), key=lambda (k,v): v["insert"]) # recompute sorted libraries by insert
    ## TO DO: create a reference
    

    if not os.path.exists("alignments"):
        os.makedirs("alignments")
    os.chdir("alignments")
    print sorted_libraries_by_insert
    sorted_libraries_by_insert =  align._align_reads(global_config, sample_config,  sorted_libraries_by_insert) # align reads
    print sorted_libraries_by_insert
    sorted_alignments_by_insert = align._merge_bam_files(global_config, sample_config, sorted_libraries_by_insert) # merge alignments
    print sorted_libraries_by_insert
    sorted_alignments_by_insert = align.picard_CGbias(global_config, sample_config,sorted_alignments_by_insert)
    print sorted_libraries_by_insert
    sorted_alignments_by_insert = align.picard_collectInsertSizeMetrics(global_config, sample_config,sorted_alignments_by_insert)
    print sorted_libraries_by_insert
    sorted_alignments_by_insert = align.picard_markDuplicates(global_config, sample_config,sorted_alignments_by_insert)
    print sorted_libraries_by_insert

    os.chdir("..")
    return 0


def _run_trimmomatic(global_config, sample_config, sorted_libraries_by_insert):
    print "running trimmomatic ..."
    mainDir = os.getcwd()
    program        = global_config["Tools"]["trimmomatic"]["bin"]
    program_folder = os.path.dirname(program)
    adaprterFile = os.path.join(program_folder, "adapters", "NexteraMP.fa")
    if not os.path.exists(adaprterFile):
        print "file {} does not exists: it must be present: you must provide it!!!".format(adaprterFile)
        return sample_config

    for library, libraryInfo in sorted_libraries_by_insert:
        read1=libraryInfo["pair1"]
        read2=libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert = libraryInfo["insert"]
        workingDir = "{}_MP_{}".format(library,insert)
        if not os.path.exists(workingDir):
            os.makedirs(workingDir)

        os.chdir(workingDir)
        currentDir = os.getcwd()
        output_read1_pair = os.path.join(currentDir, "{}_MP_1.fastq.gz".format(library))
        output_read1_sing = os.path.join(currentDir, "{}_MP_u_1.fastq.gz".format(library))
        output_read2_pair = os.path.join(currentDir, "{}_MP_2.fastq.gz".format(library))
        output_read2_sing = os.path.join(currentDir, "{}_MP_u_2.fastq.gz".format(library))

        if os.path.exists(output_read1_pair):
            print "library {} already computed: skyp this".format(library)
            libraryInfo["pair1"] = output_read1_pair
            libraryInfo["pair2"] = output_read2_pair
            os.chdir("..")
            continue

        threads = 8
        if "threads" in sample_config:
            threads = sample_config["threads"]
        
        command = ["java",  "-jar", program, "PE", "-threads", "{}".format(threads),  "-phred33",  read1, read2,  output_read1_pair ,output_read1_sing , output_read2_pair, output_read2_sing ,"ILLUMINACLIP:{}:2:30:10".format(adaprterFile), "MINLEN:30" ]
        #print ' '.join(map(str,command))
        stdOut = open("trimmomatic.stdOut", "w")
        stdErr = open("trimmomatic.stdErr", "w")
        returnValue = subprocess.call(command, stdout=stdOut, stderr=stdErr)
        returnValue = 0
        if returnValue != 0:
            print "error while running command: {}".format(command)
        else:
            libraryInfo["pair1"] = output_read1_pair
            libraryInfo["pair2"] = output_read2_pair


        os.chdir("..")
    
    
    print "trimmomatic done"
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




