import sys, os, yaml, glob
import subprocess
import pandas as pd
import align
import common



def run(global_config, sample_config):
    sorted_libraries_by_insert = common._sort_libraries_by_insert(sample_config)
    print sorted_libraries_by_insert
    run_default(global_config, sample_config,sorted_libraries_by_insert)



def run_default(global_config, sample_config,sorted_libraries_by_insert):
    print "running the entire pipeline"
    print " remove adator"
    sample_config = _run_trimmomatic(global_config, sample_config,sorted_libraries_by_insert, "01_trimmomatic")

    sorted_libraries_by_insert = sorted(sample_config["libraries"].iteritems(), key=lambda (k,v): v["insert"]) # recompute sorted libraries by insert
    
    print " align sequences"
    if not os.path.exists("02_alignments"):
        os.makedirs("02_alignments")
    os.chdir("02_alignments")
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


def _run_trimmomatic(global_config, sample_config, sorted_libraries_by_insert, step):
    print "     running trimmomatic ..."
    mainDir = os.getcwd()
    program        = global_config["Tools"]["trimmomatic"]["bin"]
    program_folder = os.path.dirname(program)
    adapterFile = os.path.join(program_folder, "adapters", "NexteraMP.fa")
    if not os.path.exists(adapterFile):
        print "file {} does not exists: it must be present: you must provide it!!!".format(adapterFile)
        return sample_config

    runningDir = os.path.join(mainDir, step)
    if not os.path.exists(runningDir):
        os.makedirs(runningDir)
    os.chdir(runningDir)
    #now I am in running dir, I need to process one by one the libraries

    for library, libraryInfo in sorted_libraries_by_insert:
        read1=libraryInfo["pair1"]
        read2=libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert = libraryInfo["insert"]
        workingOnLibraryDir = os.path.join(runningDir, "{}_MP_{}".format(library,insert))
        if not os.path.exists(workingOnLibraryDir):
            os.makedirs(workingOnLibraryDir)
        os.chdir(workingOnLibraryDir)
        output_read1_pair = os.path.join(workingOnLibraryDir, "{}_MP_1.fastq.gz".format(library))
        output_read1_sing = os.path.join(workingOnLibraryDir, "{}_MP_u_1.fastq.gz".format(library))
        output_read2_pair = os.path.join(workingOnLibraryDir, "{}_MP_2.fastq.gz".format(library))
        output_read2_sing = os.path.join(workingOnLibraryDir, "{}_MP_u_2.fastq.gz".format(library))

        if os.path.exists(output_read1_pair):
            print "library {} already computed: skyp this".format(library)
            libraryInfo["pair1"] = output_read1_pair
            libraryInfo["pair2"] = output_read2_pair
            os.chdir("..")
            continue

        threads = 8
        if "threads" in sample_config:
            threads = sample_config["threads"]
        
        command = ["java",  "-jar", program, "PE", "-threads", "{}".format(threads),  "-phred33",  read1, read2,  output_read1_pair ,output_read1_sing , output_read2_pair, output_read2_sing ,"ILLUMINACLIP:{}:2:30:10".format(adapterFile), "MINLEN:30" ]
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

        os.chdir(runningDir)
    
    os.chdir(mainDir)
    print "trimmomatic done"
    return sample_config


