import sys, os, yaml, glob
import subprocess
import common
import string
import sys
import gzip
import pandas as pd



def _align_reads(global_config, sample_config, sorted_libraries_by_insert):
    reference   = sample_config["reference"]
    reference =  build_reference_bwa(global_config,reference)
    #now prepare for alignment
    for library, libraryInfo in sorted_libraries_by_insert:
        read1       = libraryInfo["pair1"]
        read2       = libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert      = libraryInfo["insert"]
        std         = libraryInfo["std"]
        threads     = 8
        if "threads" in sample_config:
            threads  = sample_config["threads"]
        libraryInfo["alignment"] = align_bwa_mem(global_config, read1, read2, reference, threads)
    return sorted_libraries_by_insert


def _merge_bam_files(global_config, sample_config, sorted_libraries_by_insert):
    BAMfiles = {};
    reference = sample_config["reference"]
    
    samtools = "samtools"
    if "samtools" in global_config["Tools"]:
        samtools = global_config["Tools"]["samtools"]["bin"]
    elif not common.which("samtools"):
        sys.exit("error while trying to run  samtools: bwa not present in the path and not in global config, please make sure to install bwa properly")

    for library, libraryInfo in sorted_libraries_by_insert:
        read1       = libraryInfo["pair1"]
        read2       = libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert      = libraryInfo["insert"]
        std         = libraryInfo["std"]
        alignment   = libraryInfo["alignment"]
        if insert not in BAMfiles:
            BAMfiles[insert] = [alignment]
        else:
            BAMfiles[insert].append(alignment)
    
    BAMfilesMerged = {}
    for insert, insertGroup in BAMfiles.iteritems():
        dir_insert = "lib_{}".format(insert)
        if not os.path.exists(dir_insert):
            os.makedirs(dir_insert)
        os.chdir(dir_insert)
        #check if file is already present
        bamMerged = "lib_{}.bam".format(insert)
        if os.path.exists(bamMerged):
            BAMfilesMerged[insert] = os.path.abspath(bamMerged)
            os.chdir("..")
            continue # nothiing to be done for this insert

        if len(insertGroup) == 1: # only one sample file for this insert length
            cl = ["ln", "-s", insertGroup[0], bamMerged]
            returnValue = subprocess.call(cl)
            if  not returnValue == 0:
                sys.exit("error, while soft linking {}".format(insertGroup[0]))
        else:
            cl = [samtools, "merge"]
            bamMerged = "lib_{}.bam".format(insert)
            cl.append(bamMerged)
            for bamfile in insertGroup:
                cl.append(bamfile)
            returnValue = subprocess.call(cl)
            if  not returnValue == 0:
                sys.exit("error, while merging files {}".format(insertGroup))
        BAMfilesMerged[insert] = os.path.abspath(bamMerged)
        os.chdir("..")
    
    sorted_alignments_by_insert = []
    for key in sorted(BAMfilesMerged.iterkeys()):
        sorted_alignments_by_insert.append([key, BAMfilesMerged[key]])
    return sorted_alignments_by_insert



def picard_CGbias(global_config, sample_config, sorted_alignments_by_insert):
    picard = "";
    if os.environ.get('PICARD_HOME'):
        picard = os.environ.get('PICARD_HOME')
    elif "picard" in global_config["Tools"]:
        picard = global_config["Tools"]["picard"]["bin"]
    for library, BAMfile in sorted_alignments_by_insert:
        working_dir = "lib_{}".format(library)
        os.chdir(working_dir)
        output_header = os.path.basename(BAMfile).split(".bam")[0]
        command= ["java", "-Xmx16g", "-XX:PermSize=2g", "-jar", os.path.join(picard, "CollectGcBiasMetrics.jar"),  "REFERENCE_SEQUENCE={}".format(sample_config["reference"]), "INPUT={}".format(BAMfile), \
        "OUTPUT={}.collectGcBias.txt".format(output_header), "CHART_OUTPUT={}.collectGcBias.pdf".format(output_header), "ASSUME_SORTED=true", "VALIDATION_STRINGENCY=LENIENT", "TMP_DIR=/scratch"]
        print command
        returnValue = 0;
        if not os.path.exists("{}.collectGcBias.pdf".format(output_header)):
            stdOut = open("collectGcBias.stdOut", "w")
            stdErr = open("collectGcBias.stdErr", "w")
            returnValue = subprocess.call(command, stdout=stdOut, stderr=stdErr)
        if not returnValue == 0:
            print "problem running collectGCBias"
        os.chdir("..")
    return sorted_alignments_by_insert



def picard_collectInsertSizeMetrics(global_config, sample_config, sorted_alignments_by_insert):
    picard = "";
    if os.environ.get('PICARD_HOME'):
        picard = os.environ.get('PICARD_HOME')
    elif "picard" in global_config["Tools"]:
        picard = global_config["Tools"]["picard"]["bin"]
    for library, BAMfile in sorted_alignments_by_insert:
        working_dir = "lib_{}".format(library)
        os.chdir(working_dir)
        output_header = os.path.basename(BAMfile).split(".bam")[0]
        histWide = library * 2
        command= ["java", "-Xmx16g", "-XX:PermSize=2g", "-jar", os.path.join(picard, "CollectInsertSizeMetrics.jar"), "INPUT={}".format(BAMfile), "MINIMUM_PCT=0",
        "HISTOGRAM_FILE={}.collectInsertSize.pdf".format(output_header), "OUTPUT={}.collectInsertSiz.txt".format(output_header), "HISTOGRAM_WIDTH={}".format(histWide), "VALIDATION_STRINGENCY=LENIENT", "TMP_DIR=/scratch"]
        print command
        returnValue = 0;
        if not os.path.exists("{}.collectInsertSize.pdf".format(output_header)):
            stdOut = open("collectInsertSize.stdOut", "w")
            stdErr = open("collectInsertSize.stdErr", "w")
            returnValue = subprocess.call(command, stdout=stdOut, stderr=stdErr)
        if not returnValue == 0:
            print "problem running CollectInsertSizeMetrics"
        os.chdir("..")
    return sorted_alignments_by_insert


def picard_markDuplicates(global_config, sample_config, sorted_alignments_by_insert):
    picard = "";
    if os.environ.get('PICARD_HOME'):
        picard = os.environ.get('PICARD_HOME')
    elif "picard" in global_config["Tools"]:
        picard = global_config["Tools"]["picard"]["bin"]
    for library, BAMfile in sorted_alignments_by_insert:
        working_dir = "lib_{}".format(library)
        os.chdir(working_dir)
        output_header = os.path.basename(BAMfile).split(".bam")[0]
        command= ["java", "-Xmx16g", "-XX:PermSize=3g", "-jar", os.path.join(picard, "MarkDuplicates.jar"), "INPUT={}".format(BAMfile), "OUTPUT={}_noDup.bam".format(output_header), \
          "METRICS_FILE={0}.markDuplicates.txt".format(output_header), "ASSUME_SORTED=true", "VALIDATION_STRINGENCY=LENIENT", "TMP_DIR=/scratch"]
        print command
        returnValue = 0;
        if not os.path.exists("{}_noDup.bam".format(output_header)):
            stdOut = open("removeDup.stdOut", "w")
            stdErr = open("removeDup.stdErr", "w")
            returnValue = subprocess.call(command, stdout=stdOut, stderr=stdErr)
        if not returnValue == 0:
            print "problem running MarkDuplicates"
        os.chdir("..")
    return sorted_alignments_by_insert


def build_reference_bwa(global_config, reference):
    #build the reference if not available
    program = "bwa"
    if "bwa" in global_config["Tools"]:
        program = global_config["Tools"]["bwa"]["bin"]
    elif not common.which("bwa"):
        sys.exit("error while trying to run  bwa index: bwa not present in the path and not in global config, please make sure to install bwa properly")
    # check if reference provided exisists
    reference = os.path.abspath(reference)
    if not os.path.exists(reference):
        sys.exit("error, reference file {} does not exists".format(reference))
    # check if bwa index already created
    current_dir           = os.getcwd()
    path_name, base_name  = os.path.split(reference)
    bwa_index_folder      = os.path.join(path_name, "bwa")

    #if needed create directory
    if not os.path.exists(bwa_index_folder):
        os.makedirs(bwa_index_folder)
    os.chdir(bwa_index_folder)

    # if needed soft link the reference
    if not os.path.exists(base_name):
        returnValue = subprocess.call(["ln", "-s", reference, base_name])
        if not returnValue == 0:
            sys.exit("error while trying to soft link reference sequence")

    reference = os.path.join(path_name, "bwa", base_name) # now I have a soflinked copy
    # now check if index alredy build or not
    if not os.path.exists("{}.bwt".format(reference)): # then create the index sequence
        bwa_stdOut = open("bwa_index.stdOut", "w")
        bwa_stdErr = open("bwa_index.stdErr", "w")
        print "creating reference"
        returnValue = subprocess.call([program, "index", reference], stdout=bwa_stdOut, stderr=bwa_stdErr)
        if  not returnValue == 0:
            sys.exit("error, while indexing reference file {} with bwa index".format(reference))

    #extra control to avoid problem with nuexpected return value
    if not os.path.exists("{}.bwt".format(reference)):
        sys.exit("bwa index failed")

    os.chdir(current_dir)

    return reference




def align_bwa_mem(global_config, read1, read2, reference, threads):
    aligner = "bwa"
    if "bwa" in global_config["Tools"]:
        aligner = global_config["Tools"]["bwa"]["bin"]
    elif not common.which("bwa"):
        sys.exit("error while trying to run  bwa mem: bwa not present in the path and not in global config, please make sure to install bwa properly")
    
    samtools = "samtools"
    if "samtools" in global_config["Tools"]:
        samtools = global_config["Tools"]["samtools"]["bin"]
    elif not common.which("samtools"):
        sys.exit("error while trying to run  samtools: bwa not present in the path and not in global config, please make sure to install bwa properly")

    # extract base name, assume the only dot is
    libraryBase = os.path.basename(read1).split(".fastq")[0]
    if read2:
        libraryBase = "{}_{}".format(libraryBase.split("_")[0], libraryBase.split("_")[1])
    if not os.path.exists(libraryBase):
        os.makedirs(libraryBase)
    os.chdir(libraryBase)
    mappingBase = "{}_to_{}".format(libraryBase, os.path.basename(reference).split(".fasta")[0])
    BAMsorted   = "{}.bam".format(mappingBase)
    BAMunsorted = "{}.unsorted.bam".format(mappingBase)
    SAMMapped   = "{}.unsorted.sam".format(mappingBase)
    if os.path.exists(os.path.abspath(BAMsorted)):
        BAMsorted = os.path.abspath(BAMsorted)
        os.chdir("..")
        return BAMsorted

    #if not os.path.exists(SAMMapped) and not os.path.exists(BAMunsorted):
    #    with open(SAMMapped, "w") as fh:
    #            stdErr = open("bwa.stdErr", "w")
    #            returnValue = subprocess.call([aligner, "mem", "-M", "-t", "{}".format(threads), reference, read1, read2], stdout=fh, stderr=stdErr)
    #            if  not returnValue == 0:
    #                sys.exit("error, while aligning reads with bwa mem")
    #
    #if not os.path.exists(BAMunsorted):
    #    with open(BAMunsorted, "w") as fh:
    #            stdErr = open("bwa.stdErr", "w")
    #            returnValue = subprocess.call([samtools, "view", "-S", "-F", "4", "-b", SAMMapped], stdout=fh, stderr=stdErr)
    #           if  not returnValue == 0:
    #                sys.exit("error, while while converting sam to bam")
    
    bwa_mem_command       = [aligner, "mem", "-M", "-t", "{}".format(threads), reference, read1, read2]
    samtools_view_command = [samtools, "view", "-b", "-S",  "-u",  "-"]

    if not os.path.exists(BAMunsorted):
        command = "{} | {} > {}".format(" ".join(bwa_mem_command), " ".join(samtools_view_command), BAMunsorted)
        bwa_stdOut       = open("bwa.stdOut", "w")
        bwa_stdErr       = open("bwa.stdErr", "w")
        subprocess.call(command, shell=True, stdout=bwa_stdOut, stderr=bwa_stdErr)

    samtools_sort_command = [samtools, "sort", "-@", "{}".format(threads), "-m" , "1G", BAMunsorted,  mappingBase]
    command = " ".join(samtools_sort_command)
    if not os.path.exists(BAMsorted):
        stdOut       = open("sam_sort.stdOut", "w")
        stdErr       = open("sam_sort.stdErr", "w")
        subprocess.call(command, shell=True, stdout=stdOut, stderr=stdErr)

    if os.path.exists(BAMsorted) and os.path.exists(BAMunsorted)  :
        subprocess.call(["rm", BAMunsorted])
    BAMsorted = os.path.abspath(BAMsorted)
    os.chdir("..")
    return BAMsorted


def compl_rev(read_file):
    originalHeader  = os.path.basename(read_file).split(".fastq")[0].split("_")
    outputFileName  = "{}rc_{}_{}.fastq.gz".format(originalHeader[0],originalHeader[1],originalHeader[2])
    outputFileName  = os.path.abspath(outputFileName)
    if os.path.exists(outputFileName):
        return outputFileName
    
    fd_in = ""
    if ".gz" in read_file:
        fd_in = gzip.open(filename=read_file, mode='rb')
    else:
        fd_in = open(filename, 'rb')


    fd_out          = gzip.open(outputFileName, 'wb')
    header = fd_in.readline().rstrip()
    while header:
        read    = fd_in.readline().rstrip()
        comment = fd_in.readline().rstrip()
        quality = fd_in.readline().rstrip()
        read    = rc(read)
        quality = quality[::-1]
        fd_out.write("{}\n".format(header))
        fd_out.write("{}\n".format(read))
        fd_out.write("{}\n".format(comment))
        fd_out.write("{}\n".format(quality))
        header  = fd_in.readline().rstrip()
    fd_in.close()
    fd_out.close()
    return outputFileName


import string

def rc(dna):
    complements = string.maketrans('acgtnACGTN', 'tgcanTGCAN')
    rcseq = dna.translate(complements)[::-1]
    return rcseq



def align_bwa(read1, read2, reference):
    aligner = "bwa"
    if "bwa" in global_config["Tools"]:
        program = global_config["Tools"]["bwa"]["bin"]
    elif not common.which("bwa"):
        sys.exit("error while trying to run  bwa aln: bwa not present in the path and not in global config, please make sure to install bwa properly")
    
    samtools = "samtools"
    if "samtools" in global_config["Tools"]:
        program = global_config["Tools"]["samtools"]["bin"]
    elif not common.which("samtools"):
        sys.exit("error while trying to run  samtools: bwa not present in the path and not in global config, please make sure to install bwa properly")
 
    libraryBase = os.path.basename(read1).split(".")[0]
    if read2:
        libraryBase = "{}_{}".format(libraryBase.split("_")[0], libraryBase.split("_")[1])
    if not os.path.exists(libraryBase):
        os.makedirs(libraryBase)
    os.chdir(libraryBase)
    mappingBase = "{}_to_{}".format(libraryBase, os.path.basename(reference).split(".fasta")[0])
    BAMsorted   = "{}.bam".format(mappingBase)
    BAMunsorted = "{}.unsorted.bam".format(mappingBase)
    if os.path.exists(os.path.abspath(BAMsorted)):
        BAMsorted = os.path.abspath(BAMsorted)
        os.chdir("..")
        return BAMsorted

    if not os.path.exists(BAMunsorted):
        bwaAln_1 = "{}_1.sai".format(mappingBase);
        bwaAln_2 = "{}_2.sai".format(mappingBase);
        if not os.path.exists("{}_1.sai".format(mappingBase)):
            with open(bwaAln_1, "w") as fh:
                returnValue = subprocess.call(["bwa", "aln", "-t", "8", reference, read1], stdout=fh)
                if  not returnValue == 0:
                    sys.exit("error, while aligning read {} against {}".format(read1, reference))
        if read2:
            if not os.path.exists("{}_2.sai".format(mappingBase)):
                with open(bwaAln_2, "w") as fh:
                    returnValue = subprocess.call(["bwa", "aln", "-t", "8", reference, read2], stdout=fh)
                    if  not returnValue == 0:
                        sys.exit("error, while aligning read {} against {}".format(read1, reference))
        with open(BAMunsorted, 'w') as fh:
            if read2:
                cl1 = ["bwa", "sampe", "-P", "-s",  reference, bwaAln_1, bwaAln_2, read1, read2]
                cl2 = ["samtools", "view", "-Shb", "-F", "4", "-"]
                p1 = subprocess.Popen(cl1, stdout=subprocess.PIPE)
                p2 = subprocess.Popen(cl2, stdin=p1.stdout, stdout=fh).communicate()
                p1.stdout.close()
            else:
                cl1 = ["bwa", "samse",  reference, bwaAln_1, read1]
                cl2 = ["samtools", "view", "-Shb", "-F", "4", "-"]
                p1 = subprocess.Popen(cl1, stdout=subprocess.PIPE)
                p2 = subprocess.Popen(cl2, stdin=p1.stdout, stdout=fh).communicate()
                p1.stdout.close()

    subprocess.call(["rm", bwaAln_1, bwaAln_2])
    BAMsortedHeader = "{}".format(mappingBase)
    subprocess.call(["samtools", "sort",  BAMunsorted, BAMsortedHeader])
    subprocess.call(["rm", BAMunsorted])
    BAMsorted = os.path.abspath(BAMsorted)
    os.chdir("..")
    return BAMsorted


def _run_pileup(global_config, bamfile):
    """
    Perform samtools pileup on a .bam-file
    """


    if "samtools" in global_config["Tools"]:
        samtools = global_config["Tools"]["samtools"]["bin"]
    elif not common.which("samtools"):
        sys.exit("error while trying to run samtools: samtools not present in the path and not in global config, please make sure to install samtools properly")

    pileupfile = bamfile.replace('.bam', '_coverage.csv')
    pileup_cmd = "{} mpileup {} | awk '{print $2, $4}' > {}".format(samtools, bamfile, pileupfile)
    p1 = subprocess.Popen(pileup_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    p1.wait()
    if p1.returncode == 0:
        return pileupfile
    else:
        print "Could not perform mpileup"
        return 1


def plot_coverage(pileupfile, samplename=None):
    """
    Plot coverage from samtools pileup output
    """
    if not samplename:
        samplename = pileupfile.split('.')[0]

    df = pd.io.parsers.read_csv(pileupfile, sep=' ', names=['pos', 'cov'])
    pl = df.plot(x='pos', y='cov')
    avg_cov = sum(df['cov'])/len(df)        #Calculate average coverage
    pl2 = plt.plot(range(len(df)), [avg_cov]*len(df), 'r--', linewidth=2)
    plt.xlim([0,len(df)])
    plt.ylabel('Coverage')
    plt.xlabel('Position')
    plt.title(samplename)
    plt.savefig(samplename)     #Save plot as .png
    return 0







