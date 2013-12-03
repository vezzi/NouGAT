import sys, os, yaml, glob
import subprocess
import common
import string
import sys
import gzip



def build_reference_bwa(reference):
    #build the reference if not available
    if not common.which("bwa"):
        sys.exit("error while trying to run  bwa index: bwa not present in the path, please make sure to install bwa properly")
    
    reference = os.path.abspath(reference)
    if not os.path.exists(reference):
        sys.exit("error, reference file {} does not exists".format(reference))
    
    if not os.path.exists("BWA_index"):
        os.makedirs("BWA_index")
    os.chdir("BWA_index")

    reference_local_copy = os.path.basename(reference)
    if not os.path.exists(reference_local_copy):
        print reference
        subprocess.call(["ln", "-s", reference, reference_local_copy])
    reference = reference_local_copy # now I have a soflinked copy

    if not os.path.exists("{}.bwt".format(reference)): # then create the index sequence
        bwa_stdOut = open("bwa_index.stdOut", "w")
        bwa_stdErr = open("bwa_index.stdErr", "w")
        returnValue = subprocess.call(["bwa", "index", reference], stdout=bwa_stdOut, stderr=bwa_stdErr)
        if  not returnValue == 0:
            sys.exit("error, while indexing reference file {} with bwa index".format(reference))
    reference = os.path.abspath(reference)
    os.chdir("..")
    if not os.path.exists("{}.bwt".format(reference)): #extra control to avoid problem with nuexpected return value
        sys.exit("bwa index failed")
    return reference


def align_bwa(read1, read2, reference):
    if not common.which("bwa"):
        sys.exit("error while trying to run  aling_bwa: bwa not present in the path, please make sure to install bwa properly")
    if not common.which("samtools"):
        sys.exit("error while trying to run aling_bwa: samtools not present in the path, please make sure to install samtools properly")
 
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

def align_bwa_mem(read1, read2, reference, threads):
    if not common.which("bwa"):
        sys.exit("error while trying to run  aling_bwa: bwa not present in the path, please make sure to install bwa properly")
    if not common.which("samtools"):
        sys.exit("error while trying to run aling_bwa: samtools not present in the path, please make sure to install samtools properly")
 
    libraryBase = os.path.basename(read1).split(".")[0]
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

    if not os.path.exists(SAMMapped):
        with open(SAMMapped, "w") as fh:
                returnValue = subprocess.call(["bwa", "mem", "-M", "-t", "{}".format(threads), reference, read1, read2], stdout=fh)
                if  not returnValue == 0:
                    sys.exit("error, while aligning reads with bwa mem")
    if not os.path.exists(BAMunsorted):
        with open(BAMunsorted, "w") as fh:
                returnValue = subprocess.call(["samtools", "view", "-S", "-F", "4", "-b", SAMMapped], stdout=fh)
                if  not returnValue == 0:
                    sys.exit("error, while while converting sam to bam")

    subprocess.call(["rm", SAMMapped])
    BAMsortedHeader = "{}".format(mappingBase)
    subprocess.call(["samtools", "sort",  BAMunsorted, BAMsortedHeader])
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
