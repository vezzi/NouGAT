import sys, os, yaml, glob
import subprocess
import string
import sys





def run(global_config, sample_config):
    tool = sample_config["operation"]["assemble"]["tool"]
    if tool == "masurca":
        _run_masurca(global_config, sample_config)
    elif tool == "soapdenovo":
        _run_soapdenovo(global_config, sample_config)
    elif tool == "abyss":
        _run_abyss(global_config, sample_config)
    elif tool == "spades":
        _run_spades(global_config, sample_config)
    elif tool == "abyss_mergePairs":
        _run_abyss_mergePairs(global_config, sample_config)
    elif tool == "cabog":
        _run_cabog(global_config, sample_config)
    else:
        print "tool {} is not yet supported".format(tool)

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


def _run_soapdenovo(global_config, sample_config):
    print "running SOAPdenovo ..."
    outputName = sample_config["output"]
    mainDir = os.getcwd()
    soapFolder = os.path.join(os.getcwd(), "soapdenovo")
    if not os.path.exists(soapFolder):
        os.makedirs(soapFolder)
    else:
        print "done (soapdenovo folder already present, assumed already run)"
        return
    os.chdir(soapFolder)
    ##create DATA directory
    sorted_libraries_by_insert = sorted(sample_config["libraries"].iteritems(), key=lambda (k,v): v["insert"])
    sorted_libraries_by_insert = prepare_folder_structure(sorted_libraries_by_insert)

    programBIN      = global_config["assemble"]["soapdenovo"]["bin"] # in masurca case there is no exectuable as a make file must be created
    program_options =global_config["assemble"]["soapdenovo"]["options"]
    kmer            = sample_config["kmer"]

    soap_config_file = open("configuration.txt", "w")
    soap_config_file.write("max_rd_len=250\n") #TODO make this a parameter in the options
    rank = 1
    for library, libraryInfo in sorted_libraries_by_insert:
        soap_config_file.write("[LIB]\n")
        read1       =libraryInfo["pair1"]
        read2       =libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert      = libraryInfo["insert"]
        std         = libraryInfo["std"]
        soap_config_file.write("avg_ins={}\n".format(insert))
        soap_config_file.write("rank={}\n".format(rank))
        rank += 1
        soap_config_file.write("map_len=30\n")
        if orientation=="innie" or orientation=="none":
            soap_config_file.write("asm_flags=3\n")
            soap_config_file.write("pair_num_cutoff=3\n")
            soap_config_file.write("reverse_seq=0\n")
            if read2 is None:
                soap_config_file.write("q={}\n".format(read1))
            else:
                soap_config_file.write("q1={}\n".format(read1))
                soap_config_file.write("q2={}\n".format(read2))
        elif orientation=="outtie":
            soap_config_file.write("asm_flags=2\n")
            soap_config_file.write("pair_num_cutoff=5\n")
            soap_config_file.write("reverse_seq=1\n")
            soap_config_file.write("q1={}\n".format(read1))
            soap_config_file.write("q2={}\n".format(read2))

    soap_config_file.close()
    soap_stdOut = open("soap.stdOut", "w")
    soap_stdErr = open("soap.stdErr", "w")
    os.makedirs(os.path.join(soapFolder, "runSOAP"))
    os.chdir("runSOAP")
    command = [programBIN , "all", "-s", "../configuration.txt", "-K", "{}".format(kmer),  "-o", "soapAssembly", "-p" , "8" ] # TODO: make everyhitn more parameter dependent
    returnValue = subprocess.call(command, stdout=soap_stdOut, stderr=soap_stdErr)
    os.chdir("..")
    if returnValue == 0:
        if(os.path.exists(os.path.join("runSOAP","soapAssembly.scafSeq"))):
            subprocess.call(["mv", os.path.join("runSOAP","soapAssembly.scafSeq"), "{}.scf.fasta".format(outputName) ])
            subprocess.call(["mv", os.path.join("runSOAP","soapAssembly.contig"), "{}.ctg.fasta".format(outputName) ])
            subprocess.call(["rm", "-r", "runSOAP"])
        else:
            print "something wrong with SOAPdenovo -> no contig file generated"
    else:
        print "SOAPdenovo terminated with an error. Please check running folder for more informations"
    os.chdir("..")
    print "SOAPdenovo exectued"
    return


def _run_masurca(global_config, sample_config):
    print "running Masurca ..."
    mainDir = os.getcwd()
    masurcaFolder = os.path.join(os.getcwd(), "masurca")
    if not os.path.exists(masurcaFolder):
        os.makedirs(masurcaFolder)
    else:
        print "done (masurca folder already present, assumed already run)"
        return
    os.chdir(masurcaFolder)
    ##create DATA directory
    sorted_libraries_by_insert = sorted(sample_config["libraries"].iteritems(), key=lambda (k,v): v["insert"])
    sorted_libraries_by_insert = prepare_folder_structure(sorted_libraries_by_insert)

    programPATH = global_config["assemble"]["masurca"]["bin"] # in masurca case there is no exectuable as a make file must be created
    program_options=global_config["assemble"]["masurca"]["options"]

    masurca_config_file = open("configuration.txt", "w")
    masurca_config_file.write("PATHS\n")
    masurca_config_file.write("JELLYFISH_PATH=" + os.path.join("{}".format(programPATH), "bin" ) + "\n" )
    masurca_config_file.write("SR_PATH=" +  os.path.join("{}".format(programPATH), "bin" ) + "\n" )
    masurca_config_file.write("CA_PATH=" +  os.path.join("{}".format(programPATH), "CA/Linux-amd64/bin") + "\n")
    masurca_config_file.write("END\n")
    
    masurca_config_file.write("\n")
    
    masurca_config_file.write("DATA\n")
    allTheLetters = string.lowercase
    libraryPE    = "p"
    libraryPEnum = 0
    libraryMP    = "m"
    libraryMPnum = 0
#TODO: single ended reads
    for library, libraryInfo in sorted_libraries_by_insert:
        read1=libraryInfo["pair1"]
        read2=libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert      = libraryInfo["insert"]
        std         = libraryInfo["std"]
        if orientation=="innie":
            if read2 is not None:
                configurationLine = "PE = {}{} {} {} {} {}".format(libraryPE, allTheLetters[libraryPEnum], insert, std, read1, read2)
                masurca_config_file.write("{}\n".format(configurationLine))
                libraryPEnum+=1 #TODO: check when more than 21 PE libraries ae specified
        elif orientation=="outtie":
            configurationLine = "JUMP = {}{} {} {} {} {}".format(libraryMP, allTheLetters[libraryMPnum], insert, std, read1, read2)
            masurca_config_file.write("{}\n".format(configurationLine))
            libraryMPnum += 1  #TODO: check when more than 21 PE libraries ae specified
    masurca_config_file.write("END\n")
    
    masurca_config_file.write("\n")
    
    masurca_config_file.write("PARAMETERS\n")
    #this is k-mer size for deBruijn graph values between 25 and 101 are supported, auto will compute the optimal size based on the read data and GC content
    masurca_config_file.write("GRAPH_KMER_SIZE=auto\n")
    #set this to 1 for Illumina-only assemblies and to 0 if you have 2x or more long (Sanger, 454) reads
    masurca_config_file.write("USE_LINKING_MATES=1\n")
    #this parameter is useful if you have too many jumping library mates. Typically set it to 60 for bacteria and something large (300) for mammals
    masurca_config_file.write("LIMIT_JUMP_COVERAGE = 60\n")
    #these are the additional parameters to Celera Assembler.  do not worry about performance, number or processors or batch sizes -- these are computed automatically. for mammals do not set cgwErrorRate above 0.15!!!
    masurca_config_file.write("CA_PARAMETERS = ovlMerSize=30 cgwErrorRate=0.25 ovlMemory=4GB\n")
    #minimum count k-mers used in error correction 1 means all k-mers are used.  one can increase to 2 if coverage >100
    masurca_config_file.write("KMER_COUNT_THRESHOLD = 1\n")
    #auto-detected number of cpus to use
    masurca_config_file.write("NUM_THREADS= 8\n")
    #this is mandatory jellyfish hash size
    masurca_config_file.write("JF_SIZE=1000000000\n")
    #this specifies if we do (1) or do not (0) want to trim long runs of homopolymers (e.g. GGGGGGGG) from 3' read ends, use it for high GC genomes
    masurca_config_file.write("DO_HOMOPOLYMER_TRIM=0\n")
    masurca_config_file.write("END\n")
    masurca_config_file.write("\n")

    masurca_config_file.close()

    masurca_stdOut = open("masurca.stdOut", "w")
    masurca_stdErr = open("masurca.stdErr", "w")
    command = [os.path.join(programPATH,"bin/runSRCA.pl") , "configuration.txt"]
    subprocess.call(command, stdout=masurca_stdOut, stderr=masurca_stdErr)

    if not os.path.exists("assemble.sh"):
        print "MaSuRCA: assemble.sh not created. Unknown failure"
        return
    command = ["./assemble.sh"]
    returnValue = subprocess.call(command, stdout=masurca_stdOut, stderr=masurca_stdErr)
    if returnValue == 0:
        if os.path.exists(os.path.join("runABySS","{}-contigs.fa".format(outputName))):
            subprocess.call(["cp", os.path.join("runABySS","{}-contigs.fa".format(outputName)), "{}.ctg.fasta".format(outputName) ])
            subprocess.call(["cp", os.path.join("runABySS","{}-scaffolds.fa".format(outputName)), "{}.scaf.fasta".format(outputName) ])
            subprocess.call(["rm", "-r", "runABySS"])
        else:
            print "something wrong with MaSuRCA -> no contig file generated"
    else:
        print "MaSuRCA terminated with an error. Please check running folder for more informations"
    os.chdir("..")
    print "MaSuRCA executed"
    return





def _run_abyss(global_config, sample_config):
    print "running abyss ..."
    ##TO DO ::: avoid to load module
    ##subprocess.call(["module","load","abyss/1.3.5"])
    outputName = sample_config["output"]
    mainDir = os.getcwd()
    abyssFolder = os.path.join(os.getcwd(), "abyss")
    if not os.path.exists(abyssFolder):
        os.makedirs(abyssFolder)
    else:
        print "done (abyss folder already present, assumed already run)"
        #return
   
    os.chdir(abyssFolder)
    ##create DATA directory
    sorted_libraries_by_insert = sorted(sample_config["libraries"].iteritems(), key=lambda (k,v): v["insert"])
    sorted_libraries_by_insert = prepare_folder_structure(sorted_libraries_by_insert)

    
    abyss_stdOut = open("abyss.stdOut", "a")
    abyss_stdErr = open("abyss.stdErr", "a")
    program=global_config["assemble"]["abyss"]["bin"]
    program_options=global_config["assemble"]["abyss"]["options"]
    
    command = ""
    command += "{} ".format(program)
    for option in program_options:
            command += "{} ".format(option)
    
    libraries = {}
    for library, libraryInfo in sorted_libraries_by_insert:
        read1       = libraryInfo["pair1"]
        read2       = libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert      = libraryInfo["insert"]
        std         = libraryInfo["std"]
        if orientation=="innie" or orientation=="none":
            if read2 is None:
                if "se" not in libraries: # check if this is the first time I insert a se file
                    libraries["se"] = "se=\'"
                libraries["se"] = libraries["se"] + read1
            else:
                if not "lib" in libraries:
                    libraries["lib"] = {}
                libName = "pe{}".format(insert)
                if not libName in libraries["lib"]:
                    libraries["lib"][libName] = ""
                libraries["lib"][libName] += libraries["lib"][libName] + "{} {} ".format(read1, read2)
        else:
            if not "mp" in libraries:
                libraries["mp"] = {}
            libName = "mp{}".format(insert)
            if not libName in libraries["mp"]:
                libraries["mp"][libName] = ""
            libraries["mp"][libName] += libraries["mp"][libName] + "{} {} ".format(read1, read2)
    #now create the command
    command += "name={} ".format(outputName)
    librariesSE       = ""
    librariesPE       = ""
    librariesMP       = ""
    if "se" in libraries:
        libraries["se"] = libraries["se"] + "\'"
        librariesSE = libraries["se"]
    if "lib" in libraries:
        lib="lib=\'"
        for libPE, libPEreads in libraries["lib"].iteritems():
            lib = lib + "{} ".format(libPE)
            librariesPE += " {}=\'{}\' ".format(libPE,libPEreads)
        lib=lib + "\' "
        command += "{} ".format(lib)
    if "mp" in libraries:
        mp="mp=\'"
        for libMP, libMPreads in libraries["mp"].iteritems():
            mp = mp + "{} ".format(libMP)
            librariesMP += " {}=\'{}\' ".format(libMP,libMPreads)
        mp=mp + "\' "
        command += "{} ".format(mp)

    command += "{} ".format(librariesSE)
    command += "{} ".format(librariesPE)
    command += "{} ".format(librariesMP)

    os.makedirs(os.path.join(abyssFolder, "runABySS"))
    os.chdir("runABySS")
    returnValue = subprocess.call(command, stdout=abyss_stdOut, stderr=abyss_stdErr, shell=True)
    os.chdir("..")
    if returnValue == 0:
        if os.path.exists(os.path.join("runABySS","{}-contigs.fa".format(outputName))):
            subprocess.call(["cp", os.path.join("runABySS","{}-contigs.fa".format(outputName)), "{}.ctg.fasta".format(outputName) ])
            subprocess.call(["cp", os.path.join("runABySS","{}-scaffolds.fa".format(outputName)), "{}.scaf.fasta".format(outputName) ])
            subprocess.call(["rm", "-r", "runABySS"])
        else:
            print "something wrong with ABySS -> no contig file generated"
    else:
        print "ABySS terminated with an error. Please check running folder for more informations"
    os.chdir("..")
    print "abyss executed"
    return





def _run_abyss_mergePairs(global_config, sample_config):
    print "running abyss-mergepairs ..."
    mainDir = os.getcwd()
    abyssFolder = os.path.join(os.getcwd(), "abyss_mergePairs")
    if not os.path.exists(abyssFolder):
        os.makedirs(abyssFolder)
    else:
        print "done (abyss_mergePairs folder already present, assumed already run)"
        #return
   
    os.chdir(abyssFolder)
##create DATA directory
    sorted_libraries_by_insert = sorted(sample_config["libraries"].iteritems(), key=lambda (k,v): v["insert"])
    sorted_libraries_by_insert = prepare_folder_structure(sorted_libraries_by_insert)

    program=global_config["assemble"]["abyss_mergePairs"]["bin"]
    program_options=global_config["assemble"]["abyss_mergePairs"]["options"]
    
    command = []
    command.append(program)
    for option in program_options:
            command.append(option)
    
    libraries = {}
    for library, libraryInfo in sorted_libraries_by_insert:
        read1       = libraryInfo["pair1"]
        read2       = libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert      = libraryInfo["insert"]
        std         = libraryInfo["std"]
        outputNameArray  = read1.split('/')[-1].split('_')
        outputName = "{}_{}".format(outputNameArray[0], outputNameArray[1])
        
        if orientation=="innie":
            if read2 is not None:
                currentCommand = command;
                currentCommand.append('-o')
                currentCommand.append(outputName)
                currentCommand.append(read1)
                currentCommand.append(read2)
                print currentCommand
                abyss_stdOut = open("mergePairs_{}.stdOut".format(outputName), "a")
                abyss_stdErr = open("mergePairs_{}.stdErr".format(outputName), "a")
                
                subprocess.call(command, stdout=abyss_stdOut, stderr=abyss_stdErr)
                command_mv = ["mv", "mergePairs_{}.stdErr".format(outputName), "{}.txt".format(outputName)]
                subprocess.call(command_mv)

    os.chdir("..")
    
    return






def _run_spades(global_config, sample_config):
    print "running spades ..."
    outputName = sample_config["output"]
    mainDir = os.getcwd()
    assemblerFolder = os.path.join(os.getcwd(), "spades")
    if not os.path.exists(assemblerFolder):
        os.makedirs(assemblerFolder)
    else:
        print "done (spades folder already present, assumed already run)"
        #return
   
    os.chdir(assemblerFolder)
    ##create DATA directory
    sorted_libraries_by_insert = sorted(sample_config["libraries"].iteritems(), key=lambda (k,v): v["insert"])
    sorted_libraries_by_insert = prepare_folder_structure(sorted_libraries_by_insert)

    spades_stdOut = open("spades.stdOut", "a")
    spades_stdErr = open("spades.stdErr", "a")

    program=global_config["assemble"]["spades"]["bin"]
    program_options=global_config["assemble"]["spades"]["options"]
    
    command = ""
    command += "{} ".format(program)
    for option in program_options:
        command += "{} ".format(option)

    #creates the command on-the-fly
    peLibrary = 1
    mpLibrary = 1
    for library, libraryInfo in sorted_libraries_by_insert:
        read1       = libraryInfo["pair1"]
        read2       = libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert      = libraryInfo["insert"]
        std         = libraryInfo["std"]
        if orientation=="innie" or orientation=="none":
            if read2 is None:
                command += "--pe{}-s {} ".format(peLibrary, read1)
            else:
                command += "--pe{}-1 {} --pe{}-2 {} ".format(peLibrary, read1, peLibrary, read2)
            peLibrary += 1
        elif orientation=="outtie":
            command += "--mp{}-1 {} --mp{}-2 {} ".format(mpLibrary, read1, mpLibrary, read2)
            mpLibrary += 1
        else:
            print "orientation{} not supported.... why the program didnot failed earlier?".format(orientation)

    command += "-o {} ".format(outputName)

    returnValue = subprocess.call(command, stdout=spades_stdOut, stderr=spades_stdErr, shell=True)
    if returnValue == 0:
        if os.path.exists(os.path.join("outputName","contigs.fasta")):
            subprocess.call(["cp", os.path.join("outputName","contigs.fasta"),  "{}.ctg.fasta".format(outputName)])
            subprocess.call(["cp", os.path.join("outputName","scaffoldss.fasta"),  "{}.scf.fasta".format(outputName)])
            subprocess.call(["rm", "-r", outputName])
        else:
            print "something wrong with SPADES -> no contig file generated"
    else:
        print "SPADES terminated with an error. Please check running folder for more informations"

    os.chdir("..")
    return



def _run_cabog(global_config, sample_config):
    print "running CABOG ..."
    outputName = sample_config["output"]
    mainDir = os.getcwd()
    assemblerFolder = os.path.join(os.getcwd(), "cabog")
    if not os.path.exists(assemblerFolder):
        os.makedirs(assemblerFolder)
    else:
        print "done (cabog folder already present, assumed already run)"
        return
    os.chdir(assemblerFolder)
    ##create DATA directory
    sorted_libraries_by_insert = sorted(sample_config["libraries"].iteritems(), key=lambda (k,v): v["insert"])
    sorted_libraries_by_insert = prepare_folder_structure(sorted_libraries_by_insert)


    programBIN = global_config["assemble"]["cabog"]["bin"] # in masurca case there is no exectuable as a make file must be created
    program_options=global_config["assemble"]["cabog"]["options"]

    sys.path.insert(0, programBIN)
    libraries = 1
    cabog_stdOut = open("cabog.stdOut", "w")
    cabog_stdErr = open("cabog.stdErr", "w")
    for library, libraryInfo in sorted_libraries_by_insert:
        command_fastqToCA = os.path.join(programBIN, "fastqToCA")
        read1=libraryInfo["pair1"]
        read2=libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert      = libraryInfo["insert"]
        std         = libraryInfo["std"]
        command_fastqToCA  += " -libraryname "
        command_fastqToCA  += " {}_{}".format(outputName, libraries)
        command_fastqToCA  += " -insertsize "
        command_fastqToCA  += " {} {} ".format(insert,std)
        command_fastqToCA  += " -technology "
        command_fastqToCA  += " illumina "
        command_fastqToCA  += " -type "
        command_fastqToCA  += " illumina "
        
        if orientation=="innie" or orientation=="none" :
            command_fastqToCA  += " -innie "
            if read2 is None:
                command_fastqToCA  += " -reads "
                command_fastqToCA  += " {} ".format(read1)
            else:
                command_fastqToCA  += " -mates "
                command_fastqToCA  += " {},{} ".format(read1, read2)
        elif orientation=="outtie":
            command_fastqToCA  += " -outtie "
            command_fastqToCA  += " -mates "
            command_fastqToCA  += " {},{} ".format(read1, read2)
        command_fastqToCA  += " > "
        command_fastqToCA  += " {}_{}.frg ".format(outputName, libraries)
        subprocess.call(command_fastqToCA, stderr=cabog_stdErr, shell=True)
        libraries += 1
    
    command_runCA = os.path.join(programBIN, "runCA")
    command_runCA += "  -d runCABOGfolder -p {} *frg".format(outputName)
    returnValue = subprocess.call(command_runCA, stdout=cabog_stdOut, stderr=cabog_stdErr, shell=True)
    if returnValue == 0:
        #assembly succed, remove files and save assembly
        if os.path.exists(os.path.join("runCABOGfolder","9-terminator", "{}.ctg.fasta".format(outputName))):
            subprocess.call(["mv", os.path.join("runCABOGfolder","9-terminator", "{}.ctg.fasta".format(outputName)), "{}.ctg.fasta".format(outputName)])
            subprocess.call(["mv", os.path.join("runCABOGfolder","9-terminator", "{}.scf.fasta".format(outputName)), "{}.scf.fasta".format(outputName)])
            subprocess.call(["rm", "-r", "runCABOGfolder"])
        else:
            print "something wrong with CABOG -> no contig file generated"
    else:
        print "CABOG terminated with an error. Please check running folder for more informations"

    os.chdir("..")
    print "CABOG exectued"
    return





