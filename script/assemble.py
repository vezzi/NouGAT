import sys, os, yaml, glob
import subprocess
import string





def run(global_config, sample_config, sorted_libraries_by_insert):
    tool = sample_config["operation"]["assemble"]["tool"]
    if tool == "masurca":
        _run_masurca(global_config, sample_config, sorted_libraries_by_insert)
    elif tool == "soapdenovo":
        _run_soapdenovo(global_config, sample_config, sorted_libraries_by_insert)


def _run_soapdenovo(global_config, sample_config, sorted_libraries_by_insert):
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
    programBIN = global_config["assemble"]["soapdenovo"]["bin"] # in masurca case there is no exectuable as a make file must be created
    program_options=global_config["assemble"]["soapdenovo"]["options"]

    soap_config_file = open("configuration.txt", "w")
    soap_config_file.write("max_rd_len=150\n") #TODO make this a parameter in the options
    rank = 1
    for library, libraryInfo in sorted_libraries_by_insert:
        soap_config_file.write("[LIB]\n")
        read1=libraryInfo["pair1"]
        read2=libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert      = libraryInfo["insert"]
        std         = libraryInfo["std"]
        soap_config_file.write("avg_ins={}\n".format(insert))
        soap_config_file.write("rank={}\n".format(rank))
        rank += 1
        soap_config_file.write("map_len=30\n")
        
        if orientation=="innie":
            soap_config_file.write("asm_flags=3\n")
            soap_config_file.write("pair_num_cutoff=3\n")
            soap_config_file.write("reverse_seq=0\n")
            soap_config_file.write("q1={}\n".format(read1))
            soap_config_file.write("q2={}\n".format(read2))
        elif orientation=="outtie":
            soap_config_file.write("asm_flags=2\n")
            soap_config_file.write("pair_num_cutoff=5\n")
            soap_config_file.write("reverse_seq=1\n")
            soap_config_file.write("q1={}\n".format(read1))
            soap_config_file.write("q2={}\n".format(read2))

    soap_config_file.close()

    masurca_stdOut = open("soap.stdOut", "w")
    masurca_stdErr = open("soap.stdErr", "w")
    command = [programBIN , "all", "-s", "configuration.txt", "-K", "53",  "-o", "soapAssembly", "-p" , "8" ] # TODO: make everyhitn more parameter dependent
    subprocess.call(command, stdout=masurca_stdOut, stderr=masurca_stdErr)

    
    if(os.path.exists("soapAssembly.scafSeq")):
        os.rename("soapAssembly.scafSeq", "{}.fasta".format(outputName))
    else:
         print "something wrong with soapdenovo!!!"

    os.chdir("..")
    print "soap succesfully exectued"
    return


def _run_masurca(global_config, sample_config, sorted_libraries_by_insert):
    print "running Masurca ..."
    mainDir = os.getcwd()
    masurcaFolder = os.path.join(os.getcwd(), "masurca")
    if not os.path.exists(masurcaFolder):
        os.makedirs(masurcaFolder)
    else:
        print "done (masurca folder already present, assumed already run)"
        return
    os.chdir(masurcaFolder)
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

    for library, libraryInfo in sorted_libraries_by_insert:
        read1=libraryInfo["pair1"]
        read2=libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert      = libraryInfo["insert"]
        std         = libraryInfo["std"]
        if orientation=="innie":
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

    if(os.path.exists("assemble.sh")):
        command = ["assemble.sh"]
        subprocess.call(command, stdout=masurca_stdOut, stderr=masurca_stdErr)

    os.chdir("..")
    print "masurca succesfully exectued"
    return




def _run_abyss(global_config, sample_config, sorted_libraries_by_insert):
    print "running abyss ..."
    outputName = sample_config["output"]
    mainDir = os.getcwd()
    abyssFolder = os.path.join(os.getcwd(), "abyss")
    if not os.path.exists(abyssFolder):
        os.makedirs(abyssFolder)
    else:
        print "done (abyss folder already present, assumed already run)"
        #return
   
    os.chdir(abyssFolder)
    abyss_stdOut = open("abyss.stdOut", "a")
    abyss_stdErr = open("abyss.stdErr", "a")

    program=global_config["assemble"]["abyss"]["bin"]
    program_options=global_config["assemble"]["abyss"]["options"]
    command = [program]
    for option in program_options:
            command.append(option)
    libraries     = {}

    for library, libraryInfo in sorted_libraries_by_insert:
        read1       = libraryInfo["pair1"]
        read2       = libraryInfo["pair2"]
        orientation = libraryInfo["orientation"]
        insert      = libraryInfo["insert"]
        std         = libraryInfo["std"]
        if orientation=="innie":
            if read2 is None:
                if not libraries["se"]: # check if this is the first time I insert a se file
                    libraries["se"] = "se=\'"
                libraries["se"] = libraries["se"] + read1
            else:
                if not libraries["lib"]:
                    libraries["lib"] = {}
                libName = "pe{}".format(insert)
                if not libraries["lib"][libName]:
                    libraries["lib"][libName] = "{}=\'".format(libName)
                libraries["lib"][libName] += libraries["lib"][libName] + read1
                libraries["lib"][libName] += libraries["lib"][libName] + read1
        else:
                if not libraries["mp"]:
                    libraries["mp"] = {}
                libName = "mp{}".format(insert)
                if not libraries["mp"][libName]:
                    libraries["mp"][libName] = "{}=\'".format(libName)
                libraries["mp"][libName] += libraries["mp"][libName] + read1
                libraries["mp"][libName] += libraries["mp"][libName] + read1
        #now creat the command
        commad.append("k=54");
        command.append("name={}".format(outputName))

        if libraries["se"]:
            libraries["se"] = libraries["se"] + "\'"
            command.append(libraries["se"])
        if libraries["lib"]:
            lib="lib=\'"
            for libPE in libraries["lib"]:
                lib = lib + " {}".format(libPE)
            
        for library, libraryInfo in sorted_libraries_by_insert:
            read1       = libraryInfo["pair1"]
            read2       = libraryInfo["pair2"]
            orientation = libraryInfo["orientation"]
            insert      = libraryInfo["insert"]
            

        
        

            command.append(read1)
            
                command.append(read2)
        
    subprocess.call(command, stdout=abyss_stdOut, stderr=abyss_stdErr)
    os.chdir("..")
    
    return





