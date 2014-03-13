import sys, os, yaml, glob
import subprocess
import argparse
from collections import Counter
import pandas as pd
import csv


def main(args):
    workingDir     = os.getcwd()
    assembly       = args.assembly
    minCtgLength   = args.min_contig
    
    assembly        = _build_new_reference(assembly, minCtgLength)
    if args.only_reference > 0:
        return
    ## now generate stats
    (contigsLengthDict, contigsSequencesDict) = _compute_assembly_stats(assembly, args.genome_size)
    if args.only_stats > 0:
        return

    #problems = find_problems_in_maps(args.opgen_report, contigsLengthDict, contigsSequencesDict)

    produce_consensus(args.opgen_report, contigsLengthDict, contigsSequencesDict, args.output)
    
    return



def produce_consensus(opgen_report, contigsLengthDict, contigsSequencesDict, output):
    print opgen_report
    # starting to parse opgen report, since the csv file is unstructured I need to parse it
    placement        = "{}.placement".format(opgen_report)
    stats            = "{}.stats".format(opgen_report)
    gaps_overlaps    = "{}.gaps".format(opgen_report)
    unplaced_contigs = "{}.unplaced".format(opgen_report)
    with open(opgen_report, 'rb') as unstructured_file:
        placement_file          = open(placement, 'w')
        stats_file              = open(stats, 'w')
        gaps_overlaps_file      = open(gaps_overlaps, 'w')
        unplaced_contigs_file   = open(unplaced_contigs, 'w')
        current_writing_file    = placement_file
        for row in unstructured_file:
            if "N50" in row:
                current_writing_file = stats_file
            if "Gaps" in row:
                current_writing_file = gaps_overlaps_file
            if "Unplaced" in row:
                current_writing_file = unplaced_contigs_file
            current_writing_file.write(row)
        current_writing_file.close()
        placement_file.close()
        stats_file.close()
        gaps_overlaps_file.close()

   # I need a list with one element for each base in my map. THe problem is that the OpGen report does not tell the length of the maps, therefore I need to find the maximum limit between the placments part and the gapped part
    Maps            = {}
    Maps_length     = {}
    ContigsToMaps   = {}
    MapsToContigs   = {}
    with open(placement, 'rb') as csvfile:
        opgenCSV = csv.reader(csvfile, delimiter='\t')
        #Chromosome, Start, End, Contig, Start, End, Orientation
        header = opgenCSV.next()
        for row in opgenCSV:
            if len(row) > 0:
                #print ', '.join(row)
                (Map, Map_Start, Map_End, Contig, Ctg_Start, Ctg_End, Orientation) = (row[0], int(row[1]), int(row[2]), row[3], int(row[4]), int(row[5]), row[6])
                Contig = Contig.split(" ")[0]
                if not ContigsToMaps.has_key(Contig):
                    ContigsToMaps[Contig] = [[Map, Map_Start, Map_End, Ctg_Start, Ctg_End, Orientation]]
                else:
                    ContigsToMaps[Contig].append([Map, Map_Start, Map_End, Ctg_Start, Ctg_End, Orientation])
                if not Maps_length.has_key(Map):
                    Maps_length[Map]   = Map_End
                    MapsToContigs[Map] = [[Map_Start, Map_End, Contig, Ctg_Start, Ctg_End, Orientation]]
                else:
                    if Map_End > Maps_length[Map]:
                        Maps_length[Map] = Map_End
                        MapsToContigs[Map].append([Map_Start, Map_End, Contig, Ctg_Start, Ctg_End, Orientation])
                    else:
                        print "{} {}".format(Map_End, Maps_length[Map])

    with open(gaps_overlaps, 'rb') as csvfile:
        opgenCSV = csv.reader(csvfile, delimiter='\t')
        #Map     Type    Map_Start   Map_End     Length
        header = opgenCSV.next()
        header = opgenCSV.next()
        for row in opgenCSV:
            if len(row) > 0:
                (Map, Type, Map_Start, Map_End, Length) = (row[0], row[1], int(row[2]), int(row[3]), int(row[4]))
                if not Maps_length.has_key(Map):
                    Maps_length[Map] = Map_End
                else:
                    if Map_End > Maps_length[Map]:
                        Maps_length[Map] = Map_End

    # now Maps_length contains the lenght of the Maps... thank you so much OpGen!!!!!
    multipleHittedPositions = 0
    rescuedBases            = 0
    conflictingBases        = 0
    for Map in Maps_length:
        Maps[Map] = ["n"] * Maps_length[Map]
        print "now working with Map {}".format(Map)
        if MapsToContigs[Map] is not None: # if this map has some contigs aligning
            for hit in MapsToContigs[Map]:
                start_on_Map = hit[0]
                end_on_Map   = hit[1]
                Ctg          = hit[2]
                start_on_Ctg = hit[3]
                end_on_Ctg   = hit[4]
                orientation  = hit[5]
                #extract the seuqence
                sequence     = contigsSequencesDict[Ctg][start_on_Ctg:end_on_Ctg]
                if orientation == '-1':
                    sequence = revcom(sequence)
                letters = list(sequence)
                index = start_on_Map
                
                for base in letters:
                    if len(Maps[Map]) <= index:
                        Maps[Map].append("n");
                    if Maps[Map][index] is not "n":
                        multipleHittedPositions += 1
                        if Maps[Map][index] is "N" and base is not "N":
                            Maps[Map][index] = base
                            rescuedBases += 1
                        elif Maps[Map][index] is not "N" and base is not "N":
                            conflictingBases += 1
                    else:
                        Maps[Map][index] = base
                    index += 1
                #print hit
    print "position in th emap covered more than once {}".format(multipleHittedPositions)
    print "rescued bases {}".format(rescuedBases)
    print "conflicting bases {}".format(conflictingBases)
    print_contigs(Maps, output)


def revcom(s):
    return complement(s[::-1])

def complement(s):
    letters      = list(s)
    complemented = []
    for base in letters:
        if base is "A" or base is "a":
            complemented.append("T")
        elif base is 'C' or base is 'c':
            complemented.append("G")
        elif base is 'G' or base is 'g':
            complemented.append("C")
        elif base is 'T' or base is 't':
            complemented.append("A")
        elif base is 'N' or base is 'n':
            complemented.append("N")
        else:
            complemented.append("N")

    return ''.join(complemented)


def print_contigs(Maps, output):
    with open("{}.fasta".format(output), "w") as final_assembly:
        for Map, sequence in Maps.iteritems():
            final_assembly.write(">{}\n".format(Map))
            final_assembly.write("{}\n".format("".join(sequence)))



def find_problems_in_maps(opgen_report, contigsLengthDict, contigsSequencesDict):
    ###TODO: this function is still to be defined.
    print opgen_report
    # startint to parse opgen report, since the csv file is unstructured I need to parse it
    placement        = "{}.placement".format(opgen_report)
    stats            = "{}.stats".format(opgen_report)
    gaps_overlaps    = "{}.gaps".format(opgen_report)
    unplaced_contigs = "{}.unplaced".format(opgen_report)
    with open(opgen_report, 'rb') as unstructured_file:
        placement_file          = open(placement, 'w')
        stats_file              = open(stats, 'w')
        gaps_overlaps_file      = open(gaps_overlaps, 'w')
        unplaced_contigs_file   = open(unplaced_contigs, 'w')
        current_writing_file    = placement_file
        for row in unstructured_file:
            if "N50" in row:
                current_writing_file = stats_file
            if "Gaps" in row:
                current_writing_file = gaps_overlaps_file
            if "Unplaced" in row:
                current_writing_file = unplaced_contigs_file
            current_writing_file.write(row)
        current_writing_file.close()
        placement_file.close()
        stats_file.close()
        gaps_overlaps_file.close()

    # I need a list with one element for each base in my map. THe problem is that the OpGen report does not tell the length of the maps, therefore I need to find the maximum limit between the placments part and the gapped part
    Maps            = {}
    Maps_length     = {}
    ContigsToMaps   = {}
    with open(placement, 'rb') as csvfile:
        opgenCSV = csv.reader(csvfile, delimiter='\t')
        #Chromosome, Start, End, Contig, Start, End, Orientation
        header = opgenCSV.next()
        for row in opgenCSV:
            if len(row) > 0:
                #print ', '.join(row)
                (Map, Map_Start, Map_End, Contig, Ctg_Start, Ctg_End, Orientation) = (row[0], int(row[1]), int(row[2]), row[3], int(row[4]), int(row[5]), row[6])
                Contig = Contig.split(" ")[0]
                if not ContigsToMaps.has_key(Contig):
                    ContigsToMaps[Contig] = [[Map, Map_Start, Map_End, Ctg_Start, Ctg_End, Orientation]]
                else:
                    ContigsToMaps[Contig].append([Map, Map_Start, Map_End, Ctg_Start, Ctg_End, Orientation])
                if not Maps_length.has_key(Map):
                    Maps_length[Map] = Map_End
                else:
                    if Map_End > Maps_length[Map]:
                        Maps_length[Map] = Map_End
                    else:
                        print Map_End + " " + Maps_length[Map]

    with open(gaps_overlaps, 'rb') as csvfile:
        opgenCSV = csv.reader(csvfile, delimiter='\t')
        #Map     Type    Map_Start   Map_End     Length
        header = opgenCSV.next()
        header = opgenCSV.next()
        for row in opgenCSV:
            if len(row) > 0:
                (Map, Type, Map_Start, Map_End, Length) = (row[0], row[1], int(row[2]), int(row[3]), int(row[4]))
                if not Maps_length.has_key(Map):
                    Maps_length[Map] = Map_End
                else:
                    if Map_End > Maps_length[Map]:
                        Maps_length[Map] = Map_End
    # now Maps_length contains the lenght of the Maps... thank you so much OpGen!!!!!

    # print contigs that align in more than one map
    ##### now decide if these contigs are suspicious (need to be broken or they are ok)
    ##### if a contig maps in two different places in the same way, then this is a real repeat. If there are two clear alignments then I must suggest to breack the contig
    contigsInMultipleMaps = 0
    contigsWithErrors     = 0
    for Contig in ContigsToMaps:
        length = len(ContigsToMaps[Contig])
        if length > 1:
            contigsInMultipleMaps += 1
            contigsCoordinates     = []
            print Contig + " --> "
            for entry in ContigsToMaps[Contig]:
                contigsCoordinates.append([entry[3],entry[4], (entry[3] - entry[3])/contigsLengthDict[Contig] ]) # save only start, end, and length of alignment
                print entry
            ### check if this contig is an erroneus contig or a real repeat
    print "number of contigs aligning in multiple maps is {}".format(contigsInMultipleMaps)






    for Map in Maps_length:
        Maps[Map] = [0] * Maps_length[Map]






def _compute_assembly_stats(assembly, genomeSize):
    stats_file_name = ".".join([assembly, "statistics", "txt"])
    contigsLengthDict   = {}
    contigsSequencesDict= {}
    contigsLength       = []
    totalLength         = 0
    numContigs          = 0
    maxContigLength     = 0
    N50                 = 0
    N80                 = 0
    numNs               = 0
    
    with open(assembly, "r") as ref_fd:
        fasta_header     = ref_fd.readline()
        header           = (fasta_header.split(">")[1]).rstrip()
        sequence         = ""
        for line in ref_fd:
            if line.startswith(">"):
                contigsLength.append(len(sequence))
                contigsLengthDict[header]    = len(sequence)
                contigsSequencesDict[header] = sequence
                totalLength              += len(sequence)
                numContigs               += 1
                counter                   = Counter(sequence)
                numNs                    += (counter['n'] + counter['N'])
                sequence                  = ""
                header                    = (line.split(">")[1]).rstrip()
            else:
                sequence+=line.rstrip()
        contigsLength.append(len(sequence))
        contigsLengthDict[header] = len(sequence)
        contigsSequencesDict[header] = sequence
        totalLength              += len(sequence)
        numContigs               += 1
        counter                   = Counter(sequence)
        numNs                    += (counter['n'] + counter['N'])

    if os.path.exists(stats_file_name):
        print "assembly stast file {} already created".format(stats_file_name)
        return (contigsLengthDict, contigsSequencesDict)


    percentageNs = float(numNs)/totalLength
    contigsLength.sort()
    contigsLength.reverse()

    
    teoN50 = genomeSize * 0.5
    teoN80 = genomeSize * 0.8
    testSum = 0
    N50 = 0
    N80 = 0
    maxContigLength   = contigsLength[0]
    for con in contigsLength:
        testSum += con
        if teoN50 < testSum:
            if N50 == 0:
                N50 = con
        if teoN80 < testSum:
            N80 = con
            break

    with open(stats_file_name, "w") as stats_file_name_fd:
        stats_file_name_fd.write("#scaffolds : {}\n".format(numContigs))
        stats_file_name_fd.write("totalLength : {}\n".format(totalLength))
        stats_file_name_fd.write("maxContigLength : {}\n".format(maxContigLength))
        stats_file_name_fd.write("N50 : {}\n".format(N50))
        stats_file_name_fd.write("N80 : {}\n".format(N80))
        stats_file_name_fd.write("percentageNs : {}\n".format(percentageNs))

    return (contigsLengthDict, contigsSequencesDict)



def _build_new_reference(assembly, minCtgLength):
    new_assembly_name    = os.path.basename(assembly)
    (basename, type, extention) = new_assembly_name.split(os.extsep)
    new_assembly_name    = basename + ".{}.{}.{}".format(type, minCtgLength, extention)
    new_assembly_name =  os.path.abspath(os.path.basename(new_assembly_name))

    if os.path.exists(new_assembly_name):
        print "assembly {} already created".format(new_assembly_name)
        return new_assembly_name
    
    contigID = 1
    with open(new_assembly_name, "w") as new_ref_fd:
        with open(assembly, "r") as ref_fd:
            fasta_header    = ref_fd.readline()
            sequence        = ""
            for line in ref_fd:
                line = line
                if line.startswith(">"):
                    if len(sequence) >= minCtgLength:
                        new_ref_fd.write(">contig{}\n".format(contigID))
                        new_ref_fd.write(sequence)
                    sequence        = ""
                    contigID        += 1
                else:
                    sequence+=line
            if len(sequence) >= minCtgLength:
                new_ref_fd.write(">contig{}".format(contigID))
                new_ref_fd.write(sequence)
    return new_assembly_name



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--assembly', help="assembly sequence", type=str)
    parser.add_argument('--genome-size', help="expected genome size", default=0, type=int)
    parser.add_argument('--min-contig', help="minimum contig length", default=2000, type=int)
    parser.add_argument('--only-reference', help="generates only the new refernce (i.e., only scaffolds longer than min-contig)", default=0, type=int)
    parser.add_argument('--only-stats', help="generates only assembly stats", default=0, type=int)

    parser.add_argument('--opgen-report', help="op gen placment report", type=str)
    parser.add_argument('--find-problems-in-maps', help="generates only the new refernce (i.e., only scaffolds longer than min-contig)", default=0, type=int)
    parser.add_argument('--produce-consensus', help="generates only the new refernce (i.e., only scaffolds longer than min-contig)", default=0, type=int)

    parser.add_argument('--output', help="aoutput header name", default="opgen_scaffolded_assembly", type=str)
    args = parser.parse_args()
    if args.genome_size == 0:
        print "genome size must be specified"
        sys.exit("error")

    main(args)



