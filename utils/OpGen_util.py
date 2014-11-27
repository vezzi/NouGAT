import sys, os, yaml, glob
import subprocess
import argparse
from collections import Counter, OrderedDict
from itertools import combinations
import pandas as pd
import csv
import re

reverse_complements = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "Y": "R",
    "R": "Y",
    "W": "W",
    "S": "S",
    "K": "M",
    "M": "K",
    "D": "H",
    "H": "D",
    "V": "B",
    "B": "V",
    "N": "N"
}
ambigous = [base for base in reverse_complements.keys() if base not in ["A", "T", "C", "G"]]

def main(args):
    workingDir = os.getcwd()
    assembly = args.assembly
    minCtgLength = args.min_contig
    if args.only_reference:
        assembly = _build_new_reference(assembly, minCtgLength)
        return

    (contigsLengthDict, contigsSequencesDict) = _compute_assembly_stats(assembly, args.genome_size)
    if args.only_stats:
        return

    placement, gaps_overlaps = parse_report(args.opgen_report) 
    if args.find_problems_in_maps:
        print("Finding problems by extracting overlaping contigs, please inspect ovl/*.fasta files!")
        find_problems_in_maps(placement, contigsLengthDict, contigsSequencesDict)
        return

    produce_consensus(placement, gaps_overlaps, contigsLengthDict, contigsSequencesDict, args.output, args.place_last)    
    return


def parse_report(opgen_report):
    # starting to parse opgen report, since the csv file is unstructured I need to parse it
    placement = "{}.placement".format(opgen_report)
    stats = "{}.stats".format(opgen_report)
    gaps_overlaps = "{}.gaps".format(opgen_report)
    unplaced_contigs = "{}.unplaced".format(opgen_report)
    with open(opgen_report, 'rb') as unstructured_file:
        placement_file = open(placement, 'w')
        stats_file = open(stats, 'w')
        gaps_overlaps_file = open(gaps_overlaps, 'w')
        unplaced_contigs_file = open(unplaced_contigs, 'w')
        current_writing_file = placement_file
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
    
    return (placement, gaps_overlaps)


def produce_consensus(placement, gaps_overlaps, contigsLengthDict, contigsSequencesDict, output, place_last):
    # I need a list with one element for each base in my map. THe problem is that
    # the OpGen report does not tell the length of the maps, therefore I need to 
    # find the maximum limit between the placments part and the gapped part
    Maps = {}
    Maps_length = {}
    ContigsToMaps = {}
    MapsToContigs = {}
    with open(placement, 'rb') as csvfile:
        opgenCSV = csv.reader(csvfile, delimiter='\t')
        #Chromosome, Start, End, Contig, Start, End, Orientation
        header = opgenCSV.next()
        for row in opgenCSV:
            if len(row) > 0:
                (Map, Map_Start, Map_End, Contig, Ctg_Start, Ctg_End, Orientation) = \
                        (row[0], int(row[1]), int(row[2]), row[3], int(row[4]), int(row[5]), row[6])
                Contig = Contig.split()[0]
                if not ContigsToMaps.has_key(Contig):
                    ContigsToMaps[Contig] = []
                
                if not Maps_length.has_key(Map):
                    Maps_length[Map] = Map_End
                    MapsToContigs[Map] = []
                elif Map_End > Maps_length[Map]:
                    Maps_length[Map] = Map_End
                
                MapsToContigs[Map].append([Map_Start, Map_End, Contig, Ctg_Start, Ctg_End, Orientation])
                ContigsToMaps[Contig].append([Map, Map_Start, Map_End, Ctg_Start, Ctg_End, Orientation])

    with open(gaps_overlaps, 'rb') as csvfile:
        opgenCSV = csv.reader(csvfile, delimiter='\t')
        #Map     Type    Map_Start   Map_End     Length
        header = opgenCSV.next()
        header = opgenCSV.next()
        for row in opgenCSV:
            if len(row) > 0:
                (Map, Type, Map_Start, Map_End, Length) = \
                        (row[0], row[1], int(row[2]), int(row[3]), int(row[4]))
                if not Maps_length.has_key(Map):
                    Maps_length[Map] = Map_End
                elif Map_End > Maps_length[Map]:
                    Maps_length[Map] = Map_End

    # now Maps_length contains the lenght of the Maps... thank you so much OpGen!!!!!
    multipleHittedPositions = 0
    rescuedBases = 0
    conflictingBases = 0

    for Map in Maps_length:
        Maps[Map] = ["n"] * Maps_length[Map] # Start with a blank canvas
        print "now working with Map {}".format(Map)
        if MapsToContigs[Map] is not None: # if this map has some contigs aligning
            reordered_tigs = [ctg for ctg in MapsToContigs[Map] if ctg[2].split("_")[0] not in place_last]
            last_tigs = [ctg for ctg in MapsToContigs[Map] if ctg[2].split("_")[0] in place_last]
            reordered_tigs.extend(last_tigs)

            for hit in reordered_tigs:
                start_on_Map = hit[0] - 1
                end_on_Map = hit[1] - 1
                Ctg = hit[2]
                start_on_Ctg = hit[3] - 1
                end_on_Ctg = hit[4] - 1
                orientation = hit[5]
                #extract the seuqence
                sequence = contigsSequencesDict[Ctg][start_on_Ctg:end_on_Ctg]
                if orientation == '-1':
                    sequence = revcom(contigsSequencesDict[Ctg][end_on_Ctg:start_on_Ctg:-1])
                letters = list(sequence)
                index = start_on_Map 

                for base in letters:
                    if len(Maps[Map]) <= index:
                        Maps[Map].append("n")
                    if Maps[Map][index] is not "n":
                        multipleHittedPositions += 1
                        if Maps[Map][index].title() in ambigous and base.title() not in ambigous:
                            rescuedBases += 1
                        elif Maps[Map][index].title() == base.title():
                            pass
                        elif Maps[Map][index].title() in ambigous and base.title() in ambigous:
                            pass
                        else:
                            conflictingBases += 1
                        Maps[Map][index] = base.title()
                    else:
                        Maps[Map][index] = base.title()
                    index += 1
    #print hit
    print "Positions in the map covered more than once {}".format(multipleHittedPositions)
    print "Rescued bases {}".format(rescuedBases)
    print "Conflicting bases {}".format(conflictingBases)
    print_contigs(Maps, output)


def revcom(s):
    return complement(s)


def complement(s):
    letters = list(s)
    complemented = []
    for base in letters:
        if base.title() in reverse_complements.keys():
            complemented.append(reverse_complements[base.title()])
        else:
            complemented.append("N")
    return ''.join(complemented)


def print_contigs(Maps, output):
    with open("{}.fasta".format(output), "w") as final_assembly:
        for Map, sequence in Maps.iteritems():
            final_assembly.write(">{}\n".format(output))
            final_assembly.write("{}\n".format("".join(sequence)))


def find_problems_in_maps(placement, contigsLengthDict, contigsSequencesDict):

    Map_contigs = OrderedDict()
    with open(placement, 'rb') as csvfile:
        opgenCSV = csv.reader(csvfile, delimiter='\t')
        #Chromosome, Start, End, Contig, Start, End, Orientation
        header = opgenCSV.next()
        for row in opgenCSV:
            if len(row) > 0:
                (Map, Map_Start, Map_End, Contig, Ctg_Start, Ctg_End, Orientation) = \
                    (row[0], int(row[1]), int(row[2]), row[3], int(row[4]), int(row[5]), row[6])
                Map_contigs[(Map_Start, Map_End)] = [Ctg_Start, Ctg_End, Contig.split()[0], Orientation]

    if not os.path.exists("ovl"):
        os.mkdir("ovl")
    for combo in combinations(Map_contigs.keys(), 2):
        l_map, r_map = combo
        l_contig = Map_contigs[l_map]
        r_contig = Map_contigs[r_map]
        if l_map[1] > r_map[0] and l_map[0] < r_map[1]:
            ovl_len = min(l_map[1], r_map[1])- r_map[0]
            l_ovl = [l_contig[1] - ovl_len, l_contig[1], l_contig[2], l_contig[3]]
            r_ovl = [r_contig[0], r_contig[0] + ovl_len, r_contig[2], r_contig[3]]

            for ovl in [l_ovl, r_ovl]:
                start_on_ctg = ovl[0]
                end_on_ctg = ovl[1]
                ctg_name = ovl[2]
                orientation = ovl[3]

                sequence = contigsSequencesDict[ctg_name][start_on_ctg:end_on_ctg]
                if orientation == '-1':
                    sequence = revcom(contigsSequencesDict[ctg_name][end_on_ctg:start_on_ctg:-1])
                ovl_name = "{}-{}_{}".format(r_map[0], l_map[1], ctg_name)
                with open("ovl/{}.fasta".format(ovl_name), "w") as overlap:
                    overlap.write(">{}\n".format(ovl_name))
                    overlap.write(sequence)
                    

def _compute_assembly_stats(assembly, genomeSize):
    stats_file_name = ".".join([assembly, "statistics", "txt"])
    contigsLengthDict = {}
    contigsSequencesDict = {}
    contigsLength = []
    totalLength = 0
    numContigs = 0
    maxContigLength = 0
    N50 = 0
    N80 = 0
    numNs = 0
    
    with open(assembly, "r") as ref_fd:
        fasta_header = ref_fd.readline()
        header = (fasta_header.split(">")[1]).rstrip()
        sequence = ""
        for line in ref_fd:
            if line.startswith(">"):
                contigsLength.append(len(sequence))
                contigsLengthDict[header] = len(sequence)
                contigsSequencesDict[header] = sequence
                totalLength += len(sequence)
                numContigs += 1
                counter = Counter(sequence)
                numNs += (counter['n'] + counter['N'])
                sequence = ""
                header = (line.split(">")[1]).rstrip()
            else:
                sequence+=line.rstrip()
        contigsLength.append(len(sequence))
        contigsLengthDict[header] = len(sequence)
        contigsSequencesDict[header] = sequence
        totalLength += len(sequence)
        numContigs  += 1
        counter = Counter(sequence)
        numNs += (counter['n'] + counter['N'])

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
    match_expr = re.compile("(\S+)\.(scf|ctg)\.(fasta|fa)$", re.IGNORECASE)
    new_assembly_name = os.path.basename(assembly)
    match_res = re.match(match_expr, new_assembly_name)
    try:
        (basename, type, extention) = match_res.groups()
    except ValueError, AttributeError:
        print("The assembly name is not valid: {}".format(new_assembly_name))
        return new_assembly_name

    new_assembly_name = basename + ".{}.{}.{}".format(type, minCtgLength, extention)
    new_assembly_name =  os.path.abspath(os.path.basename(new_assembly_name))

    if os.path.exists(new_assembly_name):
        print "assembly {} already created".format(new_assembly_name)
        return new_assembly_name
    
    contigID = 1
    with open(new_assembly_name, "w") as new_ref_fd:
        with open(assembly, "r") as ref_fd:
            fasta_header = ref_fd.readline()
            sequence = ""
            for line in ref_fd:
                line = line
                if line.startswith(">"):
                    if len(sequence) >= minCtgLength:
                        new_ref_fd.write(">contig{}\n".format(contigID))
                        new_ref_fd.write(sequence)
                    sequence = ""
                    contigID += 1
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
    parser.add_argument('--only-reference', action='store_true', default = False,
            help="generates only the new refernce (i.e., only scaffolds longer than min-contig)")
    parser.add_argument('--only-stats', action='store_true', default = False,
            help="generates only assembly stats")
    parser.add_argument('--opgen-report',help="op gen placment report", type=str)
    parser.add_argument('--find-problems-in-maps',
            help="generates only the new refernce (i.e., only scaffolds longer than min-contig)",
            default=0, type=int)
    parser.add_argument('--produce-consensus',
            help="generates only the new refernce (i.e., only scaffolds longer than min-contig)",
            default=0, type=int)
    parser.add_argument('--output', help="aoutput header name",
            default="opgen_scaffolded_assembly", type=str)
    parser.add_argument('--place-last', 
            help="Place contigs with this prefix last on the consensus sequence", type=str)
    args = parser.parse_args()
    if args.genome_size == 0:
        print "genome size must be specified"
        sys.exit("error")

    main(args)

