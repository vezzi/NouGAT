import sys, os, yaml, glob
import subprocess
import argparse
import string
import gzip
from subprocess import Popen, PIPE




def main(args):
    workingDir  = os.getcwd()

    proc = Popen(["zcat", "{}".format(args.fastq)], stdout=PIPE, bufsize=1)
    with gzip.open('{}.gz'.format(args.output), 'wb') as output_f:
        for header in iter(proc.stdout.readline, b''):
            header   = header.rstrip()
            sequence = proc.stdout.readline().rstrip()
            comment  = proc.stdout.readline().rstrip()
            quality  = proc.stdout.readline().rstrip()
        
            sequence_rc = rc(sequence)
            quality_rev = reverse(quality)

            output_f.write("{}\n".format(header))
            output_f.write("{}\n".format(sequence_rc))
            output_f.write("{}\n".format(comment))
            output_f.write("{}\n".format(quality_rev))

    proc.communicate()

    return

def reverse(quality):
    seq_rev = quality[::-1]
    return seq_rev


def rc(dna):
    complements = string.maketrans('acgtrymkbdhvnACGTRYMKBDHVN', 'tgcayrkmvhdbnTGCAYRKMVHDBN')
    rcseq = dna.translate(complements)[::-1]
    return rcseq


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reverse complent the fastq file given as input')
    parser.add_argument('--fastq', help="fastq file to reverse complement", type=str)
    parser.add_argument('--output', help="Header name for the output file", type=str)
    args = parser.parse_args()

    main(args)