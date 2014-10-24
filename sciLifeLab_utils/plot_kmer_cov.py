import sys, os, yaml, glob
import subprocess
import argparse
import pandas as pd
from matplotlib import pyplot as plt


def main(args):

    histogram = args.histogram
    min_limit = args.min_limit
    max_limit = args.max_limit
    kmer = args.kmer
    output_name = args.output_name

    Kmer_histogram = pd.io.parsers.read_csv(histogram, sep='\t', header=None)
    Kmer_coverage  = Kmer_histogram[Kmer_histogram.columns[0]].tolist()
    Kmer_count = Kmer_histogram[Kmer_histogram.columns[1]].tolist()
    Kmer_freq = [Kmer_coverage[i]*Kmer_count[i] for i in \
            range(len(Kmer_coverage))]
    #coverage peak, disregarding initial peak
    kmer_freq_peak = Kmer_freq.index(max(Kmer_freq[min_limit:max_limit]))
    kmer_freq_peak_value=max(Kmer_freq[min_limit:max_limit])

    xmax = max_limit
    ymax = kmer_freq_peak_value + (kmer_freq_peak_value*0.30)

    plt.plot(Kmer_coverage, Kmer_freq)
    plt.title("K-mer length = {}".format(kmer))
    plt.xlim((0,xmax))
    plt.ylim((0,ymax))
    plt.vlines(kmer_freq_peak, 0, kmer_freq_peak_value, colors='r',
            linestyles='--')
    plt.text(kmer_freq_peak, kmer_freq_peak_value+2000, str(kmer_freq_peak))
    plotname = "{}".format(output_name)
    plt.savefig(plotname)
    plt.clf()
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--histogram', type=str, required=True,
            help=".hit file produced by abyss")
    parser.add_argument('--min-limit', type=int, required=True,
            help="lower x-limit")
    parser.add_argument('--max-limit', type=int, required=True,
            help="upper x-limit")
    parser.add_argument('--kmer', type=int, required=True,
            help="kmer employed, if you do not know guess a number")
    parser.add_argument('--output-name', type=str, required=True,
            help="output name (.png)")

    args = parser.parse_args()
    main(args)

