import os, sys
import csv
from pathlib import Path
import argparse
import functools
import textwrap
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import numpy as np
from matplotlib import pyplot as plt

# establish arguments
parser = argparse.ArgumentParser()

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage = "lengthAnalysis.py readsFile -o outputPrefix",
                                 description = textwrap.dedent('''\
                                    Example:
                                    lengthAnalysis.py reads.fasta -o sample'

                                    '''))

parser.add_argument("readsFile", type=str,
                    help="fasta file containing trimmed insert sequences from Nanopore reads")

parser.add_argument("-o", "--outputPrefix", type=str, default="out",
                    help="string providing a prefix for output file naming")

args = parser.parse_args()

lengths = []

with open(args.readsFile) as input_handle :
    for sequence in SeqIO.parse(input_handle, "fasta") :
        lengths.append(len(sequence.seq))

hist, bin_edges = np.histogram(lengths)

print(hist)
print(bin_edges)

plt.hist(hist, bins=bin_edges, log=True)
#lengths.plot.hist(grid=True, bins="auto", rwidth=0.9,
#                   color='#607c8e')
plt.title('Histogram of Insert Lengths for ' + args.outputPrefix)
plt.xlabel('Insert Length')
plt.ylabel('Number of Reads')
plt.grid(axis='y', alpha=0.75)

plt.savefig(args.outputPrefix + '.insert-histogram.png')
