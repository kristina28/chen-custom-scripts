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

# establish arguments
parser = argparse.ArgumentParser()

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage = "trimAdapters.py readsFile -f forwardAdaptor -r reverseAdaptor -o outputPrefix",
                                 description = textwrap.dedent('''\
                                    Example:
                                    trimAdapters.py reads.fastq -f 'ACTTGCCTGTCGCTCTATCTTC' -r 'TTTCTGTTGGTGCTGATATTGC'

                                    both adaptor strings should be in 5' to 3' orientation.
                                    If not specified, the sequences above are the defaults.
                                    In the case of a barcoded adaptor, use the sequence between the barcode and insert.
                                    '''))

parser.add_argument("readsFile", type=str,
                    help="fastq file containing Nanopore reads")

parser.add_argument("-f", "--forwardAdaptor", type=str, default="ACTTGCCTGTCGCTCTATCTTC",
                    help="string specifying the forward adaptor sequence in 5' to 3' orientation")

parser.add_argument("-r", "--reverseAdaptor", type=str, default="TTTCTGTTGGTGCTGATATTGC",
                    help="string specifying the reverse adaptor sequence in 5' to 3' orientation")

parser.add_argument("-o", "--outputPrefix", type=str, default="out",
                    help="string providing a prefix for output file naming")

args = parser.parse_args()

keeps = []
forOnly = []
discards = []

errorMargin = 2
fLength = len(args.forwardAdaptor)
rLength = len(args.reverseAdaptor)
rcForAdaptor = Seq(args.forwardAdaptor).reverse_complement()
rcRevAdaptor = Seq(args.reverseAdaptor).reverse_complement()

with open(args.readsFile) as input_handle :
    for sequence in SeqIO.parse(input_handle, "fastq") :
        fAlign = pairwise2.align.localxs(sequence.seq, args.forwardAdaptor, -1, -0.5, one_alignment_only=True)
        if len(fAlign) > 0 :
            #print("adaptor1 aligned to sequence")
            #print(format_alignment(*fAlign[0]))
            #print(fAlign[0][2], fAlign[0][3], fAlign[0][4])
            if fAlign[0][2] >= fLength - errorMargin :
                rcRAlign = pairwise2.align.localxs(sequence.seq[fAlign[0][4]+1:], rcRevAdaptor, -1, -0.5, one_alignment_only=True)
                if len(rcRAlign) > 0 :
                    #print("reverse complement of adaptor2 aligned to sequence")
                    #print(format_alignment(*rcRAlign[0]))
                    #print(rcRAlign[0][2], rcRAlign[0][3], rcRAlign[0][4])
                    if rcRAlign[0][2] >= rLength - errorMargin :
                        insert = sequence[fAlign[0][4]+1:(rcRAlign[0][3]+fAlign[0][4]+1)]
                        keeps.append(insert)
                else :
                    insert = sequence[fAlign[0][4]+1:]
                    forOnly.append(insert)
        else :
            rAlign = pairwise2.align.localxs(sequence.seq, args.reverseAdaptor, -1, -0.5, one_alignment_only=True)
            if len(rAlign) > 0 :
                #print("adaptor2 aligned to sequence")
                #print(format_alignment(*rAlign[0]))
                #print(rAlign[0][2], rAlign[0][3], rAlign[0][4])
                if rAlign[0][2] >= rLength - errorMargin :
                    rcFAlign = pairwise2.align.localxs(sequence.seq[rAlign[0][4]+1:], rcForAdaptor, -1, -0.5, one_alignment_only=True)
                    if len(rcFAlign) > 0 :
                        #print("reverse complement of adaptor1 aligned to sequence")
                        #print(format_alignment(*rcFAlign[0]))
                        #print(rcFAlign[0][2], rcFAlign[0][3], rcFAlign[0][4])
                        if rcFAlign[0][2] >= fLength - errorMargin :
                            insert = sequence[rAlign[0][4]+1:(rcFAlign[0][3]+rAlign[0][4]+1)]
                            keeps.append(insert)
                    else :
                        insert = sequence[rAlign[0][4]+1:]
                        forOnly.append(insert)
            else :
                discards.append(sequence)

with open((args.outputPrefix + "." + "paired.fasta"), "w") as output_handle :
    SeqIO.write(keeps, output_handle, "fasta")

with open((args.outputPrefix + "." + "forOnly.fasta"), "w") as output_handle :
    SeqIO.write(forOnly, output_handle, "fasta")

with open((args.outputPrefix + "." + "unaligned.fasta"), "w") as output_handle :
    SeqIO.write(discards, output_handle, "fasta")
