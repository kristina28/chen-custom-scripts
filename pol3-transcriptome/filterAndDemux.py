#replaces excludePolyTails.py in the pipeline

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
import os, sys
from pathlib import Path
import argparse
import functools
import textwrap

# establish arguments
parser = argparse.ArgumentParser()

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage = "filterAndDemux.py input -b barcodes.csv -f filters.txt",
                                 description = textwrap.dedent('''\
                                    Example:
                                      barcodes.csv:
                                        Sample,AAGAAAGTTGTCGGTGTCTTTGTG
                                        Control,GAGTCTTGTGTCCCAGTTACCAGG
                                      filters.txt:
                                        AAAAAAAAAAAAAAAAAAAA

                                      filterAndDemux.py reads.fastq -b barcodes.csv -f filters.txt \
                                                        -o output -e 2 -s fastq
                                    '''))

parser.add_argument("input", type=str,
                    help="Nanopore fastq file")

parser.add_argument("-b", "--barcodes", type=str,
                    help="comma-separated file containing sample name and barcode sequence in pairs")

parser.add_argument("-f", "--filterSeqs", type=str,
                    help="text file containing list of sequences to filter out of the data")

parser.add_argument("-o", "--outputPrefix", type=str, default="out",
                    help="identifier for output annotated text file")

parser.add_argument("-e", "--errorMargin", type=int, default=2,
                    help="identifier for output annotated text file")

parser.add_argument("-s", "--seqtype", type=str, default="fastq",
                    help="identifier for output annotated text file")

args = parser.parse_args()

def filter(readsArray, target, margin) :
    target = Seq(target)
    keep = []
    discard = []
    for sequence in readsArray :
        #print("sequence is " + sequence.seq)
        #print("target is " + target)
        pf = passFilter(sequence.seq, target, margin)
        #print(pf)
        if pf :
            keep.append(sequence)
        else :
            discard.append(sequence)
    return (keep, discard)

def passFilter(sequence, target, margin) :
    target = Seq(target)
    pf = False
    score = pairwise2.align.localxs(sequence, target, -1, -1, score_only=True)
    #print("alignment score is: " + str(score))
    if score <= len(target) - margin :
        pf = True
    return pf

def demux(sequence, barcodes, margin) :
    scoreMatrix = {}
    for barcode in barcodes :
        barcode = barcode.split(",")
        barcode[1] = Seq(barcode[1])
        score = pairwise2.align.localxs(sequence.seq, barcode[1], -1, 0, score_only=True)
        rcScore = pairwise2.align.localxs(sequence.seq, barcode[1].reverse_complement(), -1, 0, score_only=True)
        #print("score is ", score, " and rcScore is ", rcScore)
        maxScore = max(score, rcScore)
        scoreMatrix[barcode[0]] = maxScore
    maxBarcode = max(scoreMatrix, key=scoreMatrix.get)
    #print(maxBarcode)
    #print(scoreMatrix[maxBarcode])
    if scoreMatrix[maxBarcode] >= len(barcode[1]) - margin :
        return maxBarcode
    else :
        return "error"

keeps = list(SeqIO.parse(args.input, args.seqtype))
targetFile = open(args.filterSeqs, 'r')
targets = targetFile.readlines()
targetFile.close()

barFile = open(args.barcodes, 'r')
barcodes = barFile.readlines()
barFile.close()

#print(keeps)
discards = []
print("there are " + str(len(keeps)) + " sequences remaining before filtering")
for target in targets :
    filtered = filter(keeps, target, args.errorMargin)
    keeps = filtered[0]
    discards = discards + filtered[1]

print("there are " + str(len(keeps)) + " sequences remaining after filtering")
#print(keeps)

demuxedFasta = {}
for barcode in barcodes :
    barcode = barcode.split(",")
    demuxedFasta[barcode[0]] = []

print(demuxedFasta)

for sequence in keeps :
    assignment = demux(sequence, barcodes, args.errorMargin)
    #print("sequence is ", sequence)
    #print("sample assignment is " + assignment)
    if assignment in demuxedFasta :
        #print(demuxedFasta[assignment])
        old = demuxedFasta[assignment]
        #print("old value is ", old)
        old.append(sequence)
        #print("new value is ", old)
        demuxedFasta[assignment] = old
        #print(demuxedFasta[assignment])
print(type(demuxedFasta))

for assignment in demuxedFasta :
    filename = str(args.outputPrefix + ".nanoreads." + assignment + "." + args.seqtype)
    print(filename)
    print("there are ", len(demuxedFasta[assignment]), " sequences in sample", assignment)
    print(demuxedFasta[assignment][0])
    with open(filename, "w") as output_handle :
        SeqIO.write(demuxedFasta[assignment], output_handle, args.seqtype)
