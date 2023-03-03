import os, sys
from pathlib import Path
import argparse
import functools
import textwrap

# establish arguments
parser = argparse.ArgumentParser()

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage = "annotateBlast.py blastOut geneNames -o outputPrefix",
                                 description = textwrap.dedent('''\
                                    Example:
                                    python annotateBlast.py blastOut.txt geneNames.txt -o sample

                                    The blast output text file should be outfmt 6 with the standard columns
                                    the geneList text file should have the names of the genes as obtained from the fasta headers
                                    of the reference RNA.fa file from NCBI, with the ">" symbol removed
                                    '''))

parser.add_argument("blastOut", type=str,
                    help="results of blast alignment of trimmed inserts against reference genome, in tab-delimited outfmt 6")

parser.add_argument("geneNames", type=str,
                    help="list of genes from reference RNA.fa files from NCBI")

parser.add_argument("-o", "--outputPrefix", type=str, default="out",
                    help="identifier for output annotated text file")

args = parser.parse_args()

#open reference and destination files from input arguments
blastF = open(args.blastOut,'r')
namesF = open(args.geneNames,'r')

#read all lines from both files into memory
blastList = blastF.readlines()
namesList = namesF.readlines()

blastF.close()
namesF.close()

annoList = []
outFile = open(args.outputPrefix + ".ann.txt",'a')

for alignment in blastList :
    blastLine = alignment.rstrip().split("\t")
    blastGene = blastLine[1]
    for rna in namesList :
        rnaAnno = rna.split(" ")
        rnaGene = rnaAnno[0]
        rnaName = " ".join(rnaAnno[5:]).rstrip()
        if blastGene == rnaGene :
            blastLine.append(rnaName)
    annoList.append(blastLine)
    outLine = "\t".join(blastLine)
    outFile.write(outLine + "\n")

outFile.close()
