# import necessary python packages
import os, sys
from pathlib import Path
import argparse
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from findPermutations import *
from outputHits import *
from patternMaker import *
import functools
import textwrap

# establish arguments
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage = "teloPermutations.py genome [--permutation] [--pattern] [--char1] [--char2] [--len1] [--len2] [--prefix] [-t] [-f] [-s] [-m] [-a] [-h]",
                                 description = textwrap.dedent('''\
                                    Permutation Example:
                                    python teloPermutations.py genome.fna
                                                --permutation "TAACCCAAGTA" --len1 200 --len2 800
                                                --char1 "g" --char2 "a" --prefix "test-perm" -t -f

                                    Pattern Example:
                                    python teloPermutations.py genome.fna
                                                --pattern "TTACTTGGG" --len1 200 --len2 800
                                                --prefix "test-perm" -t -f -s -m -a
                                    '''))
parser.add_argument("genome", type=str,
                    help="fasta file in which to search for desired permutation(s)")

parser.add_argument("--permutation", type=str,
                    help="permutation to search for, excluding restricted characters; either this or pattern are required")

parser.add_argument("--pattern", type=str,
                    help="pattern to search for, excluding restricted characters; either this or permutation are required")

parser.add_argument("--char1", type=str, default="N",
                    help="restricted character preceding permutation, defaults to N")

parser.add_argument("--char2", type=str, default="N",
                    help="restricted character following permutation, defaults to N")

parser.add_argument("--len1", type=int, default=100,
                    help="length of sequence to include preceding permutation, defaults to 100")

parser.add_argument("--len2", type=int, default=100,
                    help="length of sequence to include following permutation, defaults to 100")

parser.add_argument("--prefix", type=str, default="ASUBC",
                    help="unique ID prefix to assign to all hits from the search")

parser.add_argument("-t", "--table", action="store_true",
                    help="outputs tabular summary of identified permutations if specified")

parser.add_argument("-f", "--fasta", action="store_true",
                    help="outputs fasta of identified permutations with specified surrounding regions if specified")

parser.add_argument("-m", "--multiple", action="store_true",
                    help="for use with pattern option only; generates one multi-fasta/table for all hits for each permutation length")

parser.add_argument("-a", "--all", action="store_true",
                    help="for use with pattern option only; generates one multi-fasta/table for all hits for all permutation lengths combined")

parser.add_argument("-s", "--single", action="store_true",
                    help="for use with pattern option only; generates one multi-fasta/table for all hits for each permutation individually")

args = parser.parse_args()

if args.permutation:
    if args.pattern:
        sys.exit("Error: select either --permutation or --pattern but not both")
    else:
        hitsCount = findPermutations(args.genome,
                                args.permutation,
                                args.len1, args.len2,
                                args.char1, args.char2,
                                args.prefix, 0)
        hits = hitsCount[0]
        if args.table:
            hitsTablePerm(hits, args.prefix)
        if args.fasta:
            hitsFastaPerm(hits, args.prefix)
        sys.exit("Success! " + str(len(hits)) + " hits for permutation " + args.permutation + " have been identified")
elif args.pattern:
    patternArray = patternMaker(args.pattern)
    permutationArray = []
    for i in range(len(args.pattern)+2,1+len(args.pattern)*2):
        permutationArray.append([])
    count = 0
    for record in patternArray:
        hitsCount = findPermutations(args.genome,
                            record[0],
                            args.len1, args.len2,
                            record[1], record[2],
                            args.prefix, count)
        hits = hitsCount[0]
        count = hitsCount[1]
        if len(hits) > 0:
            print(hits[0][3])
            permutationArray[int(hits[0][3])-(len(args.pattern)+2)].append(hits)
            print(str(len(hits)) + " hits for permutation " + record[0] + " have been identified")
    if args.table:
        hitsTablePattern(permutationArray, args.prefix, args.single, args.multiple, args.all)
    if args.fasta:
        hitsFastaPattern(permutationArray, args.prefix, args.single, args.multiple, args.all)
else:
    sys.exit("Error: must select at least one of --permutation and --pattern")
