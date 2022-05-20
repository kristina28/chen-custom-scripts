import os, sys
from pathlib import Path
from Bio import SeqIO

# this function takes the array of id and sequence data from findPermutations and creates a tabular summary
def hitsTablePerm(data, prefix):
    header = "ID\tTemplate_permutation\tPermutation_nt\tMatch_sequence\tChr\tMatch_start\tMatch_end\tStrand"
    with open(os.path.basename(prefix + ".summary.txt"), 'w') as textfile:
        textfile.write(header + "\n")
        textfile.writelines("\t".join([i[x] for x in [0,1,3,2,4,7,8,9]]) + "\n" for i in data)


# this function creates a multi-fasta file for all permutation hits in the provided array
def hitsFastaPerm(data, prefix):
    with open(os.path.basename(prefix + ".fasta"), 'w') as out_handle:
        for i in data:
            out_handle.write(">%s\n%s\n" % ("|".join(i[:7]), i[10]))

def hitsTablePattern(data, prefix, single, multiple, all):
    header = "ID\tTemplate_permutation\tPermutation_nt\tMatch_sequence\tChr\tMatch_start\tMatch_end\tStrand"
    if single:
        for lengthArray in data:
            for hitsArray in lengthArray:
                length = int(hitsArray[0][3])
                perm = hitsArray[0][1][4:(4+length)]
                with open(os.path.basename(prefix + "." + perm + ".summary.txt"), 'w') as textfile:
                    textfile.write(header + "\n")
                    textfile.writelines("\t".join([i[x] for x in [0,1,3,2,4,7,8,9]]) + "\n" for i in hitsArray)
    if multiple:
        for lengthArray in data:
            length = lengthArray[0][0][3]
            with open(os.path.basename(prefix + "." + str(length) + ".summary.txt"), 'w') as textfile:
                textfile.write(header + "\n")
                for hitsArray in lengthArray:
                    textfile.writelines("\t".join([i[x] for x in [0,1,3,2,4,7,8,9]]) + "\n" for i in hitsArray)
    if all:
        with open(os.path.basename(prefix + ".all.summary.txt"), 'w') as textfile:
            textfile.write(header + "\n")
            for lengthArray in data:
                for hitsArray in lengthArray:
                    textfile.writelines("\t".join([i[x] for x in [0,1,3,2,4,7,8,9]]) + "\n" for i in hitsArray)

def hitsFastaPattern(data, prefix, single, multiple, all):
    if single:
        for lengthArray in data:
            for hitsArray in lengthArray:
                length = int(hitsArray[0][3])
                perm = hitsArray[0][1][4:(4+length)]
                with open(os.path.basename(prefix + "." + perm + ".fasta"), 'w') as out_handle:
                    for i in hitsArray:
                        out_handle.write(">%s\n%s\n" % ("|".join(i[:7]), i[10]))
    if multiple:
        for lengthArray in data:
            length = lengthArray[0][0][3]
            with open(os.path.basename(prefix + "." + length + ".fasta"), 'w') as out_handle:
                for hitsArray in lengthArray:
                    for i in hitsArray:
                        out_handle.write(">%s\n%s\n" % ("|".join(i[:7]), i[10]))
    if all:
        with open(os.path.basename(prefix + ".all.fasta"), 'w') as out_handle:
            for lengthArray in data:
                for hitsArray in lengthArray:
                    for i in hitsArray:
                        out_handle.write(">%s\n%s\n" % ("|".join(i[:7]), i[10]))
