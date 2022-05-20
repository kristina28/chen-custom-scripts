from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord

def findPermutations(genome,permutation,len1,len2,char1,char2,prefix,count):
    permArray = []
    with open(genome) as in_handle:
        for title, seq in SimpleFastaParser(in_handle):
            start = 0
            end = len(seq)
            length = len(permutation)
            while permutation.lower() in seq[start:end].lower():
                index = seq.lower().find(permutation.lower(),start,end)
                subseq = seq[(index-len1):(index+length+len2)]
                if ((seq[index-1].lower() != char1.lower()) and (seq[index+length].lower() != char2.lower())):
                    newRecord = [(prefix + "_" + f'{count:04d}'),
                                 ("[-" + char1.lower() + "]" + permutation + "[-" + char2.lower() + "]"),
                                 (seq[index-1].upper() + permutation + seq[index+length].upper()),
                                 str(len(permutation)),
                                 title.split(" ")[0],
                                 (str(index-len1) + "-" + str(index+length+len2) + "(+)"),
                                 ("L=" + str(length + len1 + len2)),
                                 str(index-len2),
                                 str(index+length+len1),
                                 "+",
                                 subseq]
                    permArray.append(newRecord)
                    count += 1
                start = index+length
            permRev = reverse_complement(permutation)
            char1Rev = reverse_complement(char2)
            char2Rev = reverse_complement(char1)
            start = 0
            while permRev.lower() in seq[start:end].lower():
                index = seq.lower().find(permRev.lower(),start,end)
                subseq = seq[(index-len2):(index+length+len1)]
                if ((seq[index-1].lower() != char1Rev.lower()) and (seq[index+length].lower() != char2Rev.lower())):
                    newRecord = [(prefix + "_" + f'{count:04d}'),
                                 ("[-" + char1.lower() + "]" + permutation + "[-" + char2.lower() + "]"),
                                 (reverse_complement(seq[index+length].upper()) + permutation + reverse_complement(seq[index-1].upper())),
                                 str(len(permutation)),
                                 title.split(" ")[0],
                                 (str(index-len2) + "-" + str(index+length+len1) + "(-)"),
                                 ("L=" + str(length + len1 + len2)),
                                 str(index-len2),
                                 str(index+length+len1),
                                 "-",
                                 reverse_complement(subseq)]
                    permArray.append(newRecord)
                    count += 1
                start = index+length
    return([permArray,count])
