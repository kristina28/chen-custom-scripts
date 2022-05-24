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
                if (((index == 0) or (seq[index-1].lower() != char1.lower())) and (((index+length == len(seq)) or seq[index+length].lower() != char2.lower()))):
                    newRecord = makeForRecord(prefix, count, char1, permutation, char2, seq, index, length, title, len1, len2, subseq)
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
                if (((index == 0) or (seq[index-1].lower() != char1Rev.lower())) and ((index+length == len(seq)) or (seq[index+length].lower() != char2Rev.lower()))):
                    newRecord = makeRevRecord(prefix, count, char1, permutation, char2, seq, index, length, title, len1, len2, subseq)
                    permArray.append(newRecord)
                    count += 1
                start = index+length
    return([permArray,count])

def makeRevRecord(prefix, count, char1, permutation, char2, seq, index, length, title, len1, len2, subseq):
    if index == 0:
        edge1 = "N"
    else:
        edge1 = reverse_complement(seq[index-1].upper())
    if index + length >= len(seq):
        edge2 = "N"
    else:
        edge2 = reverse_complement(seq[index+length].upper())
    newRecord = [(prefix + "_" + f'{count:04d}'),
                 ("[-" + char1.lower() + "]" + permutation + "[-" + char2.lower() + "]"),
                 (edge2 + permutation + edge1),
                 str(len(permutation)),
                 title.split(" ")[0],
                 (str(index-len2) + "-" + str(index+length+len1) + "(-)"),
                 ("L=" + str(len(subseq))),
                 str(index-len2),
                 str(index+length+len1),
                 "-",
                 reverse_complement(subseq)]
    return(newRecord)

def makeForRecord(prefix, count, char1, permutation, char2, seq, index, length, title, len1, len2, subseq):
    if index == 0:
        edge1 = "N"
    else:
        edge1 = seq[index-1].upper()
    if index + length >= len(seq):
        edge2 = "N"
    else:
        edge2 = seq[index+length].upper()
    newRecord = [(prefix + "_" + f'{count:04d}'),
                 ("[-" + char1.lower() + "]" + permutation + "[-" + char2.lower() + "]"),
                 (edge1 + permutation + edge2),
                 str(len(permutation)),
                 title.split(" ")[0],
                 (str(index-len1) + "-" + str(index+length+len2) + "(+)"),
                 ("L=" + str(len(subseq))),
                 str(index-len2),
                 str(index+length+len1),
                 "+",
                 subseq]
    return(newRecord)
