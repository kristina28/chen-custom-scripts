from Bio.Seq import Seq, reverse_complement

def patternMaker(pattern):
    pArray = []
    count = 2
    while count <= len(pattern):
        print(count)
        pLen = len(pattern)
        seq = reverse_complement(pattern)
        seq1 = seq + seq[0:count]
        iter = 0
        while iter < pLen:
            char1 = seq1[len(seq1) - (count+1)]
            char2 = seq1[count]
            pArray.append([seq1, char1, char2])
            seq1 = seq1[1:] + char2
            iter += 1
        count += 1
    return(pArray)
