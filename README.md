# permutation-finder
Telomere permutation finder for Chen lab

## Installation Notes:

This package was written and tested with Python 3.10.4 in a Linux environment.

I recommend setting up a dedicated conda environment as general best practice, but the required packages are fairly standard and you may already have them installed in your base environment.

## Download:

First, go to the directory into which you want to install the package. Then, download the release archive from GitHub, unpack, and run.
```
wget https://github.com/kristina28/permutation-finder/archive/refs/tags/v1.0.1.tar.gz
tar -xzvf v1.0.1.tar.gz
```
(or whatever is the latest release available)

If setup.py doesn't have executable privileges, run the following code to set that up.
```
chmod +x setup.py
```
Then, run the following lines to quickly install the code into your Python site-packages (will be environment-specific)
```
python setup.py build
python setup.py install
```
## Dependencies:

pathlib v1.0.1
```
conda install -c anaconda pathlib
```
biopython v1.79
```
conda install -c conda-forge biopython
```
## Quick Start to Running the Code:

### Option 1: Individual Permutation

In this mode, the script takes a permutation sequence and the genome of interest and locates all appearances of the permutation in the genome. Restricted characters preceding and following the permutation sequence can be provided to filter the returned hits.

If the -t flag is specified, a tabular summary of the hits is generated. If the -f flag is specified, a fasta file containing the hits is generated. The --len1 and --len2 flags provide a way to set the number of bases before and after each hit to include in the fasta file.

#### Example:
```
python teloPermutations.py genome.fa --permutation "TAACCCAAGTA" --len1 200 --len2 800 --char1 "g" --char2 "a" --prefix "telo1" -t -f
```
#### Expected Output:

**telo1.summary.txt**
```
ID	        Template_permutation	  Permutation_nt    Match_sequence	  Chr	            Match_start	  Match_end	Strand
test-perm_0000	[-g]TAACCCAAGTA[-a]	  11	            ATAACCCAAGTAC	  NW_019640800.1	302828	    303839	  +
test-perm_0001	[-g]TAACCCAAGTA[-a]	  11	            ATAACCCAAGTAC	  NW_019640811.1	117771	    118782	  +
test-perm_0002	[-g]TAACCCAAGTA[-a]	  11	            TTAACCCAAGTAC	  NW_019640833.1	373691	    374702	  -
...
```

**telo1.fasta**

```
>test-perm_0000|[-g]TAACCCAAGTA[-a]|ATAACCCAAGTAC|11|NW_019640800.1|303428-304439(+)|L=1011
catacacacacacacacacacacacacacgcgtatgcaCGCGTGGGGGGAAGGGGGGTTAACAACGAAGTCGCACACATCGTTGCAGTCGGTAGCGCGAGTACATCAAAACAATAGCAATACCTCGTCATCGTCGAAGGGCCGttggtatatatgtgtaccaCTGCATCGTGACGCCTACgccgctatatacacgcgaaatAACCCAAGTACTTGACATTCCTATACGCCTGTATGGAACATGCATTGCAGGCGCTGAGCTCGTGTGAACACAcctataatacacatacataaacaaCGCAGCCGACCGATGCGcgcacgagtatatataagcgagggCTGCAGGGCGAGAGGTGCGACATGCTCACACTCTGACACGTACAGAGCGCGAATCGATACacaaagacacacacacacatgcaccaCCGATGGCCATTGCTAATCACAGCGGACTCTAATCAGCGCCGgactagctgctgctgctgcctcgcttGCTAAAGCTCTATCTCTCTACGTACTAATTGCCATACCAAGGACCCTACTCGCAAGTTGCGAAAATCCAAATCAAGACATTTCCACTCGTGCGtccatgtgtgcgtgtgtatatacgtgttatgtaccactgctgctgctgcaatatTCGGACTATATTATTACTATCGCCGCgccgcgcaagctcgcaacaACGTATTCGGAGTATCGGTATTCGGTTAATTATGATCCTcggctcggcgcgcgcgcgcaagagagagaccGTCGGTGTTATTGCTAGAGCACAGCAGCCCCCGCCGAATAACTCGGTTCCCATGTGTTGTATCGTCGAATTCTTCGACACCCCCCTTCCTGCGTTATATTATCTATGTGCCATATTGGCACGCCGCACACATGGCAGTCCGATGTGTCCCGGCAATTACTACCATGTCCATGATCTAGCTCTCATTTACGCGGCTGAGTATGCACTGCGTACATCTCTCTATACCTACCTGTTTGTTATTAGTGCAATAACTGAATA
...
```

### Option 2: All Permutations for a Telomere Pattern

In this mode, the script takes a telomere repeat pattern and generates all possible permutations and associated restricted characters from 2bp longer than the pattern to twice the pattern length. It then finds all hits in the genome for every unique permutation, as in option 1. Tabular and fasta output can be generated as individual files for each permutation, combined files for all permutations of each unique length, or combined files for all possible permutations together. Either way, the IDs provided for each hit are unique.

#### Example:
```
python teloPermutations.py genome.fa --pattern "TTACTTGGG" --len1 200 --len2 800 --prefix "telo-pattern" -t -f -m -s -a
```
#### Expected Output:

_If -s flag is given:_

telo-pattern.CCCAAGTAACC.summary.txt and telo-pattern.CCCAAGTAACC.fasta

telo-pattern.CCAAGTAACCC.summary.txt and telo-pattern.CCAAGTAACCC.fasta

...

telo-pattern.AACCCAAGTAACCCAAGT.summary.txt and telo-pattern.AACCCAAGTAACCCAAGT.fasta

telo-pattern.ACCCAAGTAACCCAAGTA.summary.txt and telo-pattern.ACCCAAGTAACCCAAGTA.fasta

_If -m flag is given:_

telo-pattern.11.summary.txt and telo-pattern.11.fasta

telo-pattern.12.summary.txt and telo-pattern.12.fasta

...

telo-pattern.17.summary.txt and telo-pattern.17.fasta

telo-pattern.18.summary.txt and telo-pattern.18.fasta

_If -a flag is given:_

telo-pattern.all.summary.txt and telo-pattern.all.fasta
