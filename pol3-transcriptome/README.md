## simple pipeline
The most basic way to run this workflow is as a series of Python scripts.

#### 1. Filter and Demultiplex
The script filterAndDemux.py is used to separate out samples by their barcode sequences as well as filter out any reads containing an unwanted sequence, such as a poly(A) sequence for removing mRNA reads. 

usage: filterAndDemux.py input -b barcodes.csv -f filters.txt

Example:
  barcodes.csv:
    Sample,AAGAAAGTTGTCGGTGTCTTTGTG
    Control,GAGTCTTGTGTCCCAGTTACCAGG
  filters.txt:
    AAAAAAAAAAAAAAAAAAAA

  filterAndDemux.py reads.fastq -b barcodes.csv -f filters.txt -o output -e 2 -s fastq

positional arguments:
  input                 Nanopore fastq file

optional arguments:
  -h, --help            show this help message and exit
  -b BARCODES, --barcodes BARCODES
                        comma-separated file containing sample name and
                        barcode sequence in pairs
  -f FILTERSEQS, --filterSeqs FILTERSEQS
                        text file containing list of sequences to filter out
                        of the data
  -o OUTPUTPREFIX, --outputPrefix OUTPUTPREFIX
                        identifier for output annotated text file
  -e ERRORMARGIN, --errorMargin ERRORMARGIN
                        identifier for output annotated text file
  -s SEQTYPE, --seqtype SEQTYPE
                        identifier for output annotated text file
                        
#### 2. Trim Adapters
The script trimAdapters.py is used to remove adapter sequences from both ends of the Nanopore reads. This will need to be run for each demultiplexed sample individually. Because it trims on the sequence between the barcode and the insert, the identifying index information will be removed so the reads must be already demultiplexed.

usage: trimAdapters.py readsFile -f forwardAdaptor -r reverseAdaptor -o outputPrefix

Example:
trimAdapters.py reads.fastq -f 'ACTTGCCTGTCGCTCTATCTTC' -r 'TTTCTGTTGGTGCTGATATTGC'

both adaptor strings should be in 5' to 3' orientation.
If not specified, the sequences above are the defaults.
In the case of a barcoded adaptor, use the sequence between the barcode and insert.

positional arguments:
  readsFile             fastq file containing Nanopore reads

optional arguments:
  -h, --help            show this help message and exit
  -f FORWARDADAPTOR, --forwardAdaptor FORWARDADAPTOR
                        string specifying the forward adaptor sequence in 5'
                        to 3' orientation
  -r REVERSEADAPTOR, --reverseAdaptor REVERSEADAPTOR
                        string specifying the reverse adaptor sequence in 5'
                        to 3' orientation
  -o OUTPUTPREFIX, --outputPrefix OUTPUTPREFIX
                        string providing a prefix for output file naming
                        
#### 3. Length Analysis
This script can be run now, before annotation and subtraction of control reads, or it can be run following those steps, or it can be run at both times. In the SLURM workflow it is run at this step.

usage: lengthAnalysis.py readsFile -o outputPrefix

Example:
lengthAnalysis.py reads.fasta -o sample

positional arguments:
  readsFile             fasta file containing trimmed insert sequences from
                        Nanopore reads

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUTPREFIX, --outputPrefix OUTPUTPREFIX
                        string providing a prefix for output file naming
                        
This script will evaluate the sequence length of all the trimmed reads for each sample, print histogram values and generate a corresponding PNG image of the histogram. It also will need to be run for each sample in the data set.

#### 4. Blast Sequences and Annotate Reads
Before running the annotation Python script, you'll need to run some short code snippets to actually blast the input data and sort the blast output. Because we're working with the transcriptome of a single organism, it's easy to use a small local blast database built from the reference genome for that organism.
```
# run nucleotide blast on trimmed sequences
blastn -task blastn-short -db /path/to/reference/blastdb \
       -query query.fasta -out output -outfmt 6 -max_target_seqs 1 -culling_limit 1

# sort files by query IDs and gene names
sort -k2 -n output >> sorted.output
```

The output from that code can be input to annotateBlast.py. You'll also need a list of transcript names from the reference RNA for your organism so the NCBI IDs can be cross-referenced.

usage: annotateBlast.py blastOut geneNames -o outputPrefix

Example:
python annotateBlast.py blastOut.txt geneNames.txt -o sample

The blast output text file should be outfmt 6 with the standard columns
the geneList text file should have the names of the genes as obtained from the fasta headers
of the reference RNA.fa file from NCBI, with the ">" symbol removed

positional arguments:
  blastOut              results of blast alignment of trimmed inserts against
                        reference genome, in tab-delimited outfmt 6
  geneNames             list of genes from reference RNA.fa files from NCBI

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUTPREFIX, --outputPrefix OUTPUTPREFIX
                        identifier for output annotated text file
                        
#### 5. Subtract Negative Control Reads from Sample Reads
For this step you just need a few lines of bash to remove from the sample reads any annotation which occurs in the control, and to filter out the 40S ribsomal transcripts which don't have the poly(A) tail but are not transcribed by RNA Pol III. Any other known PolI transcript can be filtered out in this same way!
```
grep -vw "40S ribosomal protein" Sample.paired.fasta.blast.ann.txt > Sample.filtered.paired.fasta.blast.ann.txt

sort -k2 Control.paired.fasta.blast.ann.txt > Control.sorted.paired.fasta.blast.ann.txt
sort -k2 Sample.filtered.paired.fasta.blast.ann.txt > Sample.sorted.filtered.paired.fasta.blast.ann.txt

join -v 1 -1 2 -2 2 Sample.sorted.filtered.paired.fasta.blast.ann.txt Control.sorted.paired.fasta.blast.ann.txt > Sample.normalized.blast.ann.txt

awk '{ print $2 }' Sample.normalized.blast.ann.txt > normalized-sample-reads.txt

#use the list of sequences to subset the trimmed fasta file and generate a NTC-normalized fasta file
seqtk subseq trimmed-Sample.fa normalized-sample-reads.txt
```
At this point you can rerun the length analysis on the normalized fasta if you would like for final histogram information.

## SLURM Pipeline
If you're using the ASU slurm-based supercomputer cluster, either Agave or Sol, you can call the full-workflow.sh script and steps 1-4 above will run automatically for you. You will still need to go through step 5 manually.

An example:

./full-workflow.sh --input fastq_runid_6a85d64cdc1a33f7d9cd2ae8e7fb80b3f63ad0ed_cat.fastq \
                   --barcodes barcodes.txt --filters polya.txt --seqtype fastq \
                   --output fn --reference /data/biocore/reference_genomes/cyanidioschyzon_merolae \
                   --database c_merolae --adaptors ACTTGCCTGTCGCTCTATCTTC,TTTCTGTTGGTGCTGATATTGC
