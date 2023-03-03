#!/bin/bash

##### trim filtered and demultiplexed Nanopore reads using provided adapter sequence #####

#SBATCH -o slurm.%j.blastAnnotate.python.txt
#SBATCH -t 0-4:00
#SBATCH -c 1
#SBATCH --mem=16G

module load mamba/latest

source deactivate

source activate /data/biocore/programs/mamba-envs/nanopore-env

source params.txt
index=$SLURM_ARRAY_TASK_ID
sid=${varmap[$index]}

sampleFile="${varmap[output]}"."$sid".paired.fasta

query="$sampleFile"
output="$sampleFile".txt

blastn -task blastn-short -db "${varmap[refdir]}"/"${varmap[database]}" \
       -query "$query" -out "$output" -outfmt 6 -max_target_seqs 1 -culling_limit 1

# sort files by query IDs and gene names
sort -k2 -n "$output" >> sorted."$output"

# match blast results with gene names to provide functional information
python3 annotateBlast.py sorted."$output" c_merolae_gene_names.sorted.txt -o "$sampleFile".blast
