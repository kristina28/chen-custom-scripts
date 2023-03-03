#!/bin/bash

##### trim filtered and demultiplexed Nanopore reads using provided adapter sequence #####

#SBATCH -o slurm.%j.lengthAnalysis.python.txt
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

python3 lengthAnalysis.py "$sampleFile" -o "${varmap[output]}"."$sid" \
