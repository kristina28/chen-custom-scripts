#!/bin/bash

##### trim filtered and demultiplexed Nanopore reads using provided adapter sequence #####

#SBATCH -o slurm.%j.trim.python.txt
#SBATCH -t 0-8:00
#SBATCH -c 1
#SBATCH --mem=16G

module load mamba/latest

source deactivate

source activate /data/biocore/programs/mamba-envs/nanopore-env

source params.txt
index=$SLURM_ARRAY_TASK_ID
sid=${varmap[$index]}

sampleFile="${varmap[output]}".nanoreads."$sid"."${varmap[seqtype]}"

echo $sampleFile

python3 trimAdapters.py "$sampleFile" -o "${varmap[output]}"."$sid" \
                        -f "${varmap[adFor]}" -r "${varmap[adRev]}"
