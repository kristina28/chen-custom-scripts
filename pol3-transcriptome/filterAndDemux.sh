#!/bin/bash

#SBATCH -o slurm.%A.filterDemux.txt               # STDOUT (%A = JobId)
#SBATCH -t 0-4:00
#SBATCH -c 1
#SBATCH --mem=16G

module load mamba/latest

source deactivate

source activate /data/biocore/programs/mamba-envs/nanopore-env

source params.txt

python3 filterAndDemux.py "${varmap[input]}" \
                          --barcodes "${varmap[barcodes]}" \
                          --filterSeqs "${varmap[filters]}" \
                          --outputPrefix "${varmap[output]}" \
                          --errorMargin "${varmap[errorMargin]}" \
                          --seqtype "${varmap[seqtype]}"
