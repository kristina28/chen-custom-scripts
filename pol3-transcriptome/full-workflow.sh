#!/bin/bash

VALID_ARGS=$(getopt -o i:b:f:o:s:e:r:d:a: \
                    --long input:,barcodes:,filters:,output:,seqtype:,errorMargin:,reference:,database:,adaptors: \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

declare -A varmap
varmap[SCP_DIR]=$(pwd)
varmap[output]="out"
varmap[seqtype]="fastq"
varmap[errorMargin]=2

host=$(hostname)
if [[ "$host" == *"agave"* ]]; then
  varmap[PART]="htc"
  varmap[QOS]="normal"
elif [[ "$host" == *"sol"* ]]; then
  varmap[PART]="highmem"
  varmap[QOS]="public"
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -i | --input)
        echo "Input reads are the file '$2'"
        varmap[input]="$2"
        shift 2
        ;;
    -b | --barcodes)
        echo "Barcodes for demultiplexing are in the file '$2'"
        varmap[barcodes]="$2"
        shift 2
        ;;
    -f | --filters)
        echo "Reads containing sequences in '$2' will be filtered out"
        varmap[filters]="$2"
        shift 2
        ;;
    -o | --output)
        echo "Results files will have the prefix '$2'"
        varmap[output]="$2"
        shift 2
        ;;
    -s | --seqtype)
        echo "Input reads are in the format '$2'"
        varmap[seqtype]="$2"
        shift 2
        ;;
    -e | --errorMargin)
        echo "A passing alignment score will be no less than '$2' less than the target length"
        varmap[errorMargin]="$2"
        shift 2
        ;;
    -r | --reference)
        echo "Reference database for blast can be found in the directory '$2'"
        varmap[refdir]="$2"
        shift 2
        ;;
    -d | --database)
        echo "Reference database for blast has the prefix '$2'"
        varmap[database]="$2"
        shift 2
        ;;
    -a | --adaptors)
        echo "Comma-separated forward and reverse adaptor sequence"
        varmap[adFor]=$(echo "$2" | cut -d "," -f 1)
        varmap[adRev]=$(echo "$2" | cut -d "," -f 2)
        shift 2
        ;;
    --)
        shift;
        break
        ;;
    *)
        echo "Unexpected option: $1 - please correct."
        ;;
  esac
done

declare -p varmap > params.txt

if [ -z "${varmap[refdir]}" ]; then
  echo "Please specify a reference directory"
elif [ -z "${varmap[database]}" ]; then
  echo "Please provide your blast database prefix"
elif [ -z "${varmap[input]}" ]; then
  echo "Please specify an input file containing Nanopore reads"
elif [ -z "${varmap[barcodes]}" ]; then
  echo "Please specify an input file containing sample names and barcodes for demultiplexing"
elif [ -z "${varmap[filters]}" ]; then
  echo "Please specify an input file containing short undesired sequences to use as filters"
else
  count=$(wc -l barcodes.txt | grep -o -E '[0-9]+')
  for i in $(seq 0 $((count-1)))
  do
    id=$((i+1))
    varmap[$i]=$(awk "NR==$id" barcodes.txt | cut -d "," -f 1)
    echo "$id ${varmap[$i]}"
  done

  declare -p varmap > params.txt

  filterDemux=$(sbatch --account=$USER -p "${varmap[PART]}" -q "${varmap[QOS]}" \
                       --parsable "${varmap[SCP_DIR]}"/filterAndDemux.sh)
  echo "$filterDemux"
  trim=$(sbatch --dependency=afterok:"$filterDemux" --array=0-$((count-1)) --parsable \
               -p "${varmap[PART]}" -q "${varmap[QOS]}" --account=$USER \
               --export=ALL "${varmap[SCP_DIR]}"/trimAdapters.sh)
  length=$(sbatch --dependency=afterok:"$trim" --array=0-$((count-1)) --parsable \
               -p "${varmap[PART]}" -q "${varmap[QOS]}" --account=$USER \
               --export=ALL "${varmap[SCP_DIR]}"/lengthAnalysis.sh)
  blastAnnotate=$(sbatch --dependency=afterok:"$trim" --array=0-$((count-1)) --parsable \
               -p "${varmap[PART]}" -q "${varmap[QOS]}" --account=$USER \
               --export=ALL "${varmap[SCP_DIR]}"/blastAnnotate.sh)
fi
