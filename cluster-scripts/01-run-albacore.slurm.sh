#!/bin/bash
#SBATCH -n 28
#SBATCH --mem=96000
#SBATCH -p himemory
#SBATCH -w compute-2-0
#SBATCH --mail-type=END,FAIL
#SBATCH -J albacore 

### Strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'
set -x

### Check yer variables (probably redundant with the above but good practice)
# if [ -z ${var+x} ]; then echo "var is unset"; else echo "var is set to '$var'"; fi

### Load Modules
module load albacore/2.3.1-IGB-gcc-4.9.4-Python-3.6.1 
### Run app on file

SCRATCH=/scratch/${USER}-albacore

cd $SCRATCH

read_fast5_basecaller.py \
    -i FAST5/ \
    -t $SLURM_NPROCS \
    -s results \
    -k SQK-LSK108 \
    -f FLO-MIN106 \
    --recursive \
    -o fast5,fastq

