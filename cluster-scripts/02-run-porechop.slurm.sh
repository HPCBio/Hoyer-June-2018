#!/bin/bash
#SBATCH -n 24 
#SBATCH --mem=48000
#SBATCH -p normal 
#SBATCH --mail-type=END,FAIL
#SBATCH -J porechop

### Strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

### Check yer variables (probably redundant with the above but good practice)
# if [ -z ${var+x} ]; then echo "var is unset"; else echo "var is set to '$var'"; fi

### Load Modules
module load porechop/github-HEAD 
module load seqtk/1.2-IGB-gcc-4.9.4 

### Run app on file

NAME=${LR%.fastq.gz}

porechop-runner.py --discard_middle \
        -i $LR --threads $SLURM_NPROCS 2> $NAME.porechop.log | \
        seqtk seq -L 800 | gzip - > $NAME.qualtrim.clean.fastq.gz
