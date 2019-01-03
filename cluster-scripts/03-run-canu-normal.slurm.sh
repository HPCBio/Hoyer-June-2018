#!/bin/bash
#SBATCH -c 24 
#SBATCH --mem=72000
#SBATCH -p normal 
#SBATCH --mail-type=END,FAIL
#SBATCH -J canu

### Strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

### Check yer variables (probably redundant with the above but good practice)
# if [ -z ${var+x} ]; then echo "var is unset"; else echo "var is set to '$var'"; fi

### Load Modules
module load Canu/1.7-IGB-gcc-4.9.4-Perl-5.24.1
module load Java/1.8.0_152
module load gnuplot/5.0.6-IGB-gcc-4.9.4

### Run app on file
canu -p asm -d $SAMPLE_NAME \
    genomeSize=$GENOME_SIZE \
    useGrid=false \
    -nanopore-raw $LR 
