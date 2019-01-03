#!/bin/bash
#SBATCH -n 24 
#SBATCH --mem=72000
#SBATCH -p normal 
#SBATCH --mail-type=END,FAIL
#SBATCH -J iterative-pilon

### Strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

### Check yer variables (probably redundant with the above but good practice)
# if [ -z ${var+x} ]; then echo "var is unset"; else echo "var is set to '$var'"; fi

### Load Modules
module load Unicycler/0.4.4-IGB-gcc-4.9.4-Python-3.6.1
module load Biopython/1.68-IGB-gcc-4.9.4-Python-3.6.1  

### Run app on file
python iterative-pilon.py \
    -i $ASSEMBLY \
    -o $SAMPLE_NAME.pilon.fasta \
    -1 $FQ_R1 \
    -2 $FQ_R2 \
    -a 200 \
    -d 800 \
    -t $SLURM_NPROCS \
    --diploid \
    --pilon_path /home/apps/software/pilon/1.22-Java-1.8.0_121/pilon-1.22.jar > $SAMPLE_NAME.it-pilon.log 2>&1
