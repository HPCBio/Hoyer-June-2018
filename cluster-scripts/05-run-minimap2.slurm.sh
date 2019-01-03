#!/bin/bash
#SBATCH -n 12
#SBATCH --mem=24000
#SBATCH -p normal 
#SBATCH --mail-type=END,FAIL
#SBATCH -J minimap2

### Strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'
set -x
### Check yer variables (probably redundant with the above but good practice)
# if [ -z ${var+x} ]; then echo "var is unset"; else echo "var is set to '$var'"; fi

### Load Modules
module load minimap/2.8-IGB-gcc-4.9.4
module load SAMtools/1.5-IGB-gcc-4.9.4

### Run app on file
minimap2 -d $ASSEMBLY.mmi $ASSEMBLY

minimap2 -ax map-ont \
    -t 	$SLURM_NPROCS \
    $ASSEMBLY.mmi $LR | \
    samtools view -bS - > $SAMPLE_NAME.bam 

samtools sort -@ $SLURM_NPROCS -o $SAMPLE_NAME.nanopolish.sorted.bam $SAMPLE_NAME.bam

#bwa mem -x ont2d \
#     -t ${task.cpus} \
#     $canuAsm ${lr} | \
#     samtools sort -o ${id}.nanopolish.sorted.bam -T reads.tmp -

samtools index $SAMPLE_NAME.nanopolish.sorted.bam 
