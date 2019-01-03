#!/bin/bash
#SBATCH -n 28
#SBATCH --mem=96000
#SBATCH -p himemory
#SBATCH -w compute-2-0
#SBATCH --mail-type=END,FAIL
#SBATCH -J nanopolish-assembly

### Strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'
set -x

### Check yer variables (probably redundant with the above but good practice)
# if [ -z ${var+x} ]; then echo "var is unset"; else echo "var is set to '$var'"; fi

### Load Modules
module load SAMtools/1.5-IGB-gcc-4.9.4
module load nanopolish/0.9.0-IGB-gcc-4.9.4
module load parallel/20170622-IGB-gcc-4.9.4

### Run app on file
export NANOPOLISH_HOME=`dirname $( which nanopolish )`

SCRATCH=/scratch/${USER}-albacore

ABS_LR=$( readlink -e $LR )
BN_LR=$( basename $LR )
ABS_BAM=$( readlink -e $BAM )
BN_BAM=$( basename $BAM )
ABS_ASM=$( readlink -e $ASSEMBLY )
BN_ASM=$( basename $ASSEMBLY )

cd $SCRATCH
if [ ! -e $BN_LR ]; then
    ln -s $ABS_LR .
fi

if [ ! -e $BN_BAM ]; then
    ln -s ${ABS_BAM}* .
fi

if [ ! -e $BN_ASM ]; then
    ln -s $ABS_ASM .
fi

python $NANOPOLISH_HOME/scripts/nanopolish_makerange.py $BN_ASM | \
    parallel \
    --results nanopolish.results \
    -P 14 \
    nanopolish variants --consensus polished.{1}.fa \
    -w {1} \
    -r $BN_LR \
    -b $BN_BAM -g $BN_ASM -t 2 \
    --min-candidate-frequency 0.1 > $SAMPLE_NAME.polishing.log 2>&1

python $NANOPOLISH_HOME/scripts/nanopolish_merge.py polished.*.fa > \
    $SAMPLE_NAME.nanopolish.fasta 2> $SAMPLE_NAME.merging.log
