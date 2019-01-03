#!/bin/bash
#SBATCH -n 2
#SBATCH --mem=96000
#SBATCH -p himemory
#SBATCH -w compute-2-0
#SBATCH --mail-type=END,FAIL
#SBATCH -J nanopolish-assembly

### Strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'
### Check yer variables (probably redundant with the above but good practice)
# if [ -z ${var+x} ]; then echo "var is unset"; else echo "var is set to '$var'"; fi

### Load Modules
module load SAMtools/1.5-IGB-gcc-4.9.4
module load nanopolish/0.9.0-IGB-gcc-4.9.4

### Run app on file
export NANOPOLISH_HOME=`dirname \$( which nanopolish )`

SCRATCH=/scratch/${USER}-albacore

ABS_LR=$( readlink -e $LR )
BN_LR=$( basename $LR )
ABS_SUM=$( readlink -e $SUMMARY )
BN_SUM=$( basename $SUMMARY )

cd $SCRATCH

if [ ! -e $BN_LR ]; then
    ln -s $ABS_LR .
fi
     
#nanopolish index -s $ABS_SUM -v -d FAST5 $BN_LR > $SAMPLE_NAME.index.log 2>&1
nanopolish index -v -d FAST5 -s $ABS_SUM $BN_LR > $SAMPLE_NAME.index.log 2>&1
