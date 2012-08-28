
set -eu
set -o pipefail

vcf=$1
base=$(dirname $1)/$(basename $1 .vcf)

$SILVA_PATH/silva-preprocess $vcf $(dirname $vcf)
./src/input/standardize.py --class=0 $base.mat \
    > $base.input

