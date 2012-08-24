
set -eu
set -o pipefail

pcoord=$1
control=$2
base=$(dirname $1)/$(basename $1 .pcoord)

flt=$base.flt
mrna=$base.mrna
mat=$base.mat

if [[ ! -s $flt ]]; then
    ./src/input/synonymous.py --protein-coords -O data/refGene.pkl filter $pcoord > $flt
fi

if [[ ! -s $mrna ]]; then
    ./src/input/synonymous.py -O data/refGene.pkl annotate $flt > $mrna
fi

if [[ ! -s $mat ]]; then
    ./src/convert/mrna2mat.sh $mrna $base
fi

if [[ ! -s $base.input ]]; then
    ./src/input/standardize.py \
	--class=1 --control=$control $base.mat > $base.input
fi

#./src/convert/mat2arff.sh polymorphic 0 /data/buske/synonymous/alamut/polymorphic{.mat,}

#./src/control/preprocess.sh /data/buske/synonymous/2012-06-01_no-splice-dist/{true/true.arff,1000gp/NA10851.merged_beagle_mach.20101123.exonic-snps/NA10851.merged_beagle_mach.20101123.exonic-snps.raw.arff,NA10851}

# Update control symlinks
