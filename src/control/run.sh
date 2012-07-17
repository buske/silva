
set -eu
set -o pipefail

pcoord=$1
base=$(dirname $1)/$(basename $1 .pcoord)

mrna=$base.mrna
mat=$base.mat
if [[ ! -e $mrna ]]; then
    cat $pcoord  \
	| ./src/input/synonymous.py --protein-coords -O data/refGene.pkl filter - \
	| ./src/input/synonymous.py -O data/refGene.pkl annotate - \
	> $mrna
fi

if [[ ! -e $mat ]]; then
    ./src/convert/mrna2mat.sh $mrna $base
fi

#./src/convert/mat2arff.sh polymorphic 0 /data/buske/synonymous/alamut/polymorphic{.mat,}

#./src/control/preprocess.sh /data/buske/synonymous/2012-06-01_no-splice-dist/{true/true.arff,1000gp/NA10851.merged_beagle_mach.20101123.exonic-snps/NA10851.merged_beagle_mach.20101123.exonic-snps.raw.arff,NA10851}

# Update control symlinks
