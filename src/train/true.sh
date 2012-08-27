
set -eu
set -o pipefail

pcoord=$1
controlbase=$2
outdir=$3
basename=$(basename $pcoord .pcoord)
base=$(dirname $pcoord)/$basename


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
	--class=1 --control=$controlbase.mat $base.mat > $base.input
fi

mkdir -pv $outdir
merged=$outdir/$basename.merged.input

cp $base.input $merged
grep -v "^#" $controlbase.input >> $merged

./src/train/split_data.py $merged $outdir
./src/train/train_models.sh $outdir{,/models} forest

# Update control symlinks
