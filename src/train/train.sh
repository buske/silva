
set -eu
set -o pipefail

pcoord=$1
controlbase=$2
outdir=$(dirname $1)
base=$(basename $1 .pcoord)

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
    ./src/convert/mrna2mat.sh $mrna $outdir/$base
fi

if [[ ! -s $outdir/$base.input ]]; then
    ./src/input/standardize.py \
	--class=1 --control=$controlbase.mat $outdir/$base.mat > $outdir/$base.input
fi

cp $outdir/$base.input $outdir/$base.merged.input
grep -v "^#" $controlbase.input >> $outdir/$base.merged.input

./src/train/split_data.py $outdir/$base.merged.input $outdir/${base}_training
./src/train/train_models.sh $base.training{,/models} forest

# Update control symlinks
