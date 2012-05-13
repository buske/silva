#!/usr/bin/env bash

set -eu
set -o pipefail

function usage {
    cat <<EOF
Usage: $0 SAMPLE.arff OUTDIR

Requires: ./true.arff, and ./true.groups
Creates: SAMPLE.merged.arff, SAMPLE_LOO/,
         SAMPLE_5050/, in OUTDIR/

Existing files are preserved.
EOF
    exit 1
}

if [[ $# -ne 2 ]]; then
    usage
fi
samplefile=$1
datasetdir=$2
sample=$(basename $samplefile .arff)

if [[ ! -e true.arff || ! -e $samplefile ]]; then
    echo -e "Error: missing input files\n" >&2
    usage
fi

mkdir -pv $datasetdir
if [[ -e $datasetdir/$sample.merged.arff ]]; then
    echo "Found: $datasetdir/$sample.merged.arff" >&2
else
    echo "Creating merged, standardized data file: $datasetdir/$sample.merged.arff..." >&2
    $SILVA_PATH/util/standardize.sh true.arff $samplefile > $datasetdir/$sample.merged.arff
fi

echo "Creating 50/50 stratified train/test datasets: $datasetdir/${sample}_5050/ALL..." >&2
$SILVA_PATH/bench/split_data.py -N 25 $datasetdir/$sample.merged.arff $datasetdir/${sample}_5050/ALL

echo "Creating pseudo-LOO train/test datasets: $datasetdir/${sample}_LOO/ALL..." >&2
$SILVA_PATH/bench/split_data.py --infection --true-groups=true.groups \
    $datasetdir/$sample.merged.arff $datasetdir/${sample}_LOO/ALL

function make_subsets {
    # Files in indir must be of the form: *.(train|test).arff
    local indir=$1
    local outdir=$2
    #local outdir2=$3
    mkdir -pv $outdir #$outdir2
    echo "Creating feature-selected train/test datasets: $outdir..." >&2
    for train in $indir/*.train.arff; do 
	train=$(basename $train)
	test=$(basename $train .train.arff).test.arff
	$SILVA_PATH/util/subset_attributes.sh $indir/$train $indir/$test $outdir/$train $outdir/$test
	#$SILVA_PATH/util/subset_attributes2.sh $indir/$train $indir/$test $outdir2/$train $outdir2/$test
    done
}

make_subsets $datasetdir/${sample}_5050/ALL
#make_subsets $datasetdir/${sample}_5050/FSUB
make_subsets $datasetdir/${sample}_LOO/ALL
#make_subsets $datasetdir/${sample}_LOO/FSUB

for dir in $datasetdir/${sample}_LOO/ALL $datasetdir/${sample}_5050/ALL; do
    rankdir=$dir/scored
    echo "Running models to prioritize variants: $rankdir..." >&2
    $SILVA_PATH/models/run.sh $dir $rankdir
done
