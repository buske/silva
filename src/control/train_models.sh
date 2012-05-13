#!/usr/bin/env bash

set -eu
set -o pipefail


function usage {
    cat <<EOF
Usage: $0 DATADIR MODELDIR [MODEL...]

Expects ID.train.arff files in DATADIR/
Creates ID.model in MODELDIR/MODEL/ for each specified.
EOF
    exit 1
}

if [[ $# -lt 2 ]]; then
    usage
fi

datadir=$(readlink -e $1)
outdir=$(readlink -f $2)
shift 2

echo "Command: $0 $@" >&2

if [[ $# -gt 0 ]]; then
    models=( "$@" )
else
    models=( fld forest nnet svmmap )
fi

logdir=$outdir/log
mkdir -pv $logdir
for model in "${models[@]}"; do 
    if [[ $model == "nnet" || $model == forest* ]]; then
	$SILVA_PATH/src/control/train_model_on_dir.sh $model $datadir $outdir/$model $logdir
    else
	qsub -cwd -b y -V -e $logdir -o $logdir \
	    -l h_vmem=14G -N $model -q lunchQ \
	    "bash -x $SILVA_PATH/src/control/train_model_on_dir.sh $model $datadir $outdir/$model"
    fi
done