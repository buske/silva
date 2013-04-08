#!/usr/bin/env bash

set -eu
set -o pipefail

MODELS_DIR=$SILVA_PATH/src/models

function usage {
    cat <<EOF
Usage: $0 DATADIR OUTDIR [MODEL...]

All ID.train.input files in DATADIR are run with their
matching test files.

MODEL should be a subfolder of $MODELS_DIR/ with an appropriate
'run' script.
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
    models=( fld forest gerp nnet svmmap )
fi

logdir=$outdir/log
mkdir -pv $logdir
for model in "${models[@]}"; do 
    if [[ $model == "nnet" ]]; then
	$SILVA_PATH/src/benchmark/run_model_on_dir.sh $MODELS_DIR/$model $datadir $outdir/$model $logdir
    else
	qsub -cwd -b y -V -e $logdir -o $logdir \
	    -l h_vmem=14G -N $model -R y -l h_rt="00:08:00" \
	    "bash -x $SILVA_PATH/src/benchmark/run_model_on_dir.sh $MODELS_DIR/$model $datadir $outdir/$model"
    fi
done