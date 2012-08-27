#!/usr/bin/env bash

set -eu
set -o pipefail


# ========= CUSTOMIZABLE PARAMETERS ========= #
MODELROOT=$SILVA_PATH/src/models
n_threads=8
# =========================================== #

function usage {
    cat <<EOF
Usage: $0 MODEL DATADIR MODELDIR [LOGDIR]

Models are trained on all ID.train.input files in DATADIR
Models are written to MODELDIR/ID.model

$MODELROOT/MODEL/train should exist
If LOGDIR is set, then sub-jobs are submitted to SGE.
EOF
    exit 1
}

if [[ $# -eq 3 ]]; then
    logdir=""
elif [[ $# -eq 4 ]]; then
    logdir=$(readlink -f $4)
else
    usage
fi

echo "Command: $0 $@" >&2

model=$1
datadir=$(readlink -e $2)
outdir=$(readlink -f $3)

if [[ ! -e $MODELROOT/$model/train ]]; then
    echo "Could not find: $MODELROOT/$model/train" >&2
    exit 1
fi

mkdir -pv $outdir
function train_one {
    # Usage: $0 MODEL TRAINFILE MODELFILE
    local model=$1
    local trainfile=$2
    local modelfile=$3

    pushd $MODELROOT/$model > /dev/null
    id=$(basename $trainfile .train.input)
    if [[ -n $logdir ]]; then
	mkdir -pv $logdir
	qsub -cwd -e $logdir -o $logdir -b y -V \
	    -N "$model.$id" \
	    -l h_vmem=2G -q lunchQ \
	    "./train $modelfile $trainfile"
    else
	echo "$id" >&2
	./train $modelfile $trainfile
	i=$(expr $i + 1)
    fi
    popd > /dev/null
}

i=0
for trainfile in $datadir/*.train.input; do
    # Convert to mat
    if [[ ! -s $trainfile ]]; then
	continue
    fi
    id=$(basename $trainfile .train.input)
    modelfile=$outdir/$id.model
    if [[ ! -e $modelfile ]]; then
	train_one $model $trainfile $modelfile
    fi
    if [[ $i -ge $n_threads ]]; then
	wait
	i=0
    fi
done
wait
