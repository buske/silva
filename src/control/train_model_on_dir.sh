#!/usr/bin/env bash

set -eu
set -o pipefail


# ========= CUSTOMIZABLE PARAMETERS ========= #
ARFF2MAT=$SILVA_PATH/src/convert/arff2mat.sh
MODELROOT=$SILVA_PATH/src/models
n_threads=8
# =========================================== #

function usage {
    cat <<EOF
Usage: $0 MODEL DATADIR MODELDIR [LOGDIR]

Models are trained on all ID.train.arff files in DATADIR
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

if [[ ! -e $ARFF2MAT ]]; then
    echo "Could not find: $ARFF2MAT" >&2
    exit 1
elif [[ ! -e $MODELROOT/$model/train ]]; then
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
    id=$(basename $trainfile .train.arff.mat)
    if [[ -n $logdir ]]; then
	mkdir -pv $logdir
	qsub -cwd -e $logdir -o $logdir -b y -V \
	    -N "$model.$id" \
	    -l h_vmem=2G -q lunchQ \
	    "./train $modelfile $trainfile"
    else
	echo "$id" >&2
	$SILVA_PATH/train $modelfile $trainfile
	i=$(expr $i + 1)
    fi
    popd > /dev/null
}

i=0
for trainfile in $datadir/*.train.arff; do
    # Convert to mat
    if [[ ! -e $trainfile.mat ]]; then
	echo "Creating MAT versions of train file: $trainfile.mat..." >&2
	$ARFF2MAT $trainfile > $trainfile.mat
    fi
    id=$(basename $trainfile .train.arff)
    modelfile=$outdir/$id.model
    if [[ ! -e $modelfile ]]; then
	train_one $model $trainfile.mat $modelfile
    fi
    if [[ $i -ge $n_threads ]]; then
	wait
	i=0
    fi
done
wait
