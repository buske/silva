#!/usr/bin/env bash

set -eu
set -o pipefail


# ========= CUSTOMIZABLE PARAMETERS ========= #
ARFF2MAT=$SILVA_PATH/src/convert/arff2mat.sh
n_threads=8
# =========================================== #

function usage {
    cat <<EOF
Usage: $0 MODELDIR DATADIR OUTDIR [LOGDIR]

All ID.train.arff files in DATADIR are run with their matching test files.

MODELDIR should have 'train' and 'test' scripts in it.
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

if [[ ! -e $ARFF2MAT ]]; then
    echo "Could not find: $ARFF2MAT" >&2
    exit 1
fi

echo "Command: $0 $@" >&2

modeldir=$(readlink -e $1)
datadir=$(readlink -e $2)
outdir=$(readlink -f $3)

mkdir -pv $outdir
function run_one_example {
    # Usage: run_one_example ID TRAIN TEST
    local id=$1
    local train=$2
    local test=$3
    local model=$outdir/$id.model

    pushd $modeldir > /dev/null
    if [[ ! -s $outdir/$id.scored ]]; then
	if [[ -n $logdir ]]; then
	    mkdir -pv $logdir
	    qsub -cwd -e $logdir -o $logdir -b y -V -N "$(basename $modeldir).$id" \
		-l h_vmem=2G -q lunchQ \
		"./train $model $train && ./test $model $test > $outdir/.$id.scored && mv $outdir/.$id.scored $outdir/$id.scored && rm -f $model"
	else
	    echo "$id" >&2
	    ./train $model $train && ./test $model $test > $outdir/.$id.scored && mv $outdir/.$id.scored $outdir/$id.scored && rm -f $model &
	    i=$(expr $i + 1)
	fi
    fi
    popd > /dev/null
}

i=0
for train_file in $datadir/*.train*.arff; do
    test_file=$(echo $train_file | sed -e s/\.train/\.test/)
    
    id=${train_file##*/}
    id=${id%%.train*}
    
    if [[ ! -s $test_file ]]; then
	echo "Error: missing expected test file: $test_file" >&2
	exit 1
    fi
    # Convert to mat
    if [[ ! -s $train_file.mat || ! -s $test_file.mat ]]; then
	echo "Creating MAT versions of train/test files..." >&2
	$ARFF2MAT $train_file > $train_file.mat
	$ARFF2MAT $test_file > $test_file.mat
    fi
    
    run_one_example $id $train_file.mat $test_file.mat
    if [[ $i -ge $n_threads ]]; then
	wait
	i=0
    fi
done
wait