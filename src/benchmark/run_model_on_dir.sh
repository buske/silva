#!/usr/bin/env bash

set -eu
set -o pipefail


n_threads=8

function usage {
    cat <<EOF
Usage: $0 MODELDIR DATADIR OUTDIR [LOGDIR]

All ID.train.input files in DATADIR are run with their matching test files.

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

echo "Command: $0 $@" >&2

modeldir=$(readlink -e $1)
datadir=$(readlink -e $2)
outdir=$(readlink -m $3)

mkdir -pv $outdir
function run_one_prefix {
    # Usage: run_one_prefix prefix
    local prefix=$1

    local ids=( $(cd $datadir; find $prefix[.-]* -maxdepth 0 | sed -e 's/\..*//' | sort -u) )

    complete=1
    local id
    for id in "${ids[@]}"; do
	if [[ ! -s $datadir/$id.test.input ]]; then
	    echo "Error: missing expected test file: $datadir/$id.test.input" >&2
	    exit 1
	fi

	if [[ ! -s $outdir/$id.scored ]]; then
	    complete=
	fi
    done

    if [[ $complete ]]; then
	return 0
    fi

    pushd $modeldir > /dev/null
    if [[ -n $logdir ]]; then
	mkdir -pv $logdir
	qsub -S /bin/bash -cwd -e $logdir -o $logdir -b y -V -N "$(basename $modeldir).$prefix" \
	    -l h_vmem=2G -q lunchQ \
	    "for id in ${ids[@]}; do echo \$id; train=$datadir/\$id.train.input; test=$datadir/\$id.test.input; model=$outdir/\$id.model; ./train \$model \$train && ./test \$model \$test > $outdir/.\$id.scored && mv $outdir/.\$id.scored $outdir/\$id.scored && rm -f \$model; done"
    else
	for id in "${ids[@]}"; do
	    echo "$id" >&2
	    local model=$outdir/$id.model
	    local train=$datadir/$id.train.input
	    local test=$datadir/$id.test.input

	    ./train $model $train && ./test $model $test > $outdir/.$id.scored && mv $outdir/.$id.scored $outdir/$id.scored && rm -f $model &
	    i=$(expr $i + 1)

	    if [[ $i -ge $n_threads ]]; then
		wait
		i=0
	    fi
	done
    fi
    popd > /dev/null
}

i=0
for prefix in $(cd $datadir; find *.train.input -maxdepth 0 -type f | sed -e 's/[.-].*//' | sort -u); do
    run_one_prefix $prefix
done
wait