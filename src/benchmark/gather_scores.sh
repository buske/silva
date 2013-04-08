#!/usr/bin/env bash

set -u


dir=manuscript/benchmark
mkdir -p $dir/results

function score_one {
  local d=$1
  local name=$2

  for a in fld forest gerp nnet svmmap; do
    echo "$a"
    (echo $a; ./src/benchmark/score.py $dir/$d/scored/$a | tee >(cut -f 1 > $dir/results/.$name.rows) | grep -v "^#" | cut -f 2) > $dir/results/.$name.$a
  done
  paste $dir/results/.$name.rows $dir/results/.$name.{fld,forest,gerp,nnet,svmmap} > $dir/results/$name
}

# Method comparison split-wise
score_one NA10851_SPLIT "split"
# Method comparison infection
score_one NA10851_LOO "loo"

# Feature comparison for forest
a=forest
for f in codon_usage folding gerp motif sequence splice; do
  echo $f
  (echo "no_$f"; ./src/benchmark/score.py $dir/NA10851_${f}_FS_SPLIT/scored/$a | tee >(cut -f 1 > $dir/results/.fs.rows) | grep -v "^#" | cut -f 2) > $dir/results/.fs.$a.$f
done
paste $dir/results/.fs.rows $dir/results/.fs.$a.{codon_usage,folding,gerp,motif,sequence,splice} $dir/results/.split.gerp $dir/results/.split.$a > $dir/results/fs.$a

# Forest performance variance
#f=folding; ./src/benchmark/score.py --all $dir/NA10851_${f}_FS_LOO/scored/$a > $dir/results/loo.$f.$a.all
./src/benchmark/score.py --all $dir/NA10851_LOO/scored/$a > $dir/results/loo.$a.all
