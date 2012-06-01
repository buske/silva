#!/usr/bin/env bash

set -eu
set -o pipefail

dir="$@"




while read line; do
    grep -l "\\b$line\\b" $dir/*/*/RESULTS | awk -F"/" '{print $(NF - 2)}' | sort | uniq -c | sort -nr | perl -pe 's/\s*(\d+)\s+(.*)/\1:\2/' | xargs echo | sed -e 's/ /,/g' | paste - <(echo "$line")
done < <(awk '$1 > 0.1' $dir/*/*/RESULTS | cut -f 4-9 | sort -u) \
    | sort -t":" -k 1,1nr
