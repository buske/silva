#!/usr/bin/env bash

set -eu
set -o pipefail

#hg mv src/share/{synorder,sympri}.py
find HISTORY LICENSE README VERSION run preprocess_vcf *.* sge/ src/ tools/milk/milk/milk_version.py -type f \
	| grep -v "~$" | egrep -v "\.(pyc|pkl|gz|txt)$" \
	| xargs grep -li synorder \
	| xargs sed -e 's/SYNORDER/SYMPRI/g' -e 's/Synorder/Sympri/g' -e 's/synorder/sympri/g' -i""