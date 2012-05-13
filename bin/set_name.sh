#!/usr/bin/env bash

set -eu
set -o pipefail

#hg mv src/share/{sympri,silva}.py
find HISTORY LICENSE README VERSION run preprocess_vcf *.* sge/ src/ bin/ship.sh bin/build.sh tools/milk/milk/milk_version.py -type f \
	| grep -v "~$" | egrep -v "\.(pyc|pkl|gz|txt)$" \
	| xargs grep -li sympri \
	| xargs sed -e 's/SYMPRI/SILVA/g' -e 's/Sympri/Silva/g' -e 's/sympri/silva/g' -i""