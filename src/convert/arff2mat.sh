#!/usr/bin/env bash


# Moves class to first column
awk 'BEGIN {IGNORECASE=1;} \
/^@attribute/ {attr[a++] = $2;} \
/^@data/ {FS=","; printf "#%s",attr[a-1]; for (i=0; i<(a-1); i++) {printf "\t%s",attr[i];} print "";} \
!/^@/ && !/^$/ {printf "%s",$NF; for (i=1; i<NF; i++){printf "\t%s",$i;} print "";}' "$@"
