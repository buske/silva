#!/usr/bin/env bash

awk '! $6 || $6 >= 30;' "$@" | grep -v "INDEL" | perl -ne 'print if (/^#/ || /DP4=(\d+),(\d+),(\d+),(\d+)/ && ($1 + $2 + $3 + $4 >= 3) && ($3 > 0) && ($4 > 0))'
