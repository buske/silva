#!/usr/bin/env bash

#==================================
# Necessary environment variables.
# Feel free to modify them if you know what you are doing.

# Maximum number of threads to use
export SILVA_N_THREADS=6
# Path to the root of the SilVA directory.
export SILVA_PATH="${SILVA_PATH:-$(cd -P $(dirname $0); pwd)}"
# Path of the untar'd SilVA data directory
export SILVA_DATA="${SILVA_DATA:-$SILVA_PATH/data}"
# Control dataset to use
export SILVA_CONTROL="${SILVA_CONTROL:-$SILVA_PATH/control/NA10851}"
# Directory of pre-trained models
export SILVA_TRAINED="${SILVA_TRAINED:-$SILVA_PATH/control/models}"

# Directory to use for any temporary files. It's recommended that this point
# to somewhere on the local machine (such as /tmp).
export TMPDIR="${TMPDIR:-$(pwd)}"

# SilVA ignores variants with 1000 Genomes Project allele frequencies
# less than SILVA_AF_MIN or greater than SILVA_AF_MAX.
# Both should be in [0--1], suggested: MIN: 0 and MAX: 0.05, respectively
export SILVA_AF_MIN="${SILVA_AF_MIN:-0}"
export SILVA_AF_MAX="${SILVA_AF_MAX:-0.05}"

#==================================

version="$(cat ${SILVA_PATH}/VERSION)"

function init_message {
    cat >&2 <<EOF

SILVA $version
-----------
COMMAND:         '$0 $@'
N_THREADS:       '$SILVA_N_THREADS'
DATA:            '$SILVA_DATA'
SILVA_CONTROL:   '$SILVA_CONTROL'
SILVA_TRAINED:   '$SILVA_TRAINED'
TMPDIR:          '$TMPDIR'
EOF
}

if [[ ! -e $SILVA_DATA/refGene.pkl ]]; then
    cat >&2 <<EOF 
Annotation databases seem to be missing.
File does note exist: $SILVA_DATA/refGene.pkl

Please run the setup.sh script in this tool's root directory
or set the SILVA_DATA environment variable appropriately.
EOF
    exit 1
fi
