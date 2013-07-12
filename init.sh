#!/usr/bin/env bash

controldir=control
#==================================
# Necessary environment variables.
# Feel free to modify them if you know what you are doing.



# Uncomment the two lines below to remove RNA folding features and SilVA's
# dependency on UNAfold/ViennaRNA:
#export EXCLUDE_RNA_FOLDING=1
#controldir=control/no-folding
if [[ ! -z "${EXCLUDE_RNA_FOLDING:-}" ]]; then
    echo "EXCLUDE_RNA_FOLDING=${EXCLUDE_RNA_FOLDING}... excluding RNA folding features." >&2
fi


# Maximum number of threads to use
export SILVA_N_THREADS="${SILVA_N_THREADS:-8}"
# Path to the root of the SilVA directory.
export SILVA_PATH="${SILVA_PATH:-$(cd -P $(dirname $0); pwd)}"
# Control dataset to use
export SILVA_CONTROL="${SILVA_CONTROL:-$SILVA_PATH/$controldir/NA10851}"
# Directory of pre-trained models
export SILVA_TRAINED="${SILVA_TRAINED:-$SILVA_PATH/$controldir/models}"

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
COMMAND:         '$@'
SILVA_N_THREADS: '$SILVA_N_THREADS'
SILVA_CONTROL:   '$SILVA_CONTROL'
SILVA_TRAINED:   '$SILVA_TRAINED'
SILVA_AF_MIN:    '$SILVA_AF_MIN'
SILVA_AF_MAX:    '$SILVA_AF_MAX'
TMPDIR:          '$TMPDIR'
-----------
EOF
}

if [[ ! -e $SILVA_PATH/data/refGene.ucsc.gz ]]; then
    cat >&2 <<EOF 
Annotation databases seem to be missing.
File does note exist: $SILVA_PATH/data/refGene.ucsc.gz

Please run the setup.sh script in this tool's root directory
and make sure all necessary data files are in $SILVA_PATH/data/
EOF
    exit 1
fi
