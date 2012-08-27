#!/usr/bin/env bash

#==================================
# Necessary environment variables.
# Feel free to modify them if you know what you are doing.

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

# Python library paths
pyversion=$(python -c "import sys; print sys.version[:3]")
export PYTHONPATH="${SILVA_PATH}/lib/python$pyversion:${PYTHONPATH:-}"

# UNAfold paths
export UNAFOLD_BIN="${UNAFOLD_BIN:-${SILVA_PATH}/tools/unafold/src/hybrid-ss-min}"
export UNAFOLDDAT="${UNAFOLDDAT:-${SILVA_PATH}/tools/unafold/data}"

#==================================

version="$(cat ${SILVA_PATH}/VERSION)"

function init_message {
    cat >&2 <<EOF

SILVA $version
-----------
COMMAND:         '$0 $@'
DATA:            '$SILVA_DATA'
SILVA_CONTROL:   '$SILVA_CONTROL'
SILVA_TRAINED:   '$SILVA_TRAINED'
TMPDIR:          '$TMPDIR'
EOF
}

if [[ ! -e $UNAFOLD_BIN ]]; then
    cat >&2 <<EOF 
UNAfold does not appear to be installed.
File does note exist: $UNAFOLD_BIN

Please run the setup.sh script in this tool's root directory
or set the UNAFOLD_BIN environment variable appropriately.
EOF
    exit 1
elif [[ ! -d $UNAFOLDDAT ]]; then
    cat >&2 <<EOF 
UNAfold does not appear to be installed.
Directory does note exist: $UNAFOLDDAT

Please run the setup.sh script in this tool's root directory
or set the UNAFOLDDAT environment variable appropriately.
EOF
    exit 1
elif [[ ! -e $SILVA_DATA/refGene.pkl ]]; then
    cat >&2 <<EOF 
Annotation databases seem to be missing.
File does note exist: $SILVA_DATA/refGene.pkl

Please run the setup.sh script in this tool's root directory
or set the SILVA_DATA environment variable appropriately.
EOF
    exit 1
else
    milk_version=$(python -c "import milk; print milk.__version__")
    if [[ -z $milk_version ]]; then
	cat >&2 <<EOF 
Custom version of Python package, milk, does not appear to be properly installed.

Please run the setup.sh script in this tool's root directory
or set the PYTHON_PATH environment variable appropriately.
EOF
	exit 1
    fi
    if [[ $milk_version != silva-$version ]]; then
	cat >&2 <<EOF 
Found existing milk version: $milk_version.
This other version appears to take precedence over the custom version
included in this package.

Please resolve this, such that:
$ python -c 'import milk; print milk.__version__' 
prints "silva-$version"

The following may be sufficient:
$ cd tools/milk
$ python setup.py install
EOF
	exit 1
    fi
fi
