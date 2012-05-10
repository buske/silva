#!/usr/bin/env bash

set -eu
set -o pipefail

#==================================
# Necessary environment variables.
# Feel free to modify them if you know what you are doing.

# Path to the root of the Sympri directory.
export SYMPRI_PATH="${SYMPRI_PATH:-$(cd -P $(dirname $0); pwd)}"

# Directory to use for any temporary files. It's recommended that this point
# to somewhere on the local machine (such as /tmp).
export TMPDIR="${TMPDIR:-$(pwd)}"

# Sympri ignores variants with a 1000 Genomes Project allele frequency
# greater than of equal to this number ([0--1], suggested: 0.05)
export SYMPRI_AF_THRESH=0.05

# Python and Java library paths
export CLASSPATH="${SYMPRI_PATH}/tools/weka/weka.jar:${CLASSPATH:-}"
export PYTHONPATH="${SYMPRI_PATH}/lib/python:${PYTHONPATH:-}"

# UNAfold paths
export UNAFOLD_BIN="${UNAFOLD_BIN:-${SYMPRI_PATH}/tools/unafold/src/hybrid-ss-min}"
export UNAFOLDDAT="${UNAFOLDDAT:-${SYMPRI_PATH}/tools/unafold/data}"

# Control dataset to use
export SYMPRI_CONTROL=$SYMPRI_PATH/data/control/NA10851
#==================================

version="$(cat VERSION)"

function config_message {
    cat >&2 <<EOF
SYMPRI $version
------------

TMPDIR:           '$TMPDIR'
SYMPRI_AF_THRESH: '$SYMPRI_AF_THRESH'
SYMPRI_CONTROL:   '$SYMPRI_CONTROL'
EOF
}

if [[ ! -e $UNAFOLD_BIN ]]; then
    echo >&2 <<EOF 
UNAfold does not appear to be installed.
File does note exist: $UNAFOLD_BIN

Please run the setup.sh script in this tool's root directory.
EOF
    exit 1
elif [[ ! -d $UNAFOLDDAT ]]; then
    echo >&2 <<EOF 
UNAfold does not appear to be installed.
Directory does note exist: $UNAFOLDDAT

Please run the setup.sh script in this tool's root directory.
EOF
    exit 1
elif [[ ! -e $SYMPRI_PATH/data/refGene.pkl ]]; then
    echo >&2 <<EOF 
Annotation databases seem to be missing.
File does note exist: $SYMPRI_PATH/data/refGene.pkl

Please run the setup.sh script in this tool's root directory.
EOF
    exit 1
else
    milk_version=$(python -c "import milk; print milk.__version__")
    if [[ -z $milk_version ]]; then
	echo >&2 <<EOF 
Custom version of Python package, milk, does not appear to be properly installed.

Please run the setup.sh script in this tool's root directory.
EOF
	exit 1
    fi
    if [[ $milk_version != sympri-$version ]]; then
	echo >&2 <<EOF 
Found existing milk version: $milk_version.
This other version appears to take precedence over the custom version
included in this package.

Please resolve this, such that:
$ python -c 'import milk; print milk.__version__' 
prints out sympri-$version

The following may be sufficient:
$ cd tools/milk
$ python setup.py install
EOF
	exit 1
    fi
fi
