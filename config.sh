#!/usr/bin/env bash

set -eu
set -o pipefail

#==================================
# Necessary environment variables.
# Feel free to modify them if you know what you are doing.

# Path to the root of the Synorder directory.
export SYNORDER_PATH="${SYNORDER_PATH:-$(readlink -e $(dirname $0))}"

# Directory to use for any temporary files. It's recommended that this point
# to somewhere on the local machine (such as /tmp).
export TMPDIR="${TMPDIR:-$(pwd)}"

# Synorder ignores variants with a 1000 Genomes Project allele frequency
# greater than this number ([0--1], suggested: <=0.05)
export SYNORDER_AF_THRESH=0.05

# Python and Java library paths
export CLASSPATH="${SYNORDER_PATH}/tools/weka/weka.jar:${CLASSPATH:-}"
python_version="$(python -c "import sys; print sys.version[:3]")"
export PYTHONPATH="${SYNORDER_PATH}/lib/python${python_version}:${PYTHONPATH:-}"

# UNAfold paths
export UNAFOLD_BIN="${UNAFOLD_BIN:-${SYNORDER_PATH}/tools/unafold/src/hybrid-ss-min}"
export UNAFOLDDAT="${UNAFOLDDAT:-${SYNORDER_PATH}/tools/unafold/data}"
#==================================

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
elif [[ ! -e $SYNORDER_PATH/data/refGene.pkl ]]; then
    echo >&2 <<EOF 
Annotation databases seem to be missing.
File does note exist: $SYNORDER_PATH/data/refGene.pkl

Please run the setup.sh script in this tool's root directory.
EOF
    exit 1
elif [[ $(ls -1d $SYNORDER_PATH/lib/python${python_version}/milk-*.egg 2> /dev/null | wc -l) -eq 0 ]]; then
    echo >&2 <<EOF 
Necessary Python package has not been installed.
Directory does note exist: $SYNORDER_PATH/lib/python${python_version}/milk-*.egg

Please run the setup.sh script in this tool's root directory.
EOF
    exit 1

fi
