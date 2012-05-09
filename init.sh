#!/usr/bin/env bash

set -eu
set -o pipefail


#==================================
# Necessary environment variables.
# Feel free to modify them if you think you know what you are doing.

export SYNORDER_PATH="${SYNORDER_PATH:-$(dirname $0)}"
export TMPDIR="${TMPDIR:-.}"
export CLASSPATH="${SYNORDER_PATH}/tools/weka/weka.jar:${CLASSPATH:-}"
python_version="$(python -c "import sys; print sys.version[:3]")"
export PYTHONPATH="${SYNORDER_PATH}/lib/python${python_version}:${PYTHONPATH:-}"

export UNAFOLD_BIN="${SYNORDER_PATH}/tools/unafold/src/hybrid-ss-min"
export UNAFOLD_DATA="${SYNORDER_PATH}/tools/unafold/data"
#==================================

if [[ ! -e $UNAFOLD_BIN ]]; then
    echo >&2 <<EOF 
UNAfold does not appear to be installed.
File does note exist: $UNAFOLD_BIN

Please run the setup.sh script in this tool's root directory.
EOF
    exit 1
elif [[ ! -d $UNAFOLD_DATA ]]; then
    echo >&2 <<EOF 
UNAfold does not appear to be installed.
Directory does note exist: $UNAFOLD_DATA

Please run the setup.sh script in this tool's root directory.
EOF
    exit 1
elif [[ ! -e $SYNORDER_PATH/data/refGene.pkl.gz ]]; then
    echo >&2 <<EOF 
Annotation databases seem to be missing.
File does note exist: $SYNORDER_PATH/data/refGene.pkl.gz

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
