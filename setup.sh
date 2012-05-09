#!/usr/bin/env bash

set -eu
set -o pipefail

python_version=$(python -c "import sys; print sys.version[:3]")

function prompt {
    read -p "(Press ENTER to continue, Ctrl-C to exit) "
}

# Download unafold
if [[ ! -e tools/unafold/src/hybrid-ss-min ]]; then
    echo "Downloading UNAFold 3.8..." >&2
    prompt
    pushd tools
    if [[ ! -e unafold-3.8.tar.gz ]]; then
	wget http://mfold.rna.albany.edu/cgi-bin/UNAFold-download.cgi?unafold-3.8.tar.gz
    fi
    if [[ ! -d unafold ]]; then
	tar -xzf unafold-3.8.tar.gz
	mv unafold-3.8 unafold
    fi
    pushd unafold
    ./configure && make
    popd
    popd
fi


# Download data
if [[ ! -e data/refGene.pkl.gz ]]; then
    echo "Downloading synorder $version databases..." >&2
    prompt
    if [[ ! -e synorder_${version}_data.tar.gz ]]; then
	wget http://compbio.cs.toronto.edu/synorder/synorder_${version}_data.tar.gz
    fi
    if [[ ! -d data ]]; then
	tar -xzf synorder_${version}_data.tar.gz
    fi
fi

# Install modified version of milk
prefix="$(pwd)"
libdir="$prefix/lib/python${python_version}"
eggdir="$(ls -1d "$libdir/"milk-*.egg)"
if [[ ! -d $eggdir ]]; then
    echo "Configuring modified milk Python package..." >&2
    prompt
    mkdir -pv "$libdir"
    pushd tools/milk
    export PYTHONPATH="$libdir:$PYTHONPATH"
    python setup.py install --prefix="$prefix"
    popd
fi