#!/usr/bin/env bash

set -eu
set -o pipefail

version="$(cat VERSION)"
python_version=$(python -c "import sys; print sys.version[:3]")

function prompt {
    echo -e "$@" >&2
    read -p "(Press ENTER to continue, Ctrl-C to exit) "
}

echo "Installing Synorder $version dependencies..." >&2

# Download unafold
if [[ ! -e tools/unafold/src/hybrid-ss-min ]]; then
    prompt "\nDownloading UNAFold 3.8..."
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
if [[ ! -e data/refGene.pkl ]]; then
    prompt "\nDownloading synorder $version databases..."
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
eggdir="$(ls -1d "$libdir/"milk-*.egg 2> /dev/null || true)"
if [[ ! -d $eggdir ]]; then
    prompt "\nConfiguring modified milk Python package..."
    mkdir -pv "$libdir"
    pushd tools/milk
    export PYTHONPATH="$libdir:$PYTHONPATH"
    python setup.py install --prefix="$prefix"
    popd
fi

echo -e "\nSynorder $version successfully installed." >&2
