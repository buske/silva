#!/usr/bin/env bash

set -eu
set -o pipefail

version="$(cat VERSION)"

function prompt {
    echo -e "$@" >&2
    read -p "(Press ENTER to continue, Ctrl-C to exit) "
}

echo "Installing Sympri $version dependencies..." >&2

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
    prompt "\nDownloading sympri $version databases..."
    if [[ ! -e sympri_${version}_data.tar.gz ]]; then
	wget http://compbio.cs.toronto.edu/sympri/sympri_${version}_data.tar.gz
    fi
    if [[ ! -d data ]]; then
	tar -xzf sympri_${version}_data.tar.gz
    fi
fi


function milfail {
    cat <<EOF
There seems to have been a problem installing the custom version of the 
Python package milk. If you have a .pydistutils.cfg file, this could be
the cause. Try moving that file out of the way, rerunning this script,
and then moving it back.

Sympri $version installation failed.
EOF
    exit 1
}

# Install modified version of milk
prefix="$(pwd)"
libdir="$prefix/lib/python"
if [[ $(python -c "import milk; print milk.__version__") != sympri-$version ]]; then
    prompt "\nConfiguring modified milk Python package..."
    mkdir -pv "$libdir"
    pushd tools/milk
    export PYTHONPATH="$libdir:$PYTHONPATH"
    trap milkfail TERM EXIT
    python setup.py install --home="$prefix"
    trap - TERM EXIT
    popd
fi

echo -e "\nSympri $version successfully installed." >&2
