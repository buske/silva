#!/usr/bin/env bash

set -eu
set -o pipefail

version="$(cat VERSION)"
UNAFOLD_URL="http://mfold.rna.albany.edu/cgi-bin/UNAFold-download.cgi?unafold-3.8.tar.gz"
SYMPRI_DATA_URL="http://compbio.cs.toronto.edu/sympri/sympri-${version}_data.tar.gz"

# Check python verison
pyversion=$(python -c "import sys; print sys.version[:3]")
if [[ $pyversion != 2.6 && $pyversion != 2.7 ]]; then
    cat >&2 <<EOF
Found Python $pyversion, but requires Python 2.6 or 2.7.

Please install Python 2.6/2.7, or set your PATH such that the 'python'
command loads Python 2.6 or 2.7.
EOF
    exit 1
fi

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
	wget "$UNAFOLD_URL"
    fi
    if [[ ! -d unafold ]]; then
	tar -xzvf unafold-3.8.tar.gz
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
    datafile=sympri-${version}_data.tar.gz
    if [[ ! -e $datafile ]]; then
	#wget --progress "$SYMPRI_DATA_URL"
	ln -s /tmp/buske/sympri/dist/$datafile
    fi
    if [[ ! -d data ]]; then
	tar -xzvf $datafile
    fi
fi


function milkfail {
    cat <<EOF
There seems to have been a problem installing the custom version of the 
Python package milk. If you have a custom distutils configuration, this
could be the cause. Try moving that file out of the way, rerunning this 
script, and then moving it back.

Sympri $version installation failed.
EOF
    exit 1
}


# Install modified version of milk
prefix="$(pwd)"
if [[ $(python -c "import milk; print milk.__version__" 2> /dev/null) != sympri-$version ]]; then
    prompt "\nConfiguring modified milk Python package..."
    pushd tools/milk
    if [[ -e "$HOME/.pydistutils.cfg" ]]; then
	python setup.py install
    else
	libdir="$prefix/lib/python$pyversion"
	mkdir -pv "$libdir"
	export PYTHONPATH="$libdir:${PYTHONPATH:-}"
	trap milkfail TERM EXIT
	python setup.py install --install-lib="$libdir"
	
	trap - TERM EXIT
    fi
    popd
fi

echo -e "\nSympri $version successfully installed." >&2
