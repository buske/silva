#!/usr/bin/env bash

set -eu
set -o pipefail

version="$(cat VERSION)"
SILVA_DATA_URL="http://compbio.cs.toronto.edu/silva/release/silva-${version}_data.tar.gz"
SILVA_DATA=${SILVA_DATA:-data}

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

echo -e "\nInstalling SilVA $version dependencies..." >&2

# Install randomForest R package


# Download data
if [[ ! -e $SILVA_DATA/refGene.pkl ]]; then
    datafile=silva-${version}_data.tar.gz

    if [[ ! -e $datafile ]]; then
	prompt "\n\nFinal step:\nDownloading required SilVA $version databases ($datafile)...\nWARNING: these databases are rather large (~700MB), so this might take a while...\nIf you already downloaded this file, exit this script and set the SILVA_DATA\nenvironment variable to data directory's path."
	wget -v "$SILVA_DATA_URL"
    fi
    if [[ ! -e $SILVA_DATA/refGene.pkl ]]; then
	echo -e "\n\nUnpacking required SilVA $version databases..." >&2
	if [[ -e $datafile ]]; then
	    tar -xzvf $datafile
	else
	    echo -e "\n\nError: missing file: $datafile"
	    exit 1
	fi
    fi
fi



cat >&2 <<EOF


SilVA $version successfully installed.

To run SilVA on a VCF file, first preprocess the file with './silva-preprocess',
and then generate results from the output directory with './silva-run'
EOF

