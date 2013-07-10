#!/usr/bin/env bash

set -eu
set -o pipefail

version="$(cat VERSION)"
SILVA_DATA_URL="http://compbio.cs.toronto.edu/silva/release/silva-${version}_data.tar.gz"
SILVA_DATA=${SILVA_DATA:-data}

UNAFOLD_VERSION="3.8"
UNAFOLD_URL="http://mfold.rna.albany.edu/cgi-bin/UNAFold-download.cgi?unafold-${UNAFOLD_VERSION}.tar.gz"

VIENNA_VERSION="2.1.1"
VIENNA_URL="http://www.tbi.univie.ac.at/~ronny/RNA/packages/source/ViennaRNA-${VIENNA_VERSION}.tar.gz"

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


function get_and_make_package {
    local URL=$1
    local from=$2
    local to=$3
    
    pushd tools > /dev/null
    if [[ ! -e $from.tar.gz ]]; then
	prompt "\n\nDownloading $from..."
        wget -O $from.tar.gz "$URL"
    fi
    if [[ ! -d $to ]]; then
	prompt "\n\nUnpacking $from.tar.gz to tools/$to..."
        tar -xzvf $from.tar.gz
        mv $from $to
    fi
    pushd $to > /dev/null
    prompt "\n\nConfiguring and making $to..."
    ./configure && make
    popd > /dev/null
    popd > /dev/null
}

# Download unafold
if [[ ! -e tools/unafold/src/hybrid-ss-min ]]; then
    get_and_make_package $UNAFOLD_URL unafold-${UNAFOLD_VERSION} unafold
fi

# Make vienna
if [[ ! -e tools/vienna/install/bin/RNAfold ]]; then
    get_and_make_package $VIENNA_URL ViennaRNA-${VIENNA_VERSION} vienna
fi

# Install R randomForest package
function is_rf_missing {
    echo 'suppressMessages(require("randomForest", quietly=TRUE))' \
        | R --vanilla --quiet --slave 2>&1
}

function rf_fail {
    cat <<EOF

Problem encountered trying to automatically install the randomForest module
for R. Please make sure R is installed, in your path, and you have installed
the randomForest module before running SilVA. If you start R, randomForest
can be installed by typing install.packages("randomForest").

EOF
    exit 1
}


# Install randomForest R package
if [[ "$(is_rf_missing)" ]]; then
    prompt "Installing randomForest package for R..."
    echo 'install.packages("randomForest", repos="http://probability.ca/cran/")' \
	| R --vanilla --quiet || rf_fail
fi
if [[ "$(is_rf_missing)" ]]; then
    rf_fail
fi


# Download data
if [[ ! -e $SILVA_DATA/refGene.ucsc.gz ]]; then
    datafile=silva-${version}_data.tar.gz

    if [[ ! -e $datafile ]]; then
	prompt "\n\nFinal step:\nDownloading required SilVA $version databases ($datafile)...\nWARNING: these databases are rather large (~1GB), so this might take a while...\nIf you already downloaded and unpacked this file, exit this script and set the SILVA_DATA\nenvironment variable to path of the unpacked data directory."
	wget -v "$SILVA_DATA_URL"
    fi
    if [[ ! -e $SILVA_DATA/refGene.ucsc.gz ]]; then
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

To run SilVA on a VCF file, first preprocess the file with 
'./silva-preprocess', and then generate results from the output 
directory with './silva-run'
EOF

