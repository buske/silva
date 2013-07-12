#!/usr/bin/env bash


set -eu
set -o pipefail


######################################################################
# Uncomment the line below to exclude RNA folding (UNAfold/ViennaRNA)
# from the installation process:

#EXCLUDE_RNA_FOLDING=1

######################################################################

version="$(cat VERSION)"
SILVA_DATA_URL="http://compbio.cs.toronto.edu/silva/release/silva-${version}_data.tar.gz"

GENOME_URL="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit"


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

if [ -z "${EXCLUDE_RNA_FOLDING:-}" ]; then
    # Download unafold
    if [[ ! -e tools/unafold/src/hybrid-ss-min ]]; then
	get_and_make_package $UNAFOLD_URL unafold-${UNAFOLD_VERSION} unafold
    fi

    # Make vienna
    if [[ ! -e tools/vienna/Progs/RNAfold ]]; then
	get_and_make_package $VIENNA_URL ViennaRNA-${VIENNA_VERSION} vienna
    fi
else
    echo "EXCLUDE_RNA_FOLDING=${EXCLUDE_RNA_FOLDING}... excluding RNA folding packages." >&2
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
datafile=data/refGene.ucsc.gz
if [[ ! -e $datafile ]]; then
    tarball=silva-${version}_data.tar.gz
    if [[ ! -e $tarball ]]; then
	prompt "\n\nDownloading required SilVA $version databases ($tarball)..."
	wget -v "$SILVA_DATA_URL"
    fi
    if [[ ! -e $datafile ]]; then
	echo -e "\n\nUnpacking required SilVA $version databases into data/..." >&2
	if [[ -e $tarball ]]; then
	    tar -xzvf $tarball
	else
	    echo -e "\n\nError: missing file: $tarball"
	    exit 1
	fi
    fi
fi

# Download hg19 twobit file
mkdir -pv data
datafile=data/hg19.2bit
if [[ ! -e $datafile ]]; then
    if [[ ! -e $datafile ]]; then
	prompt "\n\nFinal step:\nDownloading hg19 2-bit file from UCSC...\n\nThis file is large (~800M) and may take a while to download.\n\nIf you already have hg19.2bit, press CTRL+C to exit this script and make a symlink to it within the data directory. Then, re-run this script to make sure everything is properly linked."
	wget -v -O $datafile $GENOME_URL
    fi
    if [[ ! -e $datafile ]]; then
	echo -e "\n\nError: missing file: $datafile"
	exit 1
    fi
fi


cat >&2 <<EOF


SilVA $version successfully installed.

To run SilVA on a VCF file, first preprocess the file with 
'./silva-preprocess', and then generate results from the output 
directory with './silva-run'
EOF

