#!/usr/bin/env bash

set -eu
set -o pipefail

function prompt {
    echo -e "$@" >&2
    read -p "(ENTER to continue)... "
}

name="silva"
version=$(cat VERSION)
echo "Building release for version: $version" >&2

prefix=${name}-${version}
builddir=build/$prefix
distbase=dist/$prefix

mkdir -pv build dist



# Tarball source in two steps!
# 1) Copy relevant source files to build
prompt "\nCopying files to build directory: $builddir"
if [[ ! -e .rsync-filter ]]; then
    echo "Error: expected .rsync-filter file in cwd" >&2
    exit 1
fi
mkdir -pv $builddir
rsync -r \
    -v \
    -u \
    --copy-links \
    --filter='. .rsync-filter' \
    . \
    $builddir/


prompt "\nReplacing <VERSION> tag and coded versions with $version within $builddir"
pushd $builddir
if [[ $version != $(cat VERSION) ]]; then
    echo "Error: version mismatch" >&2
    exit 1
fi
replace "<VERSION>" $version -- README
popd

# 2) tarball relevant source files
src=${distbase}.tar.gz
cont=y
if [[ -e $src ]]; then
    read -p "Overwrite $src (y/n)? " cont
fi
if [[ $cont == y* ]]; then
    prompt "\nTarballing dist files: $builddir to $src"
    tar -C build -hczf $src $prefix
fi

# 1b) Make RNA folding-free version
prompt "\nMake RNA-folding-free version of source: $builddir"
nofoldingdir=${builddir}-nofolding
cp -Trd $builddir ${nofoldingdir}
pushd $nofoldingdir
replace "#EXCLUDE_RNA_FOLDING" "EXCLUDE_RNA_FOLDING" \
    -- setup.sh
replace \
    "#export EXCLUDE_RNA_FOLDING" "export EXCLUDE_RNA_FOLDING" \
    "#controldir" "controldir" \
    -- init.sh
popd

# 2b) tarball relevant source files
src=${distbase}-nofolding.tar.gz
cont=y
if [[ -e $src ]]; then
    read -p "Overwrite $src (y/n)? " cont
fi
if [[ $cont == y* ]]; then
    prompt "\nTarballing dist files: $nofoldingdir to $src"
    tar -C build -hczf $src ${prefix}-nofolding
fi


# Tarball data!
data=${distbase}_data.tar.gz
cont=y
if [[ -e $data ]]; then
    read -p "Overwrite $data (y/n)? " cont
    if [[ $cont == y* ]]; then
	read -p "Are you absolutely positive (y/n)? " cont
    fi
fi
if [[ $cont == y* ]]; then
    prompt "\nTarballing data files to: $data"
    if [[ ! -d data ]]; then
	echo "Error: expected data/ directory" >&2
	exit 1
    fi
    tar --exclude="*.pkl" --exclude="*.2bit" -hczf $data data
fi

# Tarball manuscript results
data=${distbase}_manuscript.tar.gz
cont=y
if [[ -e $data ]]; then
    read -p "Overwrite $data (y/n)? " cont
    if [[ $cont == y* ]]; then
	read -p "Are you absolutely positive (y/n)? " cont
    fi
fi
if [[ $cont == y* ]]; then
    prompt "\nTarballing manuscript files to: $data"
    if [[ ! -d manuscript ]]; then
	echo "Error: expected manuscript/ directory" >&2
	exit 1
    fi
    tar -hczvf $data manuscript/{benchmark,train,*.*} manuscript/validation/{*.*,*/RESULTS} manuscript/validation/*/*.{pcoord,scored}
fi
