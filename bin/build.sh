#!/usr/bin/env bash

set -eu
set -o pipefail

function prompt {
    echo -e "$@" >&2
    read -p "(ENTER to continue)... "
}

name="silva"
old_version=$(cat VERSION)
echo "Found version: $old_version" >&2
read -p "Enter the release version number (X.Y.Z)... " new_version
if [[ $new_version != *.*.* ]]; then
    echo "Error: invalid version number" >&2
    exit 1
fi


prefix=${name}-${new_version}
builddir=build/$prefix
distbase=dist/$prefix

mkdir -pv build dist

prompt "\nReplacing version number with: $new_version"
if [[ ! -e VERSION ]]; then
    echo "Error: missing VERSION file" >&2
    exit 1
fi
echo $new_version > VERSION


prompt "\nCopying files to build directory: $builddir"
if [[ ! -e .rsync-filter ]]; then
    echo "Error: expected .rsync-filter file in cwd" >&2
    exit 1
fi
mkdir -pv $builddir
rsync -r \
    -v \
    -u \
    --links \
    --filter='. .rsync-filter' \
    . \
    $builddir/


prompt "\nReplacing <VERSION> tag and coded versions with: $new_version"
pushd $builddir
if [[ $new_version != $(cat VERSION) ]]; then
    echo "Error: version mismatch" >&2
    exit 1
fi
replace "<VERSION>" $new_version -- README
popd


src=${distbase}.tar.gz
cont=y
if [[ -e $src ]]; then
    read -p "Overwrite $src (y/n)? " cont
fi
if [[ $cont == y* ]]; then
    prompt "\nTarballing dist files: $builddir to $src"
    tar -C build -hczf $src $prefix
fi

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
    tar -hczf $data data
fi
