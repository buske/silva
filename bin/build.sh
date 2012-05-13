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


prompt "\nReplacing <VERSION> tag with: $new_version"
pushd $builddir
if [[ $new_version != $(cat VERSION) ]]; then
    echo "Error: version mismatch" >&2
    exit 1
fi
sed -e "s/<VERSION>/$new_version/g" -i"" README tools/milk/milk/milk_version.py
popd


src=${distbase}.tar.gz
if [[ ! -e $src ]]; then
    prompt "\nTarballing dist files: $builddir to $src"
    tar -C build -hczf $src $prefix
fi


data=${distbase}_data.tar.gz
if [[ ! -e $data ]]; then
    prompt "\nTarballing data files to: $data"
    if [[ ! -d data ]]; then
	echo "Error: expected data/ directory" >&2
	exit 1
    fi
    tar -hczf $data data
fi
