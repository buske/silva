#!/usr/bin/env bash

set -eu
set -o pipefail

name="sympri"
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


echo "Replacing version number with: $new_version" >&2
read -p "(ENTER to continue)... "
if [[ ! -e VERSION ]]; then
    echo "Error: missing VERSION file" >&2
    exit 1
fi
echo $new_version > VERSION


echo "Copying files to build directory: $builddir" >&2
read -p "(ENTER to continue)... "
if [[ ! -e .rsync-filter ]]; then
    echo "Error: expected .rsync-filter file in cwd" >&2
    exit 1
fi
mkdir -pv $builddir
rsync -r \
    -v \
    --links \
    --filter='. .rsync-filter' \
    . \
    $builddir/


echo "Replacing <VERSION> tag with: $new_version" >&2
read -p "(ENTER to continue)... "
pushd $builddir
if [[ $new_version != $(cat VERSION) ]]; then
    echo "Error: version mismatch" >&2
    exit 1
fi
sed -e "s/<VERSION>/$new_version/g" -i"" README
popd


echo "Tarballing dist files: $builddir" >&2
read -p "(ENTER to continue)... "
if [[ ! -d data ]]; then
    echo "Error: expected data/ directory" >&2
    exit 1
fi
src=${distbase}.tar.gz
data=${distbase}_data.tar.gz
if [[ ! -e $src ]]; then
    tar -hczf $src $stagedir
fi
if [[ ! -e $data ]]; then
    tar -hczf $data data
fi
