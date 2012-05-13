#!/usr/bin/env bash

set -eu
set -o pipefail

function prompt {
    echo -e "$@" >&2
    read -p "(ENTER to continue)... "
}

name="silva"
WWW=porta:/var/www/html/$name

read -p "Enter the version number to release (X.Y.Z)... " version
if [[ $version != *.*.* ]]; then
    echo "Error: invalid version number" >&2
    exit 1
fi

prefix=${name}-${version}
builddir=build/$prefix
distbase=dist/$prefix

if [[ ! -d $builddir || ! -e $distbase.tar.gz || ! -e ${distbase}_data.tar.gz ]]; then
    echo "Error: some input files are missing" >&2
    exit 1
else
    rsync -e ssh \
	-vu \
	--links \
	$builddir/{README,HISTORY} $distbase.tar.gz ${distbase}_data.tar.gz \
	$WWW/
fi
