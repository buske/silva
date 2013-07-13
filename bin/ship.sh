#!/usr/bin/env bash

set -eu
set -o pipefail

function prompt {
    echo -e "$@" >&2
    read -p "(ENTER to continue)... "
}

name="silva"
WWW=porta:/var/www/html/$name
releasedir=/filer/buske/$name/public_html/release

version=$(cat VERSION)
read -p "Releasing $name $version. Press ENTER to continue..."

prefix=${name}-${version}
builddir=build/$prefix
distbase=dist/$prefix

if [[ ! -d $builddir || 
      ! -e $distbase.tar.gz || 
      ! -e ${distbase}_data.tar.gz ]]; then
    echo "Error: some input files are missing" >&2
    exit 1
fi

rsync -e ssh \
    -vu \
    --links \
    $builddir/{README,HISTORY} \
    $WWW/

rsync -e ssh \
    -vu \
    --links \
    $distbase* \
    $releasedir
