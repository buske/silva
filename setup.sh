#!/usr/bin/env bash

set -eu
set -o pipefail

# Download unafold
wget http://mfold.rna.albany.edu/cgi-bin/UNAFold-download.cgi?unafold-3.8.tar.gz
tar -xzf unafold-3.8.tar.gz
cd unafold-3.8
./configure && make

# Download data


# Install modified version of milk
prefix="$(pwd)"
cd tools/milk
export PYTHONPATH="$prefix:$PYTHONPATH"
python setup.py install --prefix="$prefix"

