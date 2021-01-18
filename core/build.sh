#!/bin/sh
# This script will try to compile the source code in the subdirectory
# build-$architecture, using one processor core
# you need to have the dependencies fullfilled, see README

set -o errexit
set -o nounset

compile_threads=${1:-1}

arch="$(uname -m)"
bdir="Build-$arch"

test -d "$bdir" || mkdir "$bdir"

echo 'Running cmake'
cd "$bdir"
cmake ../

echo "Compiling source code in $bdir"

make -j "$compile_threads"
echo "Programs successfully built in $bdir, check it out!"
