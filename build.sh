#!/bin/sh
# This script will try to compile the source code in the subdirectory
# build-$architecture, using one processor core
# you need to have the dependencies fullfilled, see README

arch="$(uname -m)"
bdir="Build-$arch"

test -d "$bdir" || mkdir "$bdir"

echo 'Running cmake'
cd "$bdir"
cmake ../

echo "Compiling source code in $bdir"
make && echo "Programs successfully built in $bdir, check it out!"
