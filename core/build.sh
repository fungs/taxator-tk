#!/bin/sh
# This script will try to compile the source code in the bin subdirectory,
# using one processor core by default

set -o errexit
set -o nounset

compile_threads=${1:-1}

bdir='bin'

test -d "$bdir" || mkdir "$bdir"

echo 'Running cmake'
cd "$bdir"
cmake ../

echo "Compiling source code in $bdir"

make -j "$compile_threads"
echo "Programs successfully built in $bdir, check it out!"
