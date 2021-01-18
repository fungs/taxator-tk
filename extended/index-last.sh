#!/bin/sh
#   index-last.sh - uses LAST blastdb to build an index (aka database)
#                   for the reference data
#
#   Written in 2014 by Johannes Dröge johannes.droege@uni-duesseldorf.de
#
#   To the extent possible under law, the author(s) have dedicated all copyright
#   and related and neighboring rights to this software to the public domain
#   worldwide. This software is distributed without any warranty.
#
#   You should have received a copy of the CC0 Public Domain Dedication along
#   with this software. If not, see
#   http://creativecommons.org/publicdomain/zero/1.0/

set -o errexit

# Constants
required_programs="readlink dirname basename lastdb"
subdir=aligner-index/last
dbname=nuc

# Use binary folder
base_root="$(readlink -f "$(dirname "$0")")"
export PATH="$base_root/bin:$PATH"

# Check for required programs
for cmd in $required_programs; do
	if test -z "$(which "$cmd")"; then
		echo "'$cmd' not found in PATH."
		exit 1
	fi
done

# Parse command line
refpack="$1"
if test -z "$refpack"; then
	echo "Specify the refpack folder as first argument." 1>&2
	exit 1
fi

# Reference FASTA given?
ref_root="$(readlink -f "$refpack")"
ref_fasta="$ref_root"/refdata.fna
if test ! -e "$ref_fasta" ; then
	echo "No reference FASTA '$ref_fasta' found." 1>&2
	exit 1
fi

# Index already exists?
cpath="$ref_root"/"$subdir"
if test -d "$cpath"; then
	echo "Folder '$cpath' exists, abort." 2>&1
	exit 1
fi

# Build index
echo -n "Creating LAST index for '$ref_root'. "
mkdir -p "$cpath"
cd "$cpath"
lastdb -Q 0 -i 10 "$dbname" "$ref_fasta" 2> "$dbname".makedb.log && echo 'Success.' || echo 'Failed.' 1>&2

