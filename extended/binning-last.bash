#!/usr/bin/env bash
#   binning-last.bash - sample binning workflow using LAST
#
#   Written in 2014-2021 by Johannes Dröge code@fungs.de
#
#   To the extent possible under law, the author(s) have dedicated all copyright
#   and related and neighboring rights to this software to the public domain
#   worldwide. This software is distributed without any warranty.
#
#   You should have received a copy of the CC0 Public Domain Dedication along
#   with this software. If not, see
#   http://creativecommons.org/publicdomain/zero/1.0/

set -o errexit
set -o nounset

# Parameters
last_options_default=""             # simply use LAST defaults
taxator_speedup_default='0.5'       # 1 is fastest, 0 means using all alignments
taxator_logfile_default='/dev/null' # logging disabled
binner_logfile_default='/dev/null'  # logging disabled

# Command line arguments
refpack="${1-}"
input="${2-}"
working_project="${3-}"

# Constants
sample_name=sample
index_subdir=aligner-index/last
memory_min=5120  # minimum main memory in MB

# Load library and set path
source "${0%/*}/lib/common.sh"
progpath="$(absolutepath "$0")"
addtopath "${progpath%/*}/bin"

# System check
checkexecutables sort time gzip gunzip grep tr cut uniq tee python lastal-parallel taxator binner alignments-filter taxknife lz4
time_cmd="$(which time)"
cores_max="$(detectcores)"
memory_max="$(detectmemory)"
checkpython2 python

# Parse command line
if [ -z "$refpack" -o ! -d "$refpack" ]; then
	echo "Specify the refpack folder as first argument." 1>&2
	exit 1
fi

if [ -z "$input" -o ! -r "$input" ]; then
	echo "Specify your input FASTA file as second argument and make sure that '$input' is readable." 1>&2
	exit 1
fi

# Index built?
if [ ! -d "$refpack"/"$index_subdir" ]; then
	echo "No aligner index found in '$refpack/$index_subdir', create it first." 1>&2
	exit 1
fi

# Hardware suitable?
if [ "$memory_max" -lt "$memory_min" ]; then
	echo "Your system has $memory_max MB of memory available,"
	echo "the minimum is set to $memory_min MB."
	exit 1
fi

# Set refpack-related locations
initrefpack "$refpack"  # set variables: aligner_index, refdata, refdata_index, mapping, input, TAXATORTK_TAXONOMY_NCBI
input_filename="${input##*/}"
export TAXATORTK_TAXONOMY_NCBI  # make taxator-tk programs work with taxonomy

# Set project dir
if [ -n "$working_project" ]; then
	if ! mkdir "$working_project" 2>/dev/null; then
		if ! [ -d "$working_project" ]; then
			echo "Error: '$working_project' is not a valid output folder." 1>&2
			exit 1
		fi
	fi
else
	working_project="$(makeprojectdir "$input_filename")"
fi
echo "Project directory is '$working_project/'"
cd "$working_project"


# WORKFLOW
echo "Aligning sample against sequences in '$refdata' and assigning segments to taxa using ${cores:-$cores_max} threads."

# Align query against reference
compression_cmd='lz4' decompression_cmd='lz4 -d' $time_cmd -p -o lastal-parallel.time lastal-parallel -f 1 -X 3 -P "${cores:-$cores_max}" ${last_options:-$last_options_default} "$aligner_index" "$input" |
lastmaf2alignments-parallel -s |  # convert from MAF to tabular

# Save the alignment in a gzipped tabular format
tee >(gzip > "$sample_name".alignments.gz) |

# Assign query segments to taxa
$time_cmd -p -o taxator.time taxator \
  -a rpa \
  -g "$mapping" \
  -q "$input" \
  -f "$refdata" \
  -i "$refdata_index" \
  -p "${cores:-$cores_max}" \
  -l "${taxator_logfile:-$taxator_logfile_default}" \
  -x "${taxator_speedup:-$taxator_speedup_default}" \
  -o 1  |
sort -k1,1 > "$sample_name".gff3

echo "Assigning whole sequences."
$time_cmd -p -o binner.time binner \
  -n "$input_filename" \
  -l "${binner_logfile:-$binner_logfile_default}" \
  < "$sample_name".gff3 \
  > "$sample_name".binning

echo "Generating summary files."
binning2summary "$sample_name"

echo "Results are in '$working_project/'."
