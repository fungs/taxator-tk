#!/bin/bash
#   binning-workflow-fasta-blast.sh - sample binning workflow using BLAST
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
set -o nounset

# Parameters
blast_options_default=''               # simply use NCBI BLAST+ defaults
blast_algorithm_default='dc-megablast' # megablast: insensitive, fast (test, profile)
                                       # dc-megablast: more sensitive than megablast
                                       # blastn: sensitive, slow (many cores, low mem)
                                       # blastn-short: for short reads only
                                       # tblastx: protein space, very slow (small sample)
                                       # !note: tblastx is the executable, not parameter
taxator_speedup_default='0.5'          # 1 is fastest, 0 means using all alignments
taxator_logfile_default='/dev/null'    # logging disabled
binner_logfile_default='/dev/null'     # logging disabled

# Command line arguments
refpack="${1-}"
input="${2-}"
working_project="${3-}"

# Constants
sample_name=sample
index_subdir=aligner-index/blast

# Load library and set path
source "${0%/*}/lib/common.sh"
progpath="$(absolutepath "$0")"
addtopath "${progpath%/*}/bin"

# System check
checkexecutables ls sort gzip gunzip grep tr cut uniq tee blastn taxator binner taxknife
cores_max="$(detectcores)"

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
	echo "No aligner index found  in '$refpack/$index_subdir', create it first." 1>&2
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
blastn \
  -task "${blast_algorithm:-$blast_algorithm_default}" \
  -db "$aligner_index" \
  -outfmt '6 qseqid qstart qend qlen sseqid sstart send bitscore evalue nident length' \
  -query "$input" \
  -num_threads "${cores:-$cores_max}" \
  ${blast_options:-$blast_options_default} |
tr -d ' ' |  # fields in blast output contain spaces which must be removed

# Save the alignment in a gzipped tabular format
tee >(gzip > "$sample_name".alignments.gz) |

# Assign query segments to taxa
taxator \
  -a rpa \
  -g "$mapping" \
  -q "$input" \
  -f "$refdata" \
  -i "$refdata_index" \
  -p "${cores:-$cores_max}" \
  -l "${taxator_logfile:-$taxator_logfile_default}" \
  -x "${taxator_speedup:-$taxator_speedup_default}" \
  -o 0  |
sort -k1,1 > "$sample_name".gff3

echo "Assigning whole sequences."
binner \
  -n "$input_filename" \
  -l "${binner_logfile:-$binner_logfile_default}" \
  < "$sample_name".gff3 \
  > "$sample_name".binning

echo "Generating summary files."
binning2summary "$sample_name"

echo "Results are in '$working_project/'."

