#   common.sh - common POSIX shell code used in pipeline scripts
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

pipeline_version() {
	echo '1.5.0' # TODO: take from Git or file
}

# check if binary in path
checkexecutables() {
	for cmd in $@; do
		if [ -z "$(which "$cmd")" ]; then
			echo "'$cmd' not found in PATH." 1>&2
			return 1
		fi
	done
	return 0
}

# portable tmpdir creation
maketempdir() { mktemp -d 2>/dev/null || mktemp -d -t tmp; };  # works on Linux and Mac

# fallback routines for detecting the number of usable CPU cores
detectcores() {
	if checkexecutables nproc 2>/dev/null; then
		nproc
	else
		[ -f /proc/self/status ] || return 1
		checkexecutables grep cut tr bc wc || return 2
		grep -m 1 '^Cpus_allowed:' /proc/self/status | cut -f 2 | tr -d ',' | tr '[:lower:]' '[:upper:]' | xargs echo "obase=2; ibase=16;" | LANG='C' bc | tr -d -c '1' | wc -c
	fi
}

# available RAM in MiB (Linux only)
detectmemory() {
	if checkexecutables getconf 2>/dev/null; then
		echo "$(($(getconf PAGESIZE)*$(getconf _PHYS_PAGES)/1048576))"
		return 0
	fi
	checkexecutables awk || return 2
	awk '{if($1=="MemTotal:"){print int($2/1024 + 0.5); exit}}' /proc/meminfo
}

# turn path into absolute path
absolutepath() {
	if checkexecutables readlink 2>/dev/null; then
		readlink -f "$1"
		return 0
	fi
	(cd "$1" && echo "$PWD")
}

# add folder to PATH variable (export)
addtopath() {
	[ -d "$1" ] && export PATH="$(absolutepath "$1"):$PATH"
}

# create a new output folder with numeric suffix
makeprojectdir() {
	checkexecutables expr || return 2
	local count=0
	local input_stem="$1"
	while [ -e "$input_stem"."$count" ]; do
		count="$(expr $count + 1)"
	done
	mkdir "$input_stem"."$count"
	echo "$input_stem"."$count"
}

# check if python is 2.7+ but not 3.0+
checkpython2() {
	checkexecutables cut tr "$1" || return 2
	python_version="$($1 --version 2>&1 | cut -d ' ' -f 2)"
	IFS='.' read -r python_version_major python_version_minor python_version_patch <<-_EOF_
$python_version
_EOF_

	if [ "$python_version_major" -ne 2 -o "$python_version_minor" -lt 7 ]; then
		echo "Your Python '$python_version_major.$python_version_minor.$python_version_patch' must be at least version 2.7 but not Python 3"
		return 1
	fi
}

# set all variable for a refpack folder
initrefpack() {
	checkexecutables readlink || return 2
	refpack="$1"

	if [ -z "$refpack" ]; then
		echo 'Refpack path must be defined' 1>&2
		return 1;
	fi

	ref_root="$(readlink -f "$refpack")"
	aligner_index="$ref_root/$index_subdir/nuc"
	refdata="$ref_root/refdata.fna"
	refdata_index="$ref_root/refdata.fna.fai"
	mapping="$ref_root/mapping.tax"
	input="$(readlink -f "$input")"
	TAXATORTK_TAXONOMY_NCBI="$ref_root"/ncbi-taxonomy
}

# convert bioboxes.org binning format to simple TAB-separated
binning2taxpath() {
	checkexecutables grep taxknife || return 2
	cat $@ |
	grep -v -e '^@' -e '^#' |
	taxknife -f 2 --mode annotate -s path
}


# create a summary from TAB-separated
taxpath2taxsummary() {
	checkexecutables awk cut sort || return 2
	cat $@ |
	cut -f 2-4 |
	LC_COLLATE='C' sort -k1,1 |
	awk -F '\t' 'BEGIN{ id="%SV" }{  if( $1 != id ) { if( id != "%SV" ) { printf "%s\t%d\t%d\n", id, sum1, sum2; } id=$1; sum1=$2; sum2=$3 } else { sum1+=$2; sum2+=$3 } } END{ if(id!="%SV") printf "%s\t%d\t%d\n", id, sum1, sum2; }'
}

# create a vertical profile from TAB-separated
binning2vprofile() {
	checkexecutables grep cut taxknife sort uniq || return 2
	cat $@ |
	grep -v -e '^@' -e '^#' -e '^$' |
	cut -f 2 |
	taxknife -f 1 --mode annotate -s rank |
	LC_COLLATE='C' sort |
	uniq -c  # depth of assignment profile
}

# create all summary files from binning file
binning2summary() {
	checkexecutables tee || return 2
	local sample_name="$1"
	if [ ! -r "$sample_name".binning ]; then
		echo 'Error: invalid binning file.' 1>&2
		return 1;
	fi

	binning2vprofile "$sample_name".binning > "$sample_name".vprofile

	binning2taxpath "$sample_name".binning | tee "$sample_name".taxpath.tsv |  # human-readable binning output
	  taxpath2taxsummary | tee "$sample_name".taxsummary.tsv |  # cumulative support and amount of sequence data assigned to each taxon
	  taxsummary2krona > "$sample_name".taxsummary.html  # corresponding Krona plot
}

# calculate a hash which is unique for the topology and taxon IDs (not names)
taxonomy_version() {
	checkexecutables awk sort md5sum cut || return 2
	cat "$@" |
	awk -F '\t' '{if($1 != $3) print $1 "\t" $3}' |
	LC_COLLATE=C sort -u |
	md5sum |
	cut -d ' ' -f 1;
}
