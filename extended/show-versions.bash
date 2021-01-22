#!/usr/bin/env bash
#   show-versions.bash - report all program versions
#
#   Written in 2021 by Johannes Dröge code@fungs.de
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

# Load library and set path
source "${0%/*}/lib/common.sh"
progpath="$(absolutepath "$0")"
addtopath "${progpath%/*}/bin"

# helper for extracting version strings in Bash
extract_version() {
  ([[ "$1" =~ $2 ]] && echo "${BASH_REMATCH[1]}")
}

echo -e "This is the taxator-tk extended pipeline with the following program versions\n"

# taxator-tk extended pipeline
echo "taxator-tk pipeline scripts: v$(pipeline_version)"

# makeblastdb
regex='^makeblastdb:\ (.+)$'
checkexecutables makeblastdb && echo -n 'makeblastdb: ' && extract_version "$(makeblastdb -version | head -n 1)" "$regex"

# blastn
regex='^blastn:\ (.+)$'
checkexecutables blastn && echo -n 'blastn: ' && extract_version "$(blastn -version | head -n 1)" "$regex"

# blastp
regex='^blastp:\ (.+)$'
checkexecutables blastp && echo -n 'blastp: ' && extract_version "$(blastp -version | head -n 1)" "$regex"

# lastdb
regex='^lastdb (.+)$'
checkexecutables lastdb && echo -n 'lastdb: ' && extract_version "$(lastdb --version 2>&1)" "$regex"

# lastal
regex='^lastal (.+)$'
checkexecutables lastal && echo -n 'lastal: ' && extract_version "$(lastal --version 2>&1)" "$regex"

# LZ4
regex='^\*\*\*\ LZ4\ command\ line\ interface\ 64-bits\ (.+)\,\ by\ Yann\ Collet\ .*\*\*\*$'
checkexecutables lz4 && echo -n 'lz4: ' && extract_version "$(lz4 --version 2>&1)" "$regex"

# parallel
regex='^GNU parallel (.+)$'
checkexecutables parallel && echo -n 'parallel: ' && extract_version "$(parallel --version | head -n 1)" "$regex"

echo -e "\nPlease record and publish all program versions with your results."
