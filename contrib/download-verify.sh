#/usr/bin/env sh
#
# automatically download external software by last known URLs
# needs wget to work

set -e

checksum() {
  sha256sum | cut -f 1 -d ' '
}

tail -n+2 ./download.tsv |
 while read csum_ok filename url; do
   # download and calculate checksum
   if test -f "$filename"; then
     csum_test="$(checksum < "$filename")"
     method='existing'
   else
     csum_test="$(wget -q -O - "$url" | tee "$filename" | checksum)"
     method='download'
   fi
   
   # compare
   if test "$csum_test" != "$csum_ok"; then
     echo "'$filename': integrity check FAILED!" 1>&2
     exit 1
   else
     echo "'$filename': OK ($method)"
   fi
 done
