# Folder structure
In the `contrib` folder, you can unpack all external dependencies. If all unmodified and exact versions are used, the extended pipeline will work by pointing symbolic links into the corresponding folders here.

## Files
* `download.tsv`: file names, download URL and corresponding sha256 hashes to check file integrity for all versions used in the official release.
* `download-verify.sh`: a script to download and verify the external software packages

## Prepare
You need to unpack the downloaded files. Software which is distributed in source code must be compiled with all necessary dependencies resolved. This is a manual process.
