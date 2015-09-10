# File formats

The taxator-tk suite uses simple, text-based and, if possible, standard intermediate formats. Thereby, data can be easily transformed, edited, compressed and saved.

## Alignments

### Formatting

The input alignments to the program 'taxator' are given as TAB-separated columns with UNIX line breaks ('\n'). Comment lines start with a '#'. Empty lines are also ignored.

### Columns

The fields are defined as

1. query identifier
2. query start position (1-based, inclusive)
3. query stop position (1-based, inclusive)
4. query length
5. reference identifier
6. reference start position (1-based, inclusive)
7. reference stop position (1-based, inclusive)
8. alignment score (positive float, max. 6 digits after the decimal point, larger means better alignment)
9. e-value (float, standard or scientific notation)
10. identities (number of exactly matching positions, positive integer)
11. alignment length (positive integer)
12. alignment CIGAR code (optional, see official CIGAR definition by samtools)

### Notes

* If reference start and stop positions for nucleotide data are swapped, this denotes the reverse complement. Query position swapping is currently unsupported.
* Sequence identifiers must not contain TAB characters. Generally, space characters are allowed but discouraged as they produce problems with many aligners or alignment formats (see MAF format).
* It's ok to fill in evalues of zero if the aligner does not report any such value.
