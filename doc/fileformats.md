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

## Segment predictions

The taxonomic predictions for sequence segments that make up the output of taxator are given in General Feature Format (GFF3). Please to the [official documentation](http://www.sequenceontology.org/gff3.shtml) of this format. taxator-tk defines the following fields, some might be better understandable by considering [figure 2 in article
doi://10.1093/bioinformatics/btu745](http://bioinformatics.oxfordjournals.org/content/31/6/817/F2.expansion.html).

1. Query sequence identifier
2. Generating program/algorithm, here `taxator-tk`
3. Type of feature, here `sequence_feature`
4. Begin of feature (1-based, inclusive)
5. End of feature (1-based, inclusive)
6. Prediction score, always `0` until taxator-tk 1.4
7. Strand (not applicable), here `.`
8. Phase (not applicable), here `.`
9. Key-value attribute list with the following reserved tags:
   * `seqlen`: length of query sequence
   * `rtax`: NCBI taxon ID of the most similar alighment match (taxon S)
   * `tax`: a range of the form `low:support-high:support` or `low-high:support` where low corresponds to taxon X and high to taxon R in panel (a)
   * `ival`: a floating point interpolated relative score value in the range \[0,1\] that tell you whether the query is closer to the low (0) or high (1) taxon.

### Notes
* The calculation of the prediction score (column 6) was suspended for technical reasons prior to taxator-tk version 1.1 and should be revived until version 1.5, again.
* The repetition of some fields and tags currently wastes some disk space. However, the GFF3 file is quite small compared to alignments and can be compressed using gzip or similar. For better tracking of information in the prediction part, we might introduce feature identifiers in the future.
* rtax was added for version 1.4 to enable a nearest-neighbor classification scheme or a mixture of schemes in the consensus binning algorithm.
