# Usage

Basic concept:
The basic pipeline concept allows splitting the overall prediction into
individual sub-tasks that are implemented in separate programs. Intermediate
data is text-based and consists of lines with TAB-separated fields. At each step
of the process, the user can decide whether to save an intermediate file, e. g.
compressing it using a standard compressor like gzip, or to feed it directly
into the next program without the necessity to save the data on a disk.
A binning process can be implemented as follows:


**(1) ALIGN [=> (2) CONVERT] [=> (3) FILTER] => (4) PREDICT SEGMENTS => (5) BINNING**

## 1. **ALIGN**
Align sample sequences against a database. The aligner can be replaced
as long a conversion tool is provided that will convert the aligners
output into the intermediate alignment format used in the pipeline.
Aligners supported out-of-the-box are
- NCBI blastn (ftp://ftp.ncbi.nih.gov/blast/executables/blast+/)
- NCBI blastx (ftp://ftp.ncbi.nih.gov/blast/executables/blast+/)
- LAST lastal (http://last.cbrc.jp)

**Note:**
We recommend LAST for aligning nucleotide reads and contigs. The aligner is
extremely fast but requires a large amount of main memory. You can align a
real-sized sample within hours against all of microbial RefSeq. The current NCBI
BLAST+ suite also works but isn't as fast. However, it also runs on low-memory
systems. The BLAST index files should not be stored on slow disk storage such
as an NFS drive with high latency. The new BLAST+ has the advantage that it can
directly generate the alignments in the required format and order, making any
conversion or sorting steps superfluous (see next section).

For database construction with LAST, use only the accession to index the
sequences, which means to reduce the header in the corresponding FASTA file.
The LAST README shows how to transform this sequence file.

## 2. **CONVERT**
Transforms alignments into the universal TAB-separated format. The intermediate
format is quite simple: a tab-separated text file. We provide scripts for use
with "NCBI BLAST" (using XML output) and "LAST" (using MAF output) in the
alignment stage.

- **Native BLAST**:
With the new NCBI BLAST+ packages you can produce the taxator-tk intermediate
alignment format without conversion using the following parameter:
`-outfmt '6 qseqid qstart qend qlen sseqid sstart send bitscore evalue nident length'`

- **blastxml2alignments**:
Blast XML output converter (use "-outfmt 5" with BLAST+). Check out the options
using the "-h" option of blastn, blastx or blastall. Using this output is 
slower than directly generating the TAB-separated intermediate short format and
we highly recommend you to compress it when you want to store it.

- **lastmaf2alignments**:
LAST MAF format converter (use "-f 1"). As a limitation of the MAF format used
by LAST (http://genome.ucsc.edu/FAQ/FAQformat.html), sequence identifiers of
database and query sequences must not contain white spaces or the first white
space separated part must be uniquely trackable for you because the rest will be
discarded when reporting the alignments. LAST output is not sorted by query
identifier, so you must sort the output of the conversion script. The UNIX sort
command like sort from the GNU coreutils is convenient for doing this (set
LC_COLLATE="C" if you are using a non-english locale und Linux) although a
better implementation would only reorder the individual search batches in the
output of LAST.

**Note:**
The conversion scripts are written in python (>=2.7, not python 3). For parsing
BlastXML output, you also need a recent version of BioPython installed. BLAST
can output the alignments consecutively so they can be piped directly
into the program taxator without conversion. Note also that, depending on how
the BLAST DB was constructed, sometimes it will give the Genbank identifier (GI)
as sequence name and sometimes the full identifier. To get the correct sequence
names for your alignments (the same as in your mappings file), it is best to
reformat your DB with the correct options, or in case of XML output, try the
command line options of the conversion script blastxml2alignments. In any case,
the reported sequence names must match the FAST sequence identifiers given to
taxator and the one in the mapping file.

## 3. **FILTER**
Apply restriction filters, e.g. to exclude bad alignments. This step is optional
and ususally not needed but can help you to implement your quality policies or
to reduce the run-time.

- **alignments-filter**:
Support some filtering schemes based on absolute and relative score values or
e-values. Check out the options.

**Note:**
alignments-filter is a simple way to limit the run time of segment prediction by
fixing the number of alignments given to it because the run-time of taxator with
the RPA algorithm is in O(n) where n is the number of alignments for each
reference segment. Starting with version 1.3 of taxator-tk, there is a heuristic
filter integrated into taxator to exclude alignments not relevant for taxonomic
assignment which can be controlled by a parameter (-x) between 0 and 1 where
zero means to use all alignments and 1 to do a maxiumum filtering.


## 4. **PREDICT SEGMENTS**
Make a taxonomic classification for regions of homology. There can be multiple
such regions on a contiguous stretch of a query sequence.

- **taxator**:
Does the actual taxonomic prediction by mapping the alignment reference
sequences on the NCBI taxonomy tree. taxator-tk implements its own realignment
placement algorithm (RPA) of taxonomic classification of sequence segments that
is based on a number of pairwise alignments and strikes to approximate placement
via phylogenetic methods. This method requires some calculations but will give
conservative/reliable results with a low fraction of false predictions. As
a general and fast way, various LCA methods (similar to those in MEGAN) are also
implemented in taxator.

Note that, if you do not specify options for sequence index files for any of
query or reference, the sequences get loaded into memory. Starting from
version 1.2, sequences can be accessed directly on the disk storage by
specifying a FASTA index file (usually named like dna.fna.fai). If this file
doesn't exist, it will be generated automatically and saved for later use. For
this feature to work nicely, make sure that your disk storage is fast and has
low enough latency (e.g. generally avoid network storage). If you want to write
taxator's log file, be warned that these are quite large and can grow to several
GiBs.

## 5. **BINNING**
Use a strategy to combine one or more segment predictions and to assign entire
query sequences. Binning represents a common means to use the segment
assignments. Others could be the detection of assembly errors or to look for
unknown horizontal gene transfers.

- **binner**:
This program takes all the segment predictions of a sample and assigns whole
sequences. It provides a whole sample noise filtering step where very rare
taxa are pruned from the predictions (similar to MEGAN's min-support filter).
This only makes sense when you are processing a real sample and your sequences
are not independent of each other! Additionally, the user can specify a minimum
percentage identity for any named rank. By default, no such value is chosen but
it seems reasonable to require a (longer) sequence to be for example at least
70% identical to some known sequence of the same species. The percentage
identity is an abstract underestimate of the usual percentage identify. It is
abstract because it allows to puzzle together pieces of genomes of the same
taxon and it is an underestimate, because it does not consider gap regions in
between segments which would produce by-chance matches in a global pairwise
alignment of the query sequence. Try to set generous limits to lower the rate of
false predictions at specific taxonomic ranks. Overall, the parameters of the
current binning program seem to have very little effect on the final results,
suggesting that the algorithm is quite robust.


## Extra tools

- **taxknife**:
  Works on a given taxonomy tree. Available modes are:
  - traverse:   Traverses taxonomic ids up to a given rank.
  - annotate:  Prints the name or rank name of a given taxon.
  - tree: Generate a Newick tree.

# NCBI taxonomy

All programs that need a taxonomy to work will read the NCBI taxonomy dump
files from a directory that is defined via the environment variable
TAXATORTK_TAXONOMY_NCBI. You can download the taxonomy package
ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz but keep in mind that
the taxonomy files are modified on a regular basis. Place the files names.dmp
and nodes.dmp in a folder and specify its path in the environment variable.
An environment variable is only persistent for the current session. It must be
reset every time you open a new shell or command line terminal.

**Note:**
Make sure you have a mapping file that maps each sequence in your sequence
collection to a valid NCBI taxonomic identifier. We provide such mappings
for several refpacks. A refpack is a pre-built reference collection with a
fitting taxonomy version and a corresponding mapping from sequence identifier
to taxon ID. They can also be found in the RefSeq release catalogue
but should be checked to match the actual taxonomy. If you build your own
refpack, tt is crucial that the mapping exactly! maps the reference sequences
names, as reported by the aligner, to taxonomic IDs, as specified in the
taxonomy. It is easy to check for unmappable reference taxonomic IDs by running
the program taxknife in "annotation" mode with a list of IDs and looking for the
reported errors:

    taxknife -f 2 --mode annotate < mapping.tax >/dev/null

Invalid IDs can be found in the merged.dmp dump file or online on
http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html

Usual operation in taxator is on a pruned taxonomy with major ranks, so
transform the original mapping by running:

    taxknife -f 2 --mode traverse -r species genus family order class phylum superkingdom < mapping.tax > newmapping.tax

# Recipes

Align sequences to DATABASE using LAST and do a best-hit classification

    lastal -f 1 PATH/TO/DATABASE mysample.fna | lastmaf2alignments | sort -k1,1 | taxator -a n-best-lca -n 1 -g acc_taxid.tax

Do a BLAST-MEGAN-style classification with blastn alignments and with a top-percent
filter value of 30 % and a minimum E-value of 0.01.

    blastn -task blastn -db DATABASE --outfmt '6 qseqid qstart qend qlen sseqid sstart send bitscore evalue nident length' -query mysample.fna | taxator -a megan-lca -t 0.3 -e 0.01 -g acc_taxid.tax

If you want to save the alignments with LAST compressed and in sorted order

    lastal -f 1 DATABASE mysample.fna | lastmaf2alignments | sort -k1,1 | gzip > my.alignments.gz

and continue with a RPA prediction and saving the predictions into a file using 10 parallel threads

    zcat my.alignments.gz | taxator -a rpa -q query.fna -f ref.fna -g acc_taxid.tax -p 10 > my.predictions.unsorted.gff3

Or doing all at once without compression-decompression in BASH

    lastal -f 1 DATABASE mysample.fna | lastmaf2alignments | sort -k1,1 | tee >(gzip > my.alignments.gz) | taxator -a rpa -q query.fna -f ref.fna -g acc_taxid.tax -p 10 > my.predictions.unsorted.gff3

Using e-values with LAST is more tricky because it requires another program.
evalue filtering is discouraged with taxator-tk because RPA is robust to spurious alignments.

    lastal -f 1 DATABASE mysample.fna | lastex -z 1 DATABASE.prj mysample.prj - | lastmaf2alignments | ...

Create a binning from the segment assignments (sorting fixes output disorder in taxator multi-threading output)
require same genus and below to be at least 60 % identical to closest homolog.

    sort -k1,1 my.predictions.unsorted.gff3 | binner -i genus:0.6 > my.tax

Putting together above commands one could write the following simple binning workflow (in BASH)

    lastal -f 1 DATABASE my.fna | lastmaf2alignments | sort -k1,1 | tee <(gzip > my.alignments.gz) | taxator -a rpa -g acc_taxid.tax -q query.fna -v query.fna.fai -f ref.fna -i ref.fna.fai -p 10 | sort -k1,1 > my.predictions.gff3
    binner < my.predictions.gff3 -i genus:0.6 > my.tax

Show the corresponding predictions with taxon names

    taxknife -f 2 --mode annotate -s name < my.tax

Or if you want to show the names on the phylum level

    taxknife -f 2 --mode traverse -r phylum < my.tax | taxknife -f 2 --mode annotate -s name

To get a nice rank profile

    taxknife -f 2 --mode annotate -s rank < my.tax | sort | uniq -c

Clean up your mapping so it only points to major ranks (recommended for taxator)

    taxknife -f 2 --mode traverse -r species genus family order class phylum superkingdom < mapping.tax > mapping_cleaned.tax

Restrict run time by using only top 50 alignments

    zcat my.alignments.gz | alignments-filter -b 50 | taxator ...

A typical efficient pipeline with BLAST+ and taxator on 4 CPU cores

    blastn -task blastn -db DATABASE --outfmt '6 qseqid qstart qend qlen sseqid sstart send bitscore evalue nident length' -query mysample.fna -num_threads 2 | taxator -g acc_taxid.tax -q query.fna -v query.fna.fai -f ref.fna -i ref.fna.fai -p 3 | sort -k1,1 > my.predictions.gff3
    binner < my.predictions.gff3 > my.tax

# Taxonomic placement algorithms

- **rpa**: Use the LCA of taxa of a dynamic set of reference segments as determined by
  the ordering of pairwise scores of alignments.
- **lca**: Use the Lowest Common Ancestor (LCA) of all reference taxa of the alignments.
- **megan-lca**: Use the LCA of taxa of only the x % best matching references.
- **ic-megan-lca**: Like megan-lca but considers "unclassified" nodes in taxonomy.

# Tips

- Use taxator with FASTA indices, e.g. on a solid state drive
- Avoid spaces in the sequence identifiers (compatability problems with many aligners)
- Use short sequence identifiers for smaller data files
- Adjust the number of alignments as input to your sample sizes and make a test
  run with different numbers to ensure not loosing many details.
- Make sure to use the C locale or sorting might not give the correct sorting.
  You can force correct sorting by setting the environment variable like this:
  "export LC_COLLATE=C", ensure there is enough availabe space under /tmp
  when using sort with large data. You can use the -T parameter of sort (see
  sort manpage). Use sort from recent GNU coreutils (speed improvements).
  Sorting can also be parallelized.

# Technical aspects

The software implementation and design try to be as slick, fast and flexible as
possible to provide a real tool kit for building an analysis pipeline or to be
integrated in such. The following considerations were taken:

- All tree operations, especially the lowest common ancestor finding
  are implemented in constant time.
- The RPA algorithm depends on sequence random access. To ensure fast
  subsequences retrieval, data can be held in memory or accessed on-disk via
  indexing.
