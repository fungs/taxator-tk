# taxator-tk extended version

## Hardware configuration
The software package taxator-tk performs best if the reference data files are
stored on a local disk, not on a network drive, and these can occupy up to
several hundreds of GiB. The aligners usually need some additional space to
store their index. The index of LAST is considerably large (~6x the FASTA size).
The result files created by the binning workflow can also grow up to several
hundreds of MB as they by default save the alignments for later use. In total,
we recommend an absolute minimum of 30 GB of available local storage space.
Additionally, for some of the workflows which generate temporary data, there
should be enough available space in $TMPDIR (or /tmp, if this variable is not
set).


## Software configuration
We require a Linux system with a 2.6 kernel. When using the LAST aligner, you
need Python (>= 2.7, not 3) and a recent version of Perl. The workflow scripts
use some common UNIX tools which should be installed on every system.
The plotting routines using Krona require that Perl is installed on your system.


## Using refpacks
The provided workflows require a refpack for alignment and sequence assignment.
A refpack contains a FASTA file, NCBI taxonomy files and a mapping file which
maps sequence IDs to taxon IDs. We provide example refpacks which we construct
from different sources, please consult the official download site for download
links. To use a refpack, please unpack the refpack in a location with fast disk
access, for instance using

    tar -C "${datapath}" -xJvf "/path/to/${refpack}"

To use your own refpack, you must provide a mapping file linking sequence
IDs to valid taxon IDs. This can be difficult as the NCBI taxonomy is always
changing (see taxator-tk source README for details). Please send us an email if
you require a different or more recent refpack and have problems building it.

Create an index (use BLAST for general purpose and LAST for high-performance,
read the following paragraph)

    "${installpath}/${ttk}/index-blast.bash" "${datapath}/${refpack}"

Note: taxator-tk supports BLAST due to its ability to run on hardware with
limited main memory, like laptop computers. It is quite fast with the
(discontiguous) megablast algorithm and a locall stored index. Consider BLAST
with "megablast" algorithm for much faster but less sensitive alignment if you
want to quickly explore sample taxa and you are not interested so much into the
individual sequences. Use "blastn" algorithm instead of "discontiguous megablast"
for most extensive alignment but this will take a while for large data. See the
tips on how to change the algorithm used. Prefer LAST on large-memory systems
(20 GiB and up) as it provides both good speed and a sensitivity comparable to
the "blastn" algorithm. Better sensitivity in the local alignment will result in
more potential homologs as input to the taxator program and therefore to a larger
fraction of the sample to be assigned to taxa by the workflow. However, building
a LAST index is currently unfeasible for large refpacks, because it can take
weeks or months.

## Usage
We provide shell scripts that implement simple workflows. You can use these
scripts as templates for your own. Simply provide the binning workflow script
with the refpack location as the first argument and the FASTA file that you want
to assign as the second argument. Result files will be written in the current
working directory.

    "${installpath}/${ttk}/binning-blast.bash "${datapath}/${refpack}" contigs.fna output/

Tip: Make sure that the FASTA sequence IDs don't contain whitespace characters
like spaces. These cannot be handled correctly by many alignment programs. To
strip the identifier and to check for unique identifiers, you can use the
supplied fasta-strip-identifier command.

Tip: To change default parameters such as the blast algorithm, prepend the
variable like `blast_algorithm=megablast binning-blast.bash ...`.

Tip: The scripts detect the number of available CPU cores and uses all cores
by default. Use the Linux command "taskset" to control the number of threads
like `taskset -c 1-10 binning-blast.bash ...`. This has the advantage that
the workflow will be restricted to the first 10 cores. Alternatively, you can
prepend the variable 'cores' like `cores=4 binning-blast.bash ...`.
