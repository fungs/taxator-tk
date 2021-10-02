# taxator-tk extended version

## Hardware configuration
The software package taxator-tk performs best if the reference data files are stored on a local disk, not on a network drive, and these can occupy up to several hundreds of GiB. The aligners usually need some additional space to store their index. The index of LAST is considerably large (~6x the uncompressed FASTA size). The result files created by the binning workflow can also grow up to several hundreds of MB as they by default also store the alignments for subsequent use. In total, we recommend an absolute minimum of 30 GB of available local storage space. Additionally, for some of the workflows which generate temporary data, there should be enough available space in $TMPDIR (or /tmp, if this variable is not set).


## Software configuration
We recommend a Linux system with a minimum 2.6 kernel. When using the LAST aligner, you currently also need Python >= 2.7 and a recent version of Perl. The workflow scripts use some common UNIX tools which should be installed on every standard Linux system. The plotting routines using Krona also require Perl on your system.


## Using refpacks
The provided workflows require a refpack for alignment and sequence assignment. A refpack contains a FASTA file, NCBI taxonomy files and a mapping file which maps sequence to taxon identifiers. We provide example refpacks which we construct from different sources. Please consult the official download site for download
links. To use a refpack, please unpack the compressed refpack archive in a location with fast disk access, for instance using

    tar -C "${datapath}" -xJvf "/path/to/${refpack}"

If you want to use your own refpack, you must provide a mapping file linking sequence IDs to valid taxon IDs. This can somtimes be difficult as the NCBI taxonomy is always changing (see taxator-tk source README for details). Please send us an email if you require a different or more recent refpack and have problems building it.

Create an index (use BLAST for general purpose and LAST for high-performance, read the following paragraph)

    "${installpath}/${ttk}/index-blast.bash" "${datapath}/${refpack}"

Note: taxator-tk extended implements standard BLAST pipelines because it's general purpose and it runs on hardware with limited main memory, like laptop computers. It is quite fast with the (discontiguous) megablast algorithm and a locally stored index. Consider BLAST with "megablast" algorithm for much faster but less sensitive alignment if you want to quickly explore sample taxa and you are not interested so much into the individual sequences. Use "blastn" algorithm instead of "discontiguous megablast" for most extensive alignment but this will take a while for large data. Similar parameters apply to the protein search pipeline. Please see the tips on how to change the algorithm used. On large-memory systems (20 GiB and up), if the index building is feasible, you can also use LAST as it provides both good speed and a sensitivity comparable to the "blastn" algorithm. Better sensitivity in the local alignment will result in more potential homologs as input to the taxator program and therefore to a larger fraction of the sample to be assigned to taxa by the workflow. However, building a LAST index is currently unfeasible for large refpacks, because it can take weeks or months, even in multithreaded mode.

## Usage
We provide shell scripts that implement simple workflows. You can use these scripts as templates for your own. Simply provide the binning workflow script with the refpack location as the first argument and the FASTA file that you want to assign as the second argument. Result files will be written in the current working directory.

    "${installpath}/${ttk}/binning-blast.bash "${datapath}/${refpack}" contigs.fna output/

Tip: Make sure that the FASTA sequence IDs are short and don't contain whitespace characters (like spaces or tabs). These are handled differentyly by alignment programs and lead to problems, and they tend to enlarge the memory footprint. To strip the identifier and to check for unique identifiers, you can use the supplied `fasta-strip-identifier` command.

Tip: To change default parameters such as the blast algorithm, prepend the variable like `blast_algorithm=megablast binning-blast.bash ...`.

Tip: The scripts detect the number of available CPU cores and uses all cores by default. Use the Linux command "taskset" to control the number of threads like `taskset -c 1-10 binning-blast.bash ...`. This has the advantage that the workflow will be restricted to the first 10 cores. Alternatively, you can prepend the variable 'cores' like `cores=4 binning-blast.bash ...`.
