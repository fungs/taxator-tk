This file only lists major changes and bug fixes for each version

v1.3.3 taxator-tk
=================
* core: fix crashes in all LCA algorithms
* core: more consistent behavior of command line parameters in taxknife
* core: better error reporting

v1.3.1 taxator-tk
=================
* extended: set more options via variable with LAST and BLAST
* extended: refactor shell code for cleaner workflows and
            improve POSIX compatibility
* extended: better CPU core detection via numproc command (supports containers)
* extended: specification of output folder (backward-compatible behavior)

v1.3.0 taxator-tk
=================
* core: rewrite of RPA with better handling of input alignments and
        speed improvements
* core: C++11 compatibility
* core: exception handling
* core: bioboxes output format
* core: suppress some notifications
* core: much better scaling core algorithm due to homolog selection heuristic,
        there is no need to restrict the number of input alignments any more
* extended: CPU core auto-detection, now working with Linux taskset/cpusets
* extended: support for multi-core alignment with LAST
            (lastal-parallel wrapper script)
* extended: support of parallel search with LAST by integration of GNU parallel
            in script lastal-parallel, uses lz4 compression of temporary data
* extended: faster conversion by integration of GNU parallel in script
            lastmaf2alignments-parallel
* extended: integration of KronaTools for interactive pie chart plotting
* extended: download script for sample refpack removed, the project page will
            point to a location from where to download refpacks

v1.2.2 taxator-tk
=================
* extended: make reference data and aligner index location independent of
*           installation root (easier to run it in multiuser environments)
* extended: remove integrity check for data (should become part of the refpack)
* extended: rename and clean all scripts
* extended: change standard BLAST algorithm to "discontiguous megablast"
            for better sensitivity (increases the runtime by ~20x)

v1.2.1 taxator-tk
=================
* extended: fix bug in binning-workflow-fasta-last.sh, reference FASTA index wasn't used
* extended: replaced quickstart.pdf by a text-based README file

v1.2.0 taxator-tk
=================
* core: program taxknife substitutes name-filter and rank-filter
* core: adjust data types to support longer segments
* core: better approximation of identity scores + support values in binner
* core: optimizations in RPA procedure
* core: heuristic to avoid excessive runtime with bad alignments
* core: workaround for CMake bug with boost threads linking

v1.1.1 taxator-tk
=================
* extended: fix typo bug in create-index-last.sh
* extended: quickstart manual update

v1.1.0 taxator-tk
=================
* core: code cleanup
* core: SEQAN update with faster alignments (>2x speed improvement)
* core: indexed FASTA access (now runs on computers with much less RAM)
* core: code speedup by factor of two
* core: clean up command line options
* core: improve packaging scripts
* extended: workflow on workstations with limited memory
            (new indexed FASTA feature in taxator; with BLAST)
* extended: switch to NCBI BLAST+ (nucleotide megablast) as default,
            option to use LAST
* extended: general cleanup/improvement of included workflow scripts;
            no Python required any more in BLAST workflow

v1.0.0 taxator-tk
=================
* core: publication release, many functionality added;
  evaluations done on this version
* core: optimized for large-memory systems
* extended: sample workflow shell scripts for taxonomic binning
  includes microbial RefSeq with taxonomy
* extended: large prepared LAST index with download script

v0.3 taxator-tk
===============
* core: bug fix for noise filtering in binner
* core: tools reduced to a core set around taxator (previously predictor) and binner
  (prevously predictions-filter)
* core: multi-threading support in taxator
* core: binner supports user defined identity constrains at fixed ranks

v0.2 taxator-tk
===============
* core: predictions-filter for combination of predictions by mutiple methods/experts
  and to clean up the predictions of natural samples where we require a certain
  number of read for each taxonomic group to be considered

v0.1 taxator-tk
===============
* core: functional framework with support for various LCA methods and fast aligners
