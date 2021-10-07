# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

**Notes**
There are no recorded dates and only record type `Changed` prior to v1.4.0.
Also, no version tags existis prior to version 1.0.0. Pipeline scripts were only
added to this repository starting from version 1.4.0.

## [1.4.1] - 2021-10-03
### Fixed
* core: compile with recent compiler and Boost CMake detection

## [1.4.0] - 2021-01-26
### Added
* all: reorganized code and included extended pipeline scripts and documentation

### Changed
* core: backport of code changes from 1.5 branch
* core: add version flags to binaries
* core: seqan update to v2.4.0
* core: tree.hh update to v3.7
* extended: correctly detect Python2 version with two digits in patch number
* extended: update KronaTools to v2.7.1
* extended: update GNU parallel to v20201222
* extended: align with LZA v1.9.3
* extended: align with NCBI Blast 2.11.0+
* extended: align with last-align v1170
* extended: parallel index building and ambiguous characters with last aligner

## [1.3.3]
### Changed
* core: fix crashes in all LCA algorithms
* core: more consistent behavior of command line parameters in taxknife
* core: better error reporting

## [1.3.1]
### Changed
* extended: set more options via variable with LAST and BLAST
* extended: refactor shell code for cleaner workflows and
            improve POSIX compatibility
* extended: better CPU core detection via numproc command (supports containers)
* extended: specification of output folder (backward-compatible behavior)

## [1.3.0]
### Changed
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

## [1.2.2]
### Changed
* extended: make reference data and aligner index location independent of
*           installation root (easier to run it in multiuser environments)
* extended: remove integrity check for data (should become part of the refpack)
* extended: rename and clean all scripts
* extended: change standard BLAST algorithm to "discontiguous megablast"
            for better sensitivity (increases the runtime by ~20x)

## [1.2.1]
### Changed
* extended: fix bug in binning-workflow-fasta-last.sh, reference FASTA index wasn't used
* extended: replaced quickstart.pdf by a text-based README file

## [1.2.0]
### Changed
* core: program taxknife substitutes name-filter and rank-filter
* core: adjust data types to support longer segments
* core: better approximation of identity scores + support values in binner
* core: optimizations in RPA procedure
* core: heuristic to avoid excessive runtime with bad alignments
* core: workaround for CMake bug with boost threads linking

## [1.1.1]
### Changed
* extended: fix typo bug in create-index-last.sh
* extended: quickstart manual update

## [1.1.0]
### Changed
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

## [1.0.0]
### Changed
* core: publication release, many functionality added;
  evaluations done on this version
* core: optimized for large-memory systems
* extended: sample workflow shell scripts for taxonomic binning
  includes microbial RefSeq with taxonomy
* extended: large prepared LAST index with download script

## [0.3]
### Changed
* core: bug fix for noise filtering in binner
* core: tools reduced to a core set around taxator (previously predictor) and binner
  (prevously predictions-filter)
* core: multi-threading support in taxator
* core: binner supports user defined identity constrains at fixed ranks

## [0.2]
### Changed
* core: predictions-filter for combination of predictions by mutiple methods/experts
  and to clean up the predictions of natural samples where we require a certain
  number of read for each taxonomic group to be considered

## [0.1]
### Changed
* core: functional framework with support for various LCA methods and fast aligners

[Unreleased]: https://github.com/fungs/taxator-tk/compare/v1.3.3...HEAD
[1.4.0]: https://github.com/fungs/taxator-tk/compare/v1.3.3...v1.4.0
[1.3.3]: https://github.com/fungs/taxator-tk/compare/v1.3.2...v1.3.3
[1.3.1]: https://github.com/fungs/taxator-tk/compare/v1.3.1...v1.3.0
[1.3.0]: https://github.com/fungs/taxator-tk/compare/v1.2.2...v1.3.0
[1.2.2]: https://github.com/fungs/taxator-tk/compare/v1.2.1...v1.2.2
[1.2.1]: https://github.com/fungs/taxator-tk/compare/v1.2.0...v1.2.1
[1.2.0]: https://github.com/fungs/taxator-tk/compare/v1.1.1...v1.2.0
[1.1.1]: https://github.com/fungs/taxator-tk/compare/v1.1.0...v1.1.1
[1.1.0]: https://github.com/fungs/taxator-tk/compare/v1.0.0...v1.1.1
[1.0.0]: https://github.com/fungs/taxator-tk/releases/tag/v1.0.0
