# taxator-tk
A set of programs for the taxonomic analysis of nucleotide sequence data

## Introduction
This is the source code for the taxator toolkit. If you are looking for pre-packaged
programs or sample refpacks for taxonomic assignment, please go to the
[official download page](http://research.bifo.helmholtz-hzi.de/webapps/wa-download/).

This repository consists of two parts, the actual taxator-tk C++ code in the folder `core` and the pipeline scripts which are used to build the ready-to-use binary distribution in folder `extended`.

All issues can be filed in the
[official GitHub repository](https://github.com/fungs/taxator-tk/issues/).

## Documentation
Please see separate documentation:
* [global build and install guide](INSTALL.md)
* [changelog](CHANGELOG.md)
* [core tools compilation guide](core/BUILD.md)
* [core tools detailed usage instructions](core/USAGE.md)
* [core dependencies](core/DEPENDENCIES.md)
* [extended pipeline description](extended/README.md)
* [extended pipeline dependencies](extended/DEPENDENCIES.md)

## Code licensing
The taxator-tk source is licensend under the GPLv3 and builds on the following software components:

## Acknowledgments and contact

We acknowledge the work done by the contributors of the SeqAn and BOOST projects
as well as the author(s) of the tree.hh classes and all other anonymous
contributors. For the extended pipeline, we acknowledge all external programs.

For any questions, feedback or complaints, contact
science[at]fungs.de

Please, if you use this software and publish a paper, cite

    J. Dröge, I. Gregor, and A. C. McHardy
    Taxator-tk: precise taxonomic assignment of metagenomes by fast approximation of evolutionary neighborhoods
    Bioinformatics 2015 31: 817-824.
    doi: 10.1093/bioinformatics/btu745
