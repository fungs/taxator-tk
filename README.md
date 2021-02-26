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
* [core compilation guide](core/BUILD.md)
* [core detailed usage instructions](core/USAGE.md)
* [core dependencies](core/DEPENDENCIES.md)
* [extended pipeline description](extended/README.md)
* [extended pipeline dependencies](extended/DEPENDENCIES.md)

## FAQ

### Help, the majority of assigned taxa are unspecific, like phylum or (super)kingdom, what can I do?

### Why do taxator-tk pipelines run so long compared to others?

### Can the sequences classification be accelerated?

### Is taxator-tk made for reads or assembled sequences? Or is the program targeted towards a specific sequencing platform?

### If a sequence was assigned to a taxon by the algorithm, where can I see why?

### Can I use taxator-tk with targeted sequencing data like amplicons?

### Why do I get an error if there are spaces in the FASTA identifiers?

### I'm running out of memory. Can I run taxator-tk on my desktop machine?

### What is the advantage/disadvanage in using taxator-tk over fast k-mer or hash based sequence classification with programs like Centrifuge or Kraken?

### What is the advantage/disadvantage in using taxator-tk over marker gene phylogenetic placement with programs like MetagPhlan, CheckM, RAxML or pplacer?

### What is a refpack?

### Where can I find refpacks?

### I want to build my own refpack, why do I have mapping problems with taxonomic identifiers?

### Can I use a different taxonomy than NCBI?

### Should I use a nucleotide or protein alignment and classification pipeline?

### I want to use taxator-tk with aligner XYZ, can I?

### Can I use taxator-tk to quantify relative in-sample or between-sampe abundances?

### Why are some of the input sequences not listed in the final results files?

### I want to use a pipeline but it says that some program is missing. What can I do?

### What does the support value in the output mean?



## Code licensing
The taxator-tk source is licensend under the GPLv3 but builds on other free software components, which have their own respective licenses and are listed under dependencies.

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
