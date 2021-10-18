# taxator-tk
A set of programs for the taxonomic analysis of genetic sequences

## Introduction
This is the source code for the taxator toolkit. If you are looking for pre-packaged
programs or sample refpacks for taxonomic assignment, please go to the
[official download page](http://research.bifo.helmholtz-hzi.de/webapps/wa-download/).

This repository consists of two parts, the actual taxator-tk C++ code in the folder `core` and the pipeline scripts which are used to build the ready-to-use binary distribution in folder `extended`. The folder `contrib` contains information on how to add the external binaries which are part of the binary extended version package.

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

> 1: Help, the majority of assigned taxa are unspecific, like phylum or (super)kingdom, what can I do?

Your input sequences are likely too short or too novel, which means that they lack relative sequences in your refpack. Check that you are using a comprehensive, recent and suitable refpack for classification. For short sequences, you can still find the associated taxonomic information in the [intermediate GFF3](doc/fileformats.md) file, just not at the confidence level which is required for the final consensus assignment.

> 2: Why do taxator-tk pipelines run so long compared to other classification programs?

The classification requires both, sensitive local alignment to find distant homologs and approximate phylogenetic inference. Typically, alignment dominates the total runtime.

> 3: Can the sequences classification be accelerated?

There are many pipeline parameters that can be tuned, starting from the aligner itself, alignment length pre-filters, segment classification cutoff length or other parameters to reduce the number of pairwise alignments in the realignment placement algorithm. The provided pipelines try to provide a good balance and example.

> 4: Should taxator-tk be used with sequencing reads or assembled data?

The algorithm works for single short or long reads, but if your reads are overlapping, it does not make too much sense to apply it to unassembled data. This is, because the calculations for local alignment and inference would be repeated for overlapping regions, doing the same calculations over and over again, which leads to prohibitive runtimes. Secondly, the conservative consensus approach is most effective for longer sequences, either stemming from long read sequencing technologies or assembled contigs.

> 5: Is taxator-tk targeted towards a specific sequencing platform?

Look at FAQ#4. In general, for long-read sequences with specific error profiles such as PacBIO, you should probably use an aligner that considers the biases. Taxator-tk can be used for raw data, keeping in mind that the pairwise genetic distances will likely be over-estimated, leading to more conservative classifications.

> 6: If a sequence was assigned to a taxon by the algorithm, where can I see why?

First, look at the GFF3 output file which can also be loaded as a feature track into sequence viewers like IGV. Secondly, if enabled, the taxator logfile will give a very detailed description of the inference algorithm and the involved taxa. In protein space, it even contains the pairwise alignments.

> 7: Can I use taxator-tk with targeted sequencing data like amplicons?

You can, but to be honest, it probably does not make too much sense, because in the design of such a study, one typically already knows which gene to target. Using this information, you can do much better by phylogenetic inference and tree building. One of the main advantages of using taxator-tk is, that you neither need to know the genes analyzed nor you need to have curated and aligned reference genes for classification.

> 8: Why do I get an error with spaces in the FASTA identifier?

The FASTA format lacks a stringent format specification and different programs and local aligners interpret the identifiers differently. For instance, NCBI style adds spaces followed by meta-information. In order to work with all types of aligners, we follow the simplest implementation of the FASTA format and allow spaces in the identifiers but no extra meta-information. Use a program like [seqkit](https://github.com/shenwei356/seqkit) or the bundled `fasta-strip-identifier` to clean up sequence identifiers of the input sequences.

> 9: My computer is running out of main memory. How can I run taxator-tk on my desktop machine?

Taxator-tk was primarily designed for computers with large memory and many CPU cores in mind but we avoid holding data permanently in memory and use an on-disk FASTA index to access refpacks. If you use the extended pipelines with a low-memory aligner like Blast, it's mostly the (uncompressed) query sequences and taxonomy mapping which occupy space in main memory. So if your query sequences are particularly large, try using an on-disk FASTA index for the input sequences as well, or just split up the input sequences and run program taxator on the different pieces. You can safely concatenate the GFF3 output files for the last consensus binning step.

> 10: Can I run taxator-tk run on a compute cluster?

See FAQ#9. The core algorithm is embarassingly parallel, which means it can easily be split across compute nodes. However, you have to write your own pipeline for your compute environment, and you have to decide whether you want to do the local alignment as part of the pipeline or separately on the compute cluster.

> 11: What is the advantage/disadvanage in using taxator-tk over fast k-mer or hash based sequence classification with programs like Centrifuge or Kraken?

Algorithms using those techniques tend to assign most sequnces at species to family level. Due to phylogenetic neighborhood inference, taxator-tk handles distant homology in a more robust and less error-prone way. If this is not a priority, for instance for well sequenced environments and for human pathogens, these nearest-neighbor-style classifiers will most likely be a better fit. They are much faster.

> 12: What is the advantage/disadvantage in using taxator-tk over marker gene phylogenetic placement with programs like MetagPhlan, CheckM, RAxML or pplacer?

Taxator-tk frees you from selecting genes for analysis, identifiying them in the input sequences, building reference alignments, HMMs and phylogenies. It even handles gene-free input, as long as there is a phylogenetic signal in the data. Starting from v1.5 taxator-tk also supports protein sequences as input and the extended pipeline looks for protein coding regions in nucleotide input. On the other hand, phylogenetic placement is much more precise for specific genes and does not require a global taxonomy (which is often too rough or incorrect, compared to a gene tree).

> 13: What is a refpack?

A refpack is just a collection of reference sequences in FASTA format and a dump of the NCBI taxonomy in a specific version, which are both linked. Refpacks may contain nucleotide sequences of full genomes, partial genomes or genes or amino acid sequences for proteins.

> 14: Where can I find refpacks?

Ready to use refpacks can be downloaded along with the binary distribution of the program at the Helmholtz Centre for Infection Research. Because we cannot always keep up with new sequence collection releases or provide specialized refpacks, you are also encouraged to build your own. This is explained in the documentation or maybe you will find [this blog post](https://scienceblog.fungs.de/posts/taxator-tk-marine-refpack/) useful.

> 15: I want to build my own refpack, why do I have mapping problems with taxonomic identifiers?

The NCBI taxonomy changes frequently, replacing old identifiers by new ones. That means, old taxonomic mappings outdate quickly, when you update the taxonomy. You can always try to remap to the most recent taxonomy version, for instance using the web lookup at the [NCBI taxonomy home](https://www.ncbi.nlm.nih.gov/Taxonomy/).

> 16: Can I use a different taxonomy than NCBI?

No, currently we only support the NCBI dump file format in refpacks. In theory, you could of course bring your alternate taxonomy into the format using `names.dmp` and `nodes.dmp`.

> 17: Should I use nucleotide or protein alignment for classification?

That is a hard question. The protein alignment feature is pretty new in taxator-tk (starting at 1.5) and we don't want to give recommendations until we have more extensive experimental data. In theory, protein alignment should be more sensitive for finding distant homologs and more stable in distance estimation.

> 18: I want to use taxator-tk with aligner XYZ, can I?

Yes, just implement the [intermediate text-based and TAB-separated alignment format](doc/fileformats.md) and feed the output into `taxator`.

> 19: Can I use taxator-tk to quantify relative in-sample or between-sampe abundances (taxonomic profiling)?

In order to do so, you need to add sequence coverage information (for assembled sequences). There is no built-in functionality to do so. And you need to account for sequences in the input, that are not classified at all (see FAQ#20). Also note, that the depth of assignment is a function of the refpack taxonomic distribution and will be quite unbalanced.

> 20: Why are some of the input sequences not listed in the final results files?

The current pipelines are only based on the aligner output. If no alignment can be found, that input sequence is swallowed by the pipeline. This behavior can change in the future. The correct assignment for these sequences would be "unassigned" or taxonomy root (taxID 1).

> 21: I want to use a pipeline but it says that some program is missing. What can I do?

The pipelines combine taxator-tk core with other free software. Use the binary distribution or download and install the required programs before you run the pipeline. It is easiest to get pre-compiled software from the respective authors and place the binaries into the pipeline `bin` folder.

> 22: What is the meaning of the support value in output files?

The support value measures sequence similarity to the best reference, which supports the classification decision, in the unit number of positions. For nucleotide sequences, this is the estimated number of matches (lower bound). For amino acid sequence, it is a more abstract similarity based on the BLOSUM62 scoring matrix. In any case, dividing the sequence length by support positions yields a simple measure of percentage similarity.

> 23: I know that all my sequences come from the same organism. I want is a single classification, not one for every input sequence.

The extended pipeline scripts are made for typical metagenomic data which contain a mixture of different species. However, you can just rerun the program `binner`, which has a grouping command line parameter and is very fast, with the generated GFF3 file as input. Passing a catch-all regular expression `-g '(.*)'` will assume that all segments identified and classified in the previous calculation are from the same organism and the output will be a single classification.

## Sponsoring
<img alt="ermine" src="https://user-images.githubusercontent.com/8776981/137801205-87f0a851-bbba-4484-bd80-7ab6387cba78.png" width="100" height="100">

A license for Ermine was provided free of charge by [magicermine.com](https://magicermine.com) to make the Linux 64bit binary files for download portable.

## Code licensing
The taxator-tk source is licensend under the GPLv3 but builds on other free software components, which have their own respective licenses and are listed under dependencies.

## Acknowledgments and contact

We acknowledge the work done by the contributors of the SeqAn and BOOST projects
as well as the author(s) of the tree.hh classes and all other anonymous
contributors. For the extended pipeline, we acknowledge all the authors of all external programs, which would be too many to list them here.

For any questions, feedback or complaints, contact
science[at]fungs.de

Please, if you use this software and publish a paper, cite

    J. Dröge, I. Gregor, and A. C. McHardy
    Taxator-tk: precise taxonomic assignment of metagenomes by fast approximation of evolutionary neighborhoods
    Bioinformatics 2015 31: 817-824.
    doi: 10.1093/bioinformatics/btu745
