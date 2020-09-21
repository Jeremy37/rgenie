# rgenie

*rgenie* is an R package to analyze the sequencing output from a set of GenIE
experimental replicates.

## Installation

If you don't have devtools yet, first install devtools:
```
install.packages("devtools")
```
Then install rgenie using devtools:
```
devtools::install_github("jeremy37/rgenie")
```

## Overview

**GenIE** (genome-editing interrogation of enhancers) is an experimental method
to evaluate the effects of individual SNPs on gene transcription.

Briefly,
CRISPR-Cas9 is targeted near a SNP of interest, and typically an oligonucleotide
for homology-directed repair (HDR) will be included. With the editing done in
a pool of cells, the single nucleotide change can be achieved in a fraction of cell
chromosomes, and the remaining chromosomes will either be wild-type or will have
deletions in the region. You then extract both RNA and DNA from the pool of cells,
and do multiple replicates to amplify an amplicon from each of these (genomic DNA
and cDNA). rgenie can help to visualize the experiment and compute statistics for
the effects of different alleles on gene transcription.

The plot below shows the set of alleles measured in a set of gDNA and cDNA replicates.

![](https://github.com/Jeremy37/rgenie_example/raw/master/deletion_alleles_plot.png)

## Getting started

The [Introductory vignette](https://htmlpreview.github.io/?https://github.com/Jeremy37/rgenie/blob/master/vignettes/introduction.html) shows how to download example data and run an analysis.

The [rgenie in depth vignette](https://htmlpreview.github.io/?https://github.com/Jeremy37/rgenie/blob/master/vignettes/advanced_rgenie.html) provides more details on the parameters for some methods and how to interpret results.

## Citation

Cooper SE, Schwartzentruber J, Bello E, Coomber, EL, Bassett AR. Screening for functional transcriptional and splicing regulatory variants with GenIE. Not yet published.
