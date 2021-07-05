# rgenie

*rgenie* is an R package to analyze the sequencing output from a set of GenIE
experimental replicates. It can be used with either ATAC or RNA GenIE experiments.

## Installation

rgenie is available on CRAN.
```
install.packages("rgenie")
```

If want the latest development version, first install devtools, and then install rgenie using devtools:
```
install.packages("devtools")
devtools::install_github("jeremy37/rgenie")
```

## Overview

**GenIE** (genome-editing interrogation of enhancers) is an experimental method
to evaluate the effects of individual SNPs on either gene transcription (when
RNA is the readout) or on chromatin accessibility (when ATAC is the readout).

Briefly, CRISPR-Cas9 is targeted near a SNP of interest, and typically an
oligonucleotide for homology-directed repair (HDR) will be included. With the
editing done in a pool of cells, the single nucleotide change can be achieved in
a fraction of cell chromosomes, and the remaining chromosomes will either be
wild-type or will have deletions in the region. For an RNA readout, you then extract both RNA and
DNA from the pool of cells, and do multiple replicates to amplify an amplicon
from each of these (genomic DNA and cDNA). rgenie can help to visualize the
experiment and compute statistics for the effects of different alleles on gene
transcription. For an ATAC readout the workflow is similar, except that ATAC is
done instead of RNA extraction, and there are some differences in the
amplification protocol, as described in the papers cited below.

This plot shows the set of alleles measured in a set of gDNA and cDNA replicates.

![](https://github.com/Jeremy37/rgenie/raw/master/example_data/deletion_alleles_plot.png)

## Getting started

The [rgenie Introduction vignette](https://htmlpreview.github.io/?https://github.com/Jeremy37/rgenie/blob/master/vignettes/introduction.html) shows how to download example data and run an analysis for RNA.

The [rgenie in depth vignette](https://htmlpreview.github.io/?https://github.com/Jeremy37/rgenie/blob/master/vignettes/advanced_rgenie.html) provides more details on the parameters for some methods and how to interpret results.

The [rgenie for ATAC vignette](https://htmlpreview.github.io/?https://github.com/Jeremy37/rgenie/blob/master/vignettes/genie_atac.html) shows how to run an analysis for a GenIE-ATAC experiment.

## Citation

rgenie for ATAC:
Sarah E Cooper, Jeremy Schwartzentruber, Eve L Coomber, Andrew R Bassett, [Identification of functional regulatory variants in open chromatin using GenIE-ATAC](), *in preparation*.

rgenie:
Sarah E Cooper, Jeremy Schwartzentruber, Erica Bello, Eve L Coomber, Andrew R Bassett, [Screening for functional transcriptional and splicing regulatory variants with GenIE](https://doi.org/10.1093/nar/gkaa960), *Nucleic Acids Research*, Volume 48, Issue 22, 16 December 2020.

## Contact

Feel free to submit an issue on github, or contact me:

jeremy37 at gmail dot com
