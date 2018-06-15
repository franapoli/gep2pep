---
output: github_document
---

<!-- Grab your social icons from https://github.com/carlsednaoui/gitsocial -->
[1.2]: http://i.imgur.com/wWzX9uB.png (me on Twitter)
[1]: http://www.twitter.com/franapoli
<!-- Grab your social icons from https://github.com/carlsednaoui/gitsocial -->

## gep2pep

Pathway Expression Profiles (PEPs) are based on the expression of
pathways (defined as sets of genes) as opposed to individual
genes. This package converts gene expression profiles to PEPs and
performs enrichment analysis of both pathways and experimental
conditions, such as "Drug Set Enrichment Analysis" (finding pathways
that are consistently dysregulated by a set of drugs) and "gene2drug"
analysis (finding drugs that dysregulate a set of pathways or a single
gene).

Two papers have been published in Bioinformatics covering gep2pep
methods:

- A [paper](http://rdcu.be/pklt) about Drug Set Enrichment Analysis.
- A
  [paper](https://academic.oup.com/bioinformatics/article/34/9/1498/4721786)
  about gene2drug analysis.

Two corresponding webtools are online, which use Cmap data for both
types of analysis:

- [dsea.tigem.it](http://dsea.tigem.it)
- [gene2drug.tigem.it](http://gene2drug.tigem.it)

Gep2pep is maintained by Francesco Napolitano [![alt text][1.2]][1]


## Download and Installation

The latest stable release can be downloaded from Bioconductor at
[https://bioconductor.org/packages/release/bioc/html/gep2pep.html](https://bioconductor.org/packages/release/bioc/html/gep2pep.html). Further
instructions ar provided there.

Latest non-stable versions can be found on Github at
[https://github.com/franapoli/gep2pep](https://www.github.com/franapoli/gep2pep),
downloaded and then installed as follows:

    > install.packages("path-to-downloaded-source", repos=NULL)

## News

### v1.1.1.1

In progress version.

- Added "SGE mode", including the function `addSingleGeneSets`. This
  is to support fast computation of Kolmogorov-Smirnov statistics for
  large collections of gene sets including a single gene, which is
  useful to support gene-centric (as opposed to pathway-centric)
  analysis.
  
- `gene2pathways` now accepts a list of genes and returns all the
  pathways including either ALL of them or ANY of them according to
  the new `and` paramater.
  

### v1.1.1

- added support to deal with MsigDB release 6.1, which contains
  unconventional set categories ("ARCHIVED")

- added raw-mode to deal with large datasets. Raw mode stores PEPs to
  separate files during conversion, thus can be easily parallelized
  
- added "organism" parameter to `importMSigDB` to select sets

- added hdf5 support for large collections of PEPs

### v1

- first release
