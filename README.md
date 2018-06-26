clustree paper
=================

Data and analysis for the clustree paper

Directory structure
--------------------

* `clustree_paper.Rmd` - RMarkdown file containing the manuscript
* `data` - Data used in the analysis in the paper
* `figures` - Figures used in the main paper
* `R` - R functions used in analysis
* `style` - Document style files

R
---

This directory contains the following functions used in the analysis:

* `simulations.R` - Functions for creating, clustering and plotting simulated
  datasets

Data
----

* `barcodes.tsv` - Cell barcode information for the PBMC dataset
* `genes.tsv` - Gene information for the PBMC dataset
* `matrix.mtx` - Gene expression matrix for the PBMC dataset
* `simulations.Rds` - Simulated datasets

Figures
-------

* `figure2_algorithm.png` - Diagram of the clustering tree algorithm
* `figure2_algorithm.pdf` - Diagram of the clustering tree algorithm
* Automatically generated figures will also be placed here

Style
-----

* `gigascience.csl` - Bibliography style file
* `references.bib` - Bibliography database
* `style.docx` - Reference DOCX style file
* `style.Rmd` - Reference Rmd style file
* `style_refs.bib` - Bibliography for the reference style
