---
title: "Practical sessions part 2: scRNA"
author: "Tim Triche"
date: "October 27th, 2021"
output: 
  html_document:
    keep_md: true
vignette: >
  %\VignetteIndexEntry{Project2}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---



# Installation

Install the ExpDesign2021 package from github, if you haven't already.
If you're knitting this on your own machine, you will want to load knitr too.


```r
#install.packages("remotes")
#install.packages("BiocManager")
#BiocManager::install("VanAndelInstitute/ExpDesign2021")
library(knitr)
```

To extract just the R code, you can use knitr::knit(input, tangle=TRUE):




# Introduction

![representative scRNA plots](figure/scRNAseq.jpg)

_Extra credit assignment: figure out whose lab the above came from_

<details>
  <summary>Click for hint</summary>
  They're at the Gladstone Institute now, if memory serves. 
</details> 

For project 2 ("design your own damned experiment"), the class boldly rallied 
behind Richard Cassidy's suggestion of using single-cell RNA-seq data for an 
example dataset. This seemed like a great idea: it's big, it's messy, and it 
accompanies a [glam paper for a blue-collar experiment](https://www.nature.com/articles/s41587-020-0465-8), which is kind of cool. (The paper also covers 
single-cell vs. single-nucleus RNAseq comparisons, but we won't go there yet.)
The cell mixture model (HEK293 human cells and NIH3T3 mouse cells) is straight-
forward (mix two types of cells in two test tubes and then split).


# The Broad Institute barnyard data

_Question: why is it called a 'barnyard' experiment?_
<details>
  <summary>Click for Answer</summary>
  A traditional test of single-cell protocols is to mix cells from two species,
  and then see if the results can be easily classified into one or the other. 
  A protocol which does a good job of encapsulating cells will make it easier
  to classify each cell by species. Some 'barnyard' experiments have tossed in 
  canine or chicken cells for variety, but since it's such a common benchmark,
  people took to calling any species-mixing single-cell experiment a 'barnyard' 
  experiment, and the name stuck. 
</details> 

The library preparation methods involved in the paper may bear explanation. 
Short-read RNA sequencing experiments require cDNA as their input, since most
sequencers expect to call bases from DNA; as a consequence, one must first make
libraries of cDNA molecules from fragmented mRNA. Numerous approaches exist, 
six of which are employed here. 

![methods compared](figure/methods.png)

[The data is available from GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132044) although it's a bit unwieldy, weighing in at 300MB for the cell
mixture model counts alone. Happily, we can load up the metadata without that.

![bus format](figure/bustools.jpeg) 

The barcode, UMI, set or [BUS format](https://academic.oup.com/bioinformatics/article/35/21/4472/5487510) neatly encapsulates the process of transforming reads
resulting from a single-cell RNAseq experiment into counts of molecular barcodes
describing a cell, a molecule, and the likely genomic origin of that molecule.
It's one of the most compact representations possible for raw scRNAseq outputs, 
and the figure shows how it facilitates direct comparisons. The Broad used their
`scumi` tool instead, but the basic notion ("make everything comparable") and 
the resulting sparse matrix of counts is similar in spirit. 


## Cells


```r

# read in the cell names for the count matrix directly from GEO 
cells <- readLines(gzcon(url("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132044/suppl/GSE132044_mixture_hg19_mm10_cell.tsv.gz")))
# in English: "create a gzip connection from the URL at GEO, then read from it"

# tidy things up 
library(stringr)
library(tidyverse)

# now decode the metadata from it
tibble(name = cells) %>% # stringr::str_split takes strings and a split pattern
  mutate(experiment = str_split(name, "\\.", simplify=TRUE)[,1]) %>% # column 1
  mutate(method = str_split(name, "\\.", simplify=TRUE)[,2]) %>% # column 2 
  mutate(cell = str_split(name, "\\.", simplify=TRUE)[,3]) -> # column 3
    celltibble # the result is assigned to the object `celltibble`

# how many cells per mixture were run with each method?
with(celltibble, 
     table(method, experiment))
#>                  experiment
#> method            Mixture1 Mixture2
#>   10x-Chromium-v2     3159     3187
#>   CEL-Seq2             359      347
#>   Drop-seq            2594     3812
#>   Seq-Well            1627     1012
#>   Smart-seq2           342      343
#>   inDrops             3081     2529
#>   sci-RNA-seq          299     5023
```

Arguably, we don't need more than (say) 300 cells per method for this exercise,
and it would be nice not to demolish your computers' RAM (even though we're not
using Seurat or some piggish monstrosity like that). So we'll probably chop down
the number of cells from which you resample, at least for this project. If you
are familiar with block resampling, you already know where this is going. 

## Genes 

We can do the same type of thing for the genes involved, albeit without the 
intention of chopping them down by much if at all:


```r

# read in the gene names for the count matrix directly from GEO 
genes <- readLines(gzcon(url("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132044/suppl/GSE132044_mixture_hg19_mm10_gene.tsv.gz")))
# in English: "create a gzip connection from the URL at GEO, then read from it"

# tidy things up 
library(stringr)
library(tidyverse)

# now decode the metadata from it
tibble(name = genes) %>%  # stringr::str_split(string, pattern_to_split_on)
  mutate(ensembl = str_split(name, "_", simplify=TRUE)[,2]) %>% # ENS[MUS]G name
  mutate(genome = str_split(name, "_", simplify=TRUE)[,3]) -> # genome assembly
    genetibble # assign the result to `genetibble` 

# how many genes per assembly? 
with(genetibble, table(genome))
#> genome
#>  hg19  mm10 
#> 33354 28692
```

_Question: why didn't we use the gene symbols in columns 4:6?_
<details>
  <summary>Click for Answer</summary>
  Because some ENSEMBL gene names map to multiple RefGene/HUGO symbols. 
</details> 

With the release of Bioconductor 3.14, the project includes a tidy single cell 
experiment (data structure) package, which is great since all other single cell 
data structures kind of suck. (No, seriously, you'll find out why eventually.) 
The package is, not coincidentally, called [tidySingleCellExperiment](http://www.bioconductor.org/packages/release/bioc/vignettes/tidySingleCellExperiment/inst/doc/introduction.html):


```r

if (!require("tidySingleCellExperiment")) {
  BiocManager::install("SingleCellExperiment")
  library(SingleCellExperiment)
  BiocManager::install("tidySingleCellExperiment")
  library(tidySingleCellExperiment)
}

```


# Chopping down the data to fit in RAM

Above, I stated that maybe we don't need thousands of cells per method. 
(This is arguable; in exchange for crappy per-cell results, 10X and similar
offer more cells at the price of more dropouts and less sensitivity per cell.)
No worries, we'll just go ahead and downsample "enough" cells per method.
It turns out this happens to me often enough that I wrote some code to do it:

<details>
  <summary>Click for sample_umis() function code</summary>

```r

# adapted from a SingleCellExperiment-centric method for CITEseq
sample_umis <- function(umis, meta, block, ideal=300) {

  stopifnot(nrow(meta) == ncol(umis))
  stopifnot(length(block) == nrow(meta))

  pops <- sort(table(block))
  samplesets <- split(seq_len(nrow(meta)), block)

  keep <- integer()
  for (set in names(samplesets)) {
    sset <- samplesets[[set]]
    cells <- length(sset)
    if (cells <= ideal) {
      pct <- 100
      keep <- c(keep, sset)
      message("Kept ", cells, " cells (", pct, "%) of type ", set, ".")
    } else {
      kept <- sample(sset, size=ideal)
      pct <- round((ideal / cells) * 100)
      keep <- c(keep, kept)
      message("Kept ", ideal, " cells (", pct, "%) of type ", set, ".")
    }
  }

  pct <- round((length(keep) / ncol(umis)) * 100, 1)
  message("Kept ", length(keep), " (", pct, "%) of ", ncol(umis),
          " cells in ", length(samplesets), " blocks.")
  umis[, keep]

}
```
</details> 

_Question: can you generate block-random samples from the metadata you have?_
<details>
  <summary>Click for answer</summary>
  Easily, although you'd need to rewrite the above function to omit `umis`.
</details>


It's possible your machine or instance would crash if you did the following
(apparently nobody told IT that a laptop is a device to take RAM on a plane). 
I'm going to stick to the above "300 cells ought to be enough" and do this:

<details>
  <summary>Click for gory details</summary>

```r

# pre-downloaded MatrixMarket file from GEO
library(Matrix)

if (FALSE) {

  # put it into a column-sparse Matrix object to avoid wasting heaps of RAM
  umi_counts <- as(readMM("GSE132044_mixture_hg19_mm10_count_matrix.mtx.gz"), "dgCMatrix") # "read this as if it were already a column-sparse Matrix object"

  # never trust anyone, including yourself 
  stopifnot(nrow(umi_counts) == nrow(genetibble))
  stopifnot(ncol(umi_counts) == nrow(celltibble))

  # bolt the dimensions back onto the data
  rownames(umi_counts) <- genetibble$name
  colnames(umi_counts) <- celltibble$name

  # create a basis for blocked downsampling
  celltibble %>%
    mutate(block = paste(method, experiment, sep="_")) -> 
      blocktibble 

  # block downsample 
  downsampled <- sample_umis(umi_counts, blocktibble, blocktibble$block)
  # Kept 300 cells (9%) of type 10x-Chromium-v2_Mixture1.
  # Kept 300 cells (9%) of type 10x-Chromium-v2_Mixture2.
  # Kept 300 cells (84%) of type CEL-Seq2_Mixture1.
  # Kept 300 cells (86%) of type CEL-Seq2_Mixture2.
  # Kept 300 cells (12%) of type Drop-seq_Mixture1.
  # Kept 300 cells (8%) of type Drop-seq_Mixture2.
  # Kept 300 cells (10%) of type inDrops_Mixture1.
  # Kept 300 cells (12%) of type inDrops_Mixture2.
  # Kept 299 cells (100%) of type sci-RNA-seq_Mixture1.
  # Kept 300 cells (6%) of type sci-RNA-seq_Mixture2.
  # Kept 300 cells (18%) of type Seq-Well_Mixture1.
  # Kept 300 cells (30%) of type Seq-Well_Mixture2.
  # Kept 300 cells (88%) of type Smart-seq2_Mixture1.
  # Kept 300 cells (87%) of type Smart-seq2_Mixture2.
  # Kept 4199 (15.2%) of 27714 cells in 14 blocks.

  saveRDS(downsampled, file="barnyard_umi_sample.rds") 
  # online at https://ttriche.github.io/RDS/barnyard_umi_sample.rds

} 

```
</details>

Right then, let's get to work. You might consider running through some of the 
steps in the [tidySingleCellExperiment vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/tidySingleCellExperiment/inst/doc/introduction.html) to start with.  If your laptop or Rstudio Cloud instance still blows up, let 
me know and I'll create some smaller-er downsamples (also a familiar practice).


<details>
  <summary>Load the necessary packages</summary>


```r

# need to install from Bioconductor, 
# unless we're going to dork around with individual genes right from the start
library(SingleCellExperiment) 
library(tidySingleCellExperiment)
```
</details> 

How the experiment object is created, or (faster) loaded from a URL: 


```r

if (FALSE) { # if you want to recreate the SingleCellExperiment: 

  downsampled <- 
    readRDS(url("https://ttriche.github.io/RDS/barnyard_umi_sample.rds"))

  # create the column (cell) annotation data frame
  column_data <- as(celltibble, "DataFrame")
  rownames(column_data) <- column_data$name
  column_data <- column_data[colnames(downsampled), ]

  # create the row (gene) annotation data frame
  row_data <- as(genetibble, "DataFrame")
  rownames(row_data) <- row_data$name
  row_data <- row_data[rownames(downsampled), ]

  # create a SingleCellExperiment
  barnyard <- SingleCellExperiment(SimpleList(counts=downsampled),
                                   rowData=row_data, 
                                   colData=column_data)
  saveRDS(barnyard, file="barnyard.rds") 

}

# if you just want to load it 
tidybarnyard <- tidy(readRDS(url("https://ttriche.github.io/RDS/barnyard.rds")))
#> Loading required package: SingleCellExperiment
#> Loading required package: SummarizedExperiment
#> Loading required package: GenomicRanges
#> Loading required package: BiocGenerics
#> Error in .requirePackage(package): unable to find required package 'SingleCellExperiment'
show(tidybarnyard)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'show': object 'tidybarnyard' not found
```

I'd suggest, as a first pass, trying to sort out which are the mouse 3T3 cells,
and which are the human HEK293 cells. We can discuss this in class on Monday.

_Question: What's the easiest way to distinguish the mouse and human cells?_

<details>
  <summary>Hint:</summary>


```r

library(tidySingleCellExperiment) 
rowGenome <- rowData(tidybarnyard)$genome
UMIs <- counts(tidybarnyard) 

```
</details>

_Question: Can you think of a way to distinguish male from female cells?_

<details>
  <summary>Hint:</summary>


```r

# this requires a bit of domain knowledge, to be added

```
</details> 

_Question: Are the two answers above roughly equivalent?  Why?_


Bonus points if you can say where the two cell lines originally came from, and 
why the HEK cells don't have the usual PHI-non-compliant names from that time.

# Class exercises 

_Question: what is a UMI and why does it matter in single-cell RNAseq data?_

Load up either the UMI matrix or the SingleCellExperiment (tidy or untidy) and 
look for the `counts` assay.  (You may need to just use `counts(tidysce)`.) If 
your computer survives neither of these, you can try using summary statistics, 
but plan on getting an even smaller subsample to use for the actual Project. 

Or, just go back and look at the BUStools figure. 

_Question: what shall we investigate for project 2?_ 

I'm thinking the proportions of the mixtures, but I want to fit those first. 
Note that a t-test on proportions is not necessarily "proper"... ! 
