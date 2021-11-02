## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(devtools)
load_all("./")


## ---- loadpkgs, eval = FALSE, message = FALSE---------------------------------
## #install.packages("remotes")
## #install.packages("BiocManager")
## #BiocManager::install("VanAndelInstitute/ExpDesign2021")
## library(knitr)


## ---- tangle, eval = FALSE, message = FALSE, echo = FALSE---------------------
## # knitr::knit("project2.Rmd", tangle = TRUE)
## # [1] "project2.R"


## ---- barnyardcells, eval=TRUE------------------------------------------------

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



## ---- barnyardgenes, eval=TRUE------------------------------------------------

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



## ---- tidysinglecell, eval=FALSE----------------------------------------------
## 
## if (!require("tidySingleCellExperiment")) {
##   BiocManager::install("tidySingleCellExperiment")
##   library(tidySingleCellExperiment)
## }
## 


## ---- sample_umis-------------------------------------------------------------

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



## ---- downsampling, eval=FALSE------------------------------------------------
## 
## # pre-downloaded MatrixMarket file from GEO
## library(Matrix)
## 
## if (FALSE) {
## 
##   # put it into a column-sparse Matrix object to avoid wasting heaps of RAM
##   umi_counts <- as(readMM("GSE132044_mixture_hg19_mm10_count_matrix.mtx.gz"), "dgCMatrix") # "read this as if it were already a column-sparse Matrix object"
## 
##   # never trust anyone, including yourself
##   stopifnot(nrow(umi_counts) == nrow(genetibble))
##   stopifnot(ncol(umi_counts) == nrow(celltibble))
## 
##   # bolt the dimensions back onto the data
##   rownames(umi_counts) <- genetibble$name
##   colnames(umi_counts) <- celltibble$name
## 
##   # create a basis for blocked downsampling
##   celltibble %>%
##     mutate(block = paste(method, experiment, sep="_")) ->
##       blocktibble
## 
##   # block downsample
##   downsampled <- sample_umis(umi_counts, blocktibble, blocktibble$block)
##   # Kept 300 cells (9%) of type 10x-Chromium-v2_Mixture1.
##   # Kept 300 cells (9%) of type 10x-Chromium-v2_Mixture2.
##   # Kept 300 cells (84%) of type CEL-Seq2_Mixture1.
##   # Kept 300 cells (86%) of type CEL-Seq2_Mixture2.
##   # Kept 300 cells (12%) of type Drop-seq_Mixture1.
##   # Kept 300 cells (8%) of type Drop-seq_Mixture2.
##   # Kept 300 cells (10%) of type inDrops_Mixture1.
##   # Kept 300 cells (12%) of type inDrops_Mixture2.
##   # Kept 299 cells (100%) of type sci-RNA-seq_Mixture1.
##   # Kept 300 cells (6%) of type sci-RNA-seq_Mixture2.
##   # Kept 300 cells (18%) of type Seq-Well_Mixture1.
##   # Kept 300 cells (30%) of type Seq-Well_Mixture2.
##   # Kept 300 cells (88%) of type Smart-seq2_Mixture1.
##   # Kept 300 cells (87%) of type Smart-seq2_Mixture2.
##   # Kept 4199 (15.2%) of 27714 cells in 14 blocks.
## 
##   saveRDS(downsampled, file="barnyard_umi_sample.rds")
##   # online at https://ttriche.github.io/RDS/barnyard_umi_sample.rds
## 
## }
## 


## ---- tidysce, eval=FALSE-----------------------------------------------------
## 
## # need to install from Bioconductor,
## # unless we're going to dork around with individual genes right from the start
## library(SingleCellExperiment)
## library(tidySingleCellExperiment)


## ---- createtidysce, eval=TRUE------------------------------------------------

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
show(tidybarnyard)



## ---- bygenome, eval=FALSE----------------------------------------------------
## 
## library(tidySingleCellExperiment)
## rowGenome <- rowData(tidybarnyard)$genome
## UMIs <- counts(tidybarnyard)
## 


## ---- bysexchr, eval=FALSE----------------------------------------------------
## 
## # this requires a bit of domain knowledge, to be added
## 

