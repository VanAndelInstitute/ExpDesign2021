## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(devtools)
load_all("./")


## ---- tidysinglecell----------------------------------------------------------

if (!require("SingleCellExperiment")) {
  BiocManager::install("SingleCellExperiment")
  library(SingleCellExperiment)
}
if (!require("tidySingleCellExperiment")) {
  BiocManager::install("tidySingleCellExperiment")
  library(tidySingleCellExperiment)
}



## ---- loadpackages------------------------------------------------------------
library(SingleCellExperiment) 
library(tidySingleCellExperiment)


## ---- gettidysce--------------------------------------------------------------
# should package this... could just instantiate from a package via data(...)
if (!exists("tidybarnyard")) { 
  tidybarnyard <- readRDS(url("https://ttriche.github.io/RDS/tidybarnyard.rds"))
}


## ---- rowDataAndColumnData----------------------------------------------------

# as what is our object masquerading?
show(tidybarnyard[,0]) # "just show me information about it, with 0 cells"

# how many rows (genes) and columns (cells) are there in our tidy barnyard?
dim(tidybarnyard)

# is this equivalent to what we get from nrow() and ncol()? 
identical(dim(tidybarnyard), 
          c(nrow(tidybarnyard), ncol(tidybarnyard)))

# the assay(s) should also have this many rows and columns. Do they?
identical(dim(tidybarnyard), dim(assay(tidybarnyard)))

# there is a shortcut for the `counts` assay, since it's so common:
counts(tidybarnyard)[1:3, 1:3]

# does each row (gene) in the object have a corresponding rowData() row?
identical(rownames(tidybarnyard), rownames(rowData(tidybarnyard)))

# does each column (cell) in the object have a corresponding colData() row?
identical(colnames(tidybarnyard), rownames(colData(tidybarnyard)))

# is there a rowData row and a colData row that corresponds to each UMI count?
identical(dim(counts(tidybarnyard)), 
          c(nrow(rowData(tidybarnyard)), nrow(colData(tidybarnyard)))) 

# does that mean we can ask for a random cell with certain attributes?
as_tibble(tidybarnyard) %>%           # "turn the colData into a tibble"
  filter(method == "inDrops") %>%     # "choose rows where method == inDrops"
  slice_sample %>%                    # "randomly slice out one row"
  pull("cell") -> aCell               # "pull the column cell, assign to aCell"

# your cell:
show(aCell) 

# does that mean we can ask for a random gene with certain attributes?
as_tibble(rowData(tidybarnyard)) %>%  # "turn the rowData into a tibble"
  filter(genome == "mm10") %>%        # "choose rows where genome == mm10"
  slice_sample %>%                    # "randomly slice out one row"
  pull("name") -> aGene               # "pull the column name, assign to aGene" 

# your gene:
show(aGene) 

# how many copies of this random gene were found in this random cell? 
counts(tidybarnyard)[aGene, aCell] 

# note that the odds are good that you'll get a 0 for this random combination:
library(Matrix)                                 # for the `nnzero` function 
sparsity <- 1 - (nnzero(counts(tidybarnyard)) / # number of nonzero counts
                 prod(dim(tidybarnyard)))       # number of rows * columns 

# specifically, the chance of getting a zero is about 92.3%: 
message(round(sparsity * 100, 1), "%")

# This is fairly typical for a single-cell experiment, perhaps a bit low even.
# You could also take 10000 or so random samples to estimate the sparsity.
# For example, you could grab 100 genes from 100 cells at a time: 
as_tibble(tidybarnyard) %>%       # turn tidybarnyard's colData into a tibble
  slice_sample(n=100) %>%         # slice out 100 random rows, and then... 
  pull("cell") -> aHundredCells   # assign the "cell" column to "aHundredCells"

as_tibble(rowData(tidybarnyard)) %>% # turn tidybarnyard's rowData into a tibble
  slice_sample(n=100) %>%            # slice out 100 random rows, and then 
  pull("name") -> aHundredGenes      # assign "name" column to "aHundredGenes"

samples <- (length(aHundredCells) * length(aHundredGenes))
nonzero <- nnzero(counts(tidybarnyard)[aHundredGenes, aHundredCells])
sparsity_hat <- (samples - nonzero) / samples 
sparsity_hat # estimated sparsity

# In fact, we can use this scheme to look at sampling error:
sample_sparsity <- function(object, cells=100, genes=100) { 

  samples <- cells * genes 
  columns <- pull(slice_sample(as_tibble(object), n=cells), "cell")
  rows <- pull(slice_sample(as_tibble(rowData(object)), n=genes), "name")
  nonzero <- nnzero(counts(object)[rows, columns])
  sparsity <- (samples - nonzero) / samples
  return(sparsity) 

}

# the `replicate` function allows us to apply this many times over:
estimates <- replicate(n=100, sample_sparsity(tidybarnyard)) # "do it 100 times"

# plot estimates: 
library(ggplot2)
ggplot(tibble(estimate=estimates), aes(estimate)) + 
  geom_histogram() + 
  geom_vline(xintercept=sparsity, color="red", lwd=3) +
  theme_minimal() + 
  ggtitle("Sparsity of UMI matrix (true value in red)")

# the Central Limit Theorem lives to fight another day,
# and we have a decent idea of how to navigate our data.


## ---- UMAP, eval=FALSE--------------------------------------------------------
## 
## library(scater)
## # it is standard to log-normalize counts
## # (although it's not actually a great idea)
## tidybarnyard <- logNormCounts(tidybarnyard)
## 
## # compute UMAP embedding on the most variable genes
## tidybarnyard %>% runUMAP(ncomponents=3) -> tidyUMAP
## 
## # plot using plotly, color by method
## tidyUMAP %>%
##   plot_ly(
##     x =~ `UMAP1`,
##     y =~ `UMAP2`,
##     z =~ `UMAP3`,
##     color =~ method
##   )
## 


## ---- bygenome----------------------------------------------------------------

# identify the mouse genes:
as_tibble(rowData(tidybarnyard)) %>%  # make a tibble from the rowData
  filter(genome == "mm10") %>%        # select just the mouse mm10 genes 
  pull("name") -> mouseGenes          # assign their name column to mouseGenes

# for technical reasons, it's faster to tally expressed genes this way:
tidybarnyard$fracmouse <- 
  (colSums(counts(tidybarnyard)[mouseGenes, ] > 0)) / length(mouseGenes)

# do the same thing but with human genes:
as_tibble(rowData(tidybarnyard)) %>%  # make a tibble from the rowData
  filter(genome == "hg19") %>%        # select just the human hg19 genes 
  pull("name") -> humanGenes          # assign their name column to humanGenes

# same remarks as previously
tidybarnyard$frachuman <- 
  (colSums(counts(tidybarnyard)[humanGenes, ] > 0)) / length(humanGenes)

# why normalize by gene count?
as_tibble(rowData(tidybarnyard)) %>% select("genome") %>% table


# let's use ggplot to make sense of the results:
barnyardtibble <- as_tibble(tidybarnyard)

p <- ggplot(barnyardtibble, aes(x=fracmouse, y=frachuman)) + 
  xlab("Mouse transcripts expressed") +
  scale_x_continuous(labels = scales::percent) +
  ylab("Human transcripts expressed") +
  scale_y_continuous(labels = scales::percent) +
  geom_point(alpha=0.75, color="lightblue") +
  geom_density2d(alpha=0.5, color="blue") + 
  theme_minimal() 

# first pass at a plot: 
p + ggtitle("Barnyard experiment") 



## ---- gating------------------------------------------------------------------

# arbitrarily:
minhuman <- 0.05
minmouse <- 0.05 
# Let's label any cell that is more than 5% (0.05) of BOTH genomes as suspect.

# add a label:
barnyardtibble %>% 
  mutate(label = 
         case_when(frachuman > minhuman & fracmouse < minmouse ~ "human", 
                   fracmouse > minmouse & frachuman < minhuman ~ "mouse",
                   TRUE ~ "suspect")) -> barnyardtibble

# add our initial stab at labeling:
p <- ggplot(barnyardtibble, aes(x=fracmouse, y=frachuman, color=label)) + 
  xlab("Mouse transcripts expressed") +
  scale_x_continuous(labels = scales::percent) +
  ylab("Human transcripts expressed") +
  scale_y_continuous(labels = scales::percent) +
  geom_point(alpha=0.75) + 
  theme_minimal() 

# plot it with some gates and density lines:
p + geom_density2d(alpha=0.5, color="blue") + 
  geom_vline(xintercept=minmouse, color="red") + 
  geom_hline(yintercept=minhuman, color="red") + 
  ggtitle("First stab at gating")



## ---- doubletgate-------------------------------------------------------------

# arbitrarily:
maxhuman <- 0.21
maxmouse <- 0.19 

# make it obvious which cells are going to be gated out if we do this: 
barnyardtibble %>% 
  mutate(doublet = frachuman > maxhuman | fracmouse > maxmouse) %>%
  mutate(shading = case_when(label == "suspect" ~ 0.3,
                             doublet == TRUE ~ 0.1, 
                             TRUE ~ 0.5)) -> barnyardtibble

# start like the previous plot, but add low and high "gates" for each species:
p <- ggplot(barnyardtibble, 
            aes(x=fracmouse, y=frachuman, color=label, alpha=I(shading))) + 
  xlab("Mouse transcripts expressed") +
  scale_x_continuous(labels = scales::percent) +
  ylab("Human transcripts expressed") +
  scale_y_continuous(labels = scales::percent) +
  geom_point() + 
  geom_segment(y=0, yend=minhuman, x=minmouse, xend=minmouse, color="black") + 
  geom_segment(y=0, yend=minhuman, x=maxmouse, xend=maxmouse, color="black") + 
  geom_segment(y=minhuman, yend=minhuman, x=minmouse, xend=maxmouse, 
               color="black") + 
  geom_segment(y=minhuman, yend=minhuman, x=0, xend=minmouse, color="black") + 
  geom_segment(y=maxhuman, yend=maxhuman, x=0, xend=minmouse, color="black") + 
  geom_segment(y=minhuman, yend=maxhuman, x=minmouse, xend=minmouse, 
               color="black") + 
  theme_minimal() 
 
# plot it
p + ggtitle("aggregate doublet gating")




## ---- byMethod----------------------------------------------------------------

# reuse the plot again:
p + facet_wrap(~ method)



## ---- regressFactors----------------------------------------------------------

# we need a 0/1 outcome to perform logistic regression: 
barnyardtibble %>% mutate(classifiable = label != "suspect") -> barnyardtibble

# logistic regression in R uses the glm() or general linear model function,
# with a binomial (0/1) link, and the result can be compared like any lm(): 
fit0 <- glm(classifiable ~ 1, data=barnyardtibble, family=binomial) # null model

# add `method` to the predictors 
fit1 <- update(fit0, classifiable ~ method)

# add `method` and `experiment` to the predictors 
fit2 <- update(fit1, classifiable ~ method + experiment) 

# add `method` and `experiment` to the predictors, no intercept
fit3 <- update(fit2, classifiable ~ method + experiment + 0) 

# add `method` as the sole predictor, no intercept
fit4 <- update(fit0, classifiable ~ method + 0) 



## ---- coefficients------------------------------------------------------------

# it's better to use confidence intervals than p-values for fitting purposes,
# and it's even better yet to plot them all:
library(coefplot)

# classifiable ~ method
coefplot(fit1, trans=invlogit) + theme_minimal() 

# classifiable ~ method + experiment
coefplot(fit2, trans=invlogit) + theme_minimal() 

# classifiable ~ method + experiment, no intercept
coefplot(fit3, trans=invlogit) + theme_minimal() 

# classifiable ~ method, no intercept
coefplot(fit4, trans=invlogit) + theme_minimal() 



## ---- mixtureModel------------------------------------------------------------

# one of the greatest software packages ever written, 
# fits a Gaussian mixture model with arbitrary covariance structure and uses 
# a Bayesian penalization scheme to choose how many components exist in the mix
library(mclust)

# since these are proportional values, it makes sense to transform them: 
mfit <- Mclust(logit(barnyardtibble[, c("fracmouse","frachuman")]), 
               verbose=FALSE, G=1:3) # verbose=FALSE to avoid progress bar!
# here I have restricted Mclust to fitting, at most, 3 components (G=1:3). 
# this speeds up the process and avoids some difficult questions later on ;-)

# create a new column of the barnyard tibble with the results: 
table(mfit$classification) # it turns out that we end up with less human cells
barnyardtibble$mclass <- factor(mfit$classification)

# add mixture assignments:
p <- ggplot(barnyardtibble, 
            aes(x=fracmouse, y=frachuman, color=mclass, shape=label)) +
  xlab("Mouse transcripts expressed") +
  scale_x_continuous(labels = scales::percent) +
  ylab("Human transcripts expressed") +
  scale_y_continuous(labels = scales::percent) +
  geom_point(alpha=0.5) + 
  theme_minimal() 
 
# plot it
p + ggtitle("mixture model fit")



## ---- remix-------------------------------------------------------------------

# relabel the mixture class calls based on the previous plot
barnyardtibble %>% mutate(mclass = case_when(mclass == 1 ~ "mouse", 
                                             mclass == 3 ~ "human", 
                                             TRUE ~ "suspect"))-> barnyardtibble
barnyardtibble %>% mutate(mclassifiable = mclass != "suspect") -> barnyardtibble

# null model 
fitm0 <- glm(mclassifiable ~ 0, data=barnyardtibble, family=binomial) # random

# regress `mclassifiable` on method, no intercept:
fitm1 <- update(fitm0, mclassifiable ~ method + 0) 
coefplot(fitm1, trans=invlogit) + theme_minimal() 

# regress `mclassifiable` on method and experiment, no intercept:
fitm2 <- update(fitm0, mclassifiable ~ method + experiment + 0) 
coefplot(fitm2, trans=invlogit) + theme_minimal() 

# regress `mclassifiable` on method interacting with experiment, no intercept:
fitm3 <- update(fitm0, mclassifiable ~ method * experiment + 0) 
coefplot(fitm3, trans=invlogit) + theme_minimal() 



## ---- sample_umis-------------------------------------------------------------

# adapted from a SingleCellExperiment-centric method for CITEseq
sample_umis <- function(umis, meta, block, ideal=300) {

  # {{{
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
  # }}}

}

