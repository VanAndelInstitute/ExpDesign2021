## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(devtools)
#load_all("./")


## ---- tidysinglecell----------------------------------------------------------

if (!require("SingleCellExperiment")) {
  BiocManager::install("SingleCellExperiment")
  library(SingleCellExperiment)
}
if (!require("tidySingleCellExperiment")) {
  BiocManager::install("tidySingleCellExperiment")
  library(tidySingleCellExperiment)
}
if (!require("tidyverse")) {
  BiocManager::install("tidyverse")
  library(tidyverse)
}



## ---- loadpackages------------------------------------------------------------
library(SingleCellExperiment) 
library(tidySingleCellExperiment)


## ---- gettidysce--------------------------------------------------------------
# should package this... could just instantiate from a package via data(...)
if (!exists("tidybarnyard")) { 
  tidybarnyard <- readRDS(url("https://ttriche.github.io/RDS/tidybarnyard.rds"))
}


## ---- tidyfix-----------------------------------------------------------------

# if any column in `tidybarnyard` column data is named `cell`, rename it 
# this is currently a bug in tidySingleCellExperiment:
# https://github.com/stemangiola/tidySingleCellExperiment/issues/38
#
names(colData(tidybarnyard)) <- sub("cell",                       # pattern
                                    "barcode",                    # replacement
                                    names(colData(tidybarnyard))) # strings 


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
counts(tidybarnyard)[1:3, 1:3] # this holds the UMI counts, or a `.` for 0.

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
counts(tidybarnyard)[aGene, aCell]    # UMI counts for a given [gene, cell]. 

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



## ---- loadmclust--------------------------------------------------------------

# one of the greatest software packages ever written, 
# fits a Gaussian mixture model with arbitrary covariance structure and uses 
# a Bayesian penalization scheme to choose how many components exist in the mix
library(mclust)



## ---- mixtureModel------------------------------------------------------------

# `logit` is from the `gtools` package:
library(gtools) 
# Note to self: always compile vignettes in a fresh session with R --vanilla :-/

# since these are proportional values, it makes sense to transform them: 
mfit <- Mclust(logit(barnyardtibble[, c("fracmouse","frachuman")]), 
               verbose=FALSE, G=1:3) # verbose=FALSE to avoid progress bar!
# here I have restricted Mclust to fitting, at most, 3 components (G=1:3). 
# this speeds up the process and avoids some difficult questions later on ;-)

# create a new column of the barnyard tibble with the results: 
table(mfit$classification) # it turns out that we end up with less human cells
barnyardtibble$mclass <- factor(mfit$classification)

# plot the results
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



## ---- mixlabels---------------------------------------------------------------

# confusion matrix helps us assign correspondence
tbl <- with(barnyardtibble, table(mclass, label))
mouseclass <- which.max(tbl[, "mouse"])
humanclass <- which.max(tbl[, "human"])

# relabel the mixture assignments: 
barnyardtibble %>% 
  mutate(mixlabel = case_when(mclass == mouseclass ~ "mouse", 
                              mclass == humanclass ~ "human", 
                              TRUE ~ "suspect")
         ) -> barnyardtibble

# how did we do? 
with(barnyardtibble, table(mixlabel, label))

# specifically, do we label all the human and mouse cells confidently?
with(barnyardtibble, table(mixlabel, label))[, c("human", "mouse")]


## ---- mixlabeledplot----------------------------------------------------------

# add mixture labels to the plot:
p <- ggplot(barnyardtibble, 
            aes(x=fracmouse, y=frachuman, color=mixlabel, shape=label)) +
  xlab("Mouse transcripts expressed") +
  scale_x_continuous(labels = scales::percent) +
  ylab("Human transcripts expressed") +
  scale_y_continuous(labels = scales::percent) +
  geom_point(alpha=0.5) + 
  theme_minimal() 
 
# plot it
p + ggtitle("mixture model fit with labels")



## ---- remix-------------------------------------------------------------------

# classifiable _by mixture model_ 
barnyardtibble %>% 
  mutate(mclassifiable = mixlabel != "suspect") -> barnyardtibble

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
sample_umis <- function(umis, meta, block, ideal=300, verbose=TRUE, justnames=FALSE) {

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
      if (verbose) message("Kept ",cells," cells (",pct,"%) of type ",set,".")
    } else {
      kept <- sample(sset, size=ideal)
      pct <- round((ideal / cells) * 100)
      keep <- c(keep, kept)
      if (verbose) message("Kept ",ideal," cells (",pct,"%) of type ",set,".")
    }
  }

  pct <- round((length(keep) / ncol(umis)) * 100, 1)
  if (verbose) {
    message("Kept ", length(keep), " (", pct, "%) of ", ncol(umis),
            " cells in ", length(samplesets), " blocks.")
  }

  if (justnames) return(keep)
  else return(umis[, keep])

}


## ---- resampletidysce---------------------------------------------------------

# simple classifier: is it suspect? then it's a failure
performance_by_method <- function(subsample, atibble) {
 
  subs <- mutate(slice(atibble, subsample), # i.e., tbl %>% slice %>% mutate
                 result=if_else(mixlabel == "suspect", "failure", "success"))
  tbl <- table(select(subs, method, result)) # force of habit
  probs <- as(sweep(tbl, 1, rowSums(tbl), `/`), "matrix")
  rownames(probs) <- NULL # avoid hassles when stacking
  tibble(cbind(method=rownames(tbl), data.frame(probs)))

}

# a simple bootstrap (for efficiency, we don't bother grabbing the UMIs at all)
# per the above (and to make computing intervals easier), we'll use 200x200 
subsamples <- replicate(n=200, simplify=FALSE, # to make it easier to iterate
                        sample_umis(tidybarnyard, as_tibble(tidybarnyard), 
                                    tidybarnyard$method, ideal=200, 
                                    justnames=TRUE, verbose=FALSE))
names(subsamples) <- paste0("run", seq_along(subsamples))

# purrr::map maps a function over a list
library(purrr)

# evaluate performance with 200 bootstrap samples of (ideally) 200 cells/method
# (i.e. subsamples <- replicate(n=200, sample_umis(..., ideal=200, ...)), above)
runs <- purrr::map(subsamples, performance_by_method, atibble=barnyardtibble)

# for plotting, stack them row-by-row
results <- bind_rows(runs) 



## ---- ordered_sinaplot--------------------------------------------------------

# fire up geom_sina because it rules:
library(ggforce)

# for better resolution, you could resample batches of 500, 1000, ... 
results %>% 
  mutate(method = fct_reorder(method, desc(success), .fun='median')) %>% 
  ggplot(aes(method, success, color=method)) + 
            scale_y_continuous(labels = scales::label_percent()) + 
            ylab("classification performance") +
            geom_sina(show.legend=FALSE) + 
            coord_flip() + 
            theme_minimal() + 
            ggtitle("Cell classification results by prep, 200 cells apiece")



## ---- confidenceintervals-----------------------------------------------------

# ground rules
CI <- 0.95                                      # must be between 0 and 1 
lower <- (1 - CI) / 2                           # defined by the value of CI
upper <- 1 - lower                              # defined by the value of lower
middle <- (lower + upper) / 2                   # this is a tautology and a test
plot_title <- paste0(CI * 100, "% empirical confidence intervals, resampled")

# grouping
results %>% 
  group_by(method) %>%                          # create a grouped data frame
  summarize(lower=quantile(success, lower),     # lower limit of the CI
            middle=quantile(success, middle),   # middle (== median)
            upper=quantile(success, upper)) ->  # upper limit of the CI
    CIs                                         # assign to a new tibble, "CIs"

# reorder to plot
library(forcats) 
CIs %>% 
  mutate(method = fct_reorder(method, desc(middle))) %>% 
  ggplot(aes(x=method, ymin=lower, y=middle, ymax=upper, color=method)) + 
  scale_y_continuous(labels = scales::label_percent()) + 
  geom_pointrange(show.legend=FALSE) + 
  ylab("Classification success") + 
  coord_flip() + 
  theme_minimal() + 
  ggtitle(plot_title)



## ---- competitive-------------------------------------------------------------

# assign 1 point to the winner, 0 otherwise
# in the case of ties, split the point N ways
score_methods <- function(run) {

  top <- max(run$success)
  run$points <- as.numeric(run$success == top)
  score <- run$points / sum(run$points)
  names(score) <- run$method
  return(score)

}

# tally up the scores
runs %>% 
  purrr::map(score_methods) %>%  # map the function score_methods onto each run
  bind_rows() %>%                # bind the results as rows of a tibble 
  colSums -> scores              # compute the column sums and assign to scores

# tibble-ify for plotting 
outcome <- tibble(method=names(scores), 
                  score=scores,
                  scheme="unblocked") 

# plot it 
outcome %>% ggplot(aes(x=scheme, y=score, fill=method)) + 
            geom_col(position="fill") + 
            theme_minimal() + 
            ggtitle("Winning method, by resampling scheme")


