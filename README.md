# ExpDesign2021

[![Build Status](https://travis-ci.org/VanAndelInstitute/ExpDesign2021.png?branch=master)](https://travis-ci.org/VanAndelInstitute/ExpDesign2021)  [![codecov](https://codecov.io/gh/VanAndelInstitute/ExpDesign2021/branch/master/graph/badge.svg)](https://codecov.io/gh/VanAndelInstitute/ExpDesign2021)

## How to get started

If you're using a Mac, you may need to install [Xquartz](https://www.xquartz.org/) first; do so before proceeding. It will end like this:

[![xquartz](https://raw.githubusercontent.com/VanAndelInstitute/ExpDesign2021/main/inst/extdata/xquartz.png)](https://www.xquartz.org/)

Then pull up Rstudio or (if you're a caveman like me) R, and do the following: 

```
install.packages("remotes")
install.packages("BiocManager")
BiocManager::install("VanAndelInstitute/ExpDesign2021")
```

If you're asked which CRAN mirror to use, I suggest number (1), `cloud.r-project.org`.
If `BiocManager::install()` asks you which packages to update, I suggest (2), `CRAN packages only`.
If you have further issues with packages, you can try `install.packages("offending_package_name")` and see if that fixes the problem. 
If that doesn't work, and you're on a Mac, you may have to install some software to compile source code; [more details here](https://mac.r-project.org/tools/).
I strongly suggest that you avoid compiling your own packages unless you are writing packages yourself (even then, only if you have to).
If all goes well, you'll see something like this: 

[![installed](https://github.com/VanAndelInstitute/ExpDesign2021/raw/main/inst/extdata/InstallationResults.png)](https://github.com/VanAndelInstitute/ExpDesign2021)

The first couple of lines are checking to make sure that the `remotes` and `BiocManager` packages are up to date. 
They are not required, but if you have an older version of either package, it could be a problem, so we install the latest version to avoid this. 

## For developers

The repository includes a Makefile to facilitate some common tasks.

### Running tests

`$ make test`. Requires the [testthat](https://github.com/hadley/testthat) package. You can also specify a specific test file or files to run by adding a "file=" argument, like `$ make test file=logging`. `test_package` will do a regular-expression pattern match within the file names. See its documentation in the `testthat` package.

### Updating documentation

`$ make doc`. Requires the [roxygen2](https://github.com/klutometis/roxygen) package.
