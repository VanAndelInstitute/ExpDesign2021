#' require a package, and if it isn't installed, install it
#' 
#' @param p   the package to require/install
#' 
#' @return    whatever require() or install.packages() returns 
#'
#' @examples
#' reqInstall("BiocManager") 
#' 
#' @export
reqInstall <- function(p) {
  if (!.req(p)) install.packages(p)
  .req(p)
}


# helper
.req <- function(p) require(p, character.only=TRUE) 
