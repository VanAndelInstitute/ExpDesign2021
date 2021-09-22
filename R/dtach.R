#' detach and unload a package.  usually works the way one might expect it to.
#' 
#' @param p   the package name, unquoted like in library()
#' 
#' @return    whatever base::detach() returns, if p was attached
#' 
#' @examples 
#' library(MASS)
#' dtach(MASS) 
#' 
#' @details
#' If a package is not directly loaded, base::detach will provide helpful 
#' feedback, in the eloquent yet restrained manner of a DEC VAX compiler.
#' 
#' @export
dtach <- function(p) { 

  pkg <- as.character(sapply(match.call()[-1], deparse)[1])
  attached <- base::loadedNamespaces() 
  if (pkg %in% attached) { 
    ns <- paste("package", pkg, sep=":")
    detach(ns, unload=TRUE, character.only=TRUE)
  } else { 
    message(pkg, " is not currently attached.")
  } 

}
