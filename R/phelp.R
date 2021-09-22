#' get the list of named documentation files for a package from its namespace
#' 
#' @param   a package name(space) as an unquoted string
#' 
#' @return  an interactive list of documented topics for that package
#'
#' @examples
#' help(package="Rscripts") 
#' phelp(Rscripts) 
#' phelp(MASS) 
#' 
#' @export
#' 
phelp <- function(...) {
  help(package=as.character(sapply(match.call()[-1], deparse)[1]))
}
