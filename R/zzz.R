#' A package for the 2021 edition of Experimental Design and Biostatistics.
#'
#' @examples 
#' 
#' # phelp == 'package help'
#' phelp(ExpDesign2021) 
#' 
#' # automate fetching assignment data 
#' fetchAssignments()
#'
#' # list object sizes (lsos) and pipe the output into head (default: 6 lines)
#' lsos() %>% head 
#' 
#' @rdname ExpDesign2021 
.onAttach <- function (lib, pkgname="ExpDesign2021") {

  tidymodels_prefer(quiet=FALSE)
  conflicted::conflict_prefer("filter", "dplyr")
  conflicted::conflict_prefer("mutate", "dplyr")
  googlesheets4::gs4_deauth() 
  googlesheets4::gs4_user()

}
