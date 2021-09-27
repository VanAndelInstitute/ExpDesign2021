#' fetch the Assignments spreadsheet for massaging and analysis 
#' 
#' @return a tibble, via googlesheets4::read_sheet() 
#' 
#' @import googlesheets4
#' 
#' @export
fetchAssignments <- function(url="https://docs.google.com/spreadsheets/d/1bJw_ad0PLQmLe4RvmRzsVnrOF-9Z9Xth9GTSdwO8B4M") {

  message("Calling gs4_deauth()...") 
  gs4_deauth() 
  
  message("Fetching data from ", url, " ...")
  res <- read_sheet(url)

  message("Done.")
  return(res)

}
