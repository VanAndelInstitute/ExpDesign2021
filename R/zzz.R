.onAttach <- function (lib, pkgname="ExpDesign2021") {

  #suppressPackageStartupMessages()
  googlesheets4::gs4_deauth() 
  googlesheets4::gs4_user()
  suppressWarnings(reqInstall("tidytext")) # test a recent dependency 

}
