.onAttach <- function (lib, pkgname="ExpDesign2021") {

  googlesheets4::gs4_deauth() 
  message("googlesheets4::gs4_deauth() called to reduce annoyances.")

  suppressWarnings(reqInstall("tidytext")) # test a recent dependency 
  message("ExpDesign2021 package and dependencies loaded.")

}
