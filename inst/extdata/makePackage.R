# bootstrap functions for getting packages set up 
req <- function(p) require(p, character.only=TRUE) 
reqInstall <- function(p) { if (!req(p)) install.packages(p); req(p) }

# make sure skeletor is installed
reqInstall("skeletor") 

# create the package skeleton
skeletor(dir="../../",
         pkg="ExpDesign2021", 
         name="Triche Lab",
         email="trichelab@gmail.com",
         github="VanAndelInstitute")

# pulling the package from GitHub: 
if (FALSE) { 
  reqInstall("BiocManager") 
  BiocManager::install("VanAndelInstitute/ExpDesign2021") 
  help(package="ExpDesign2021")
}
