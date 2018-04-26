
#upgrade R
install.packages("installr") # install installr
library(installr)
install.packages("stringr")
library(stringr)
help(installr)
help(updateR)
updateR()

#get latest version of bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite()