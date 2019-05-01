# Script de création du package 
path =  "/Users/ThomasHames/Documents/Sciences_Actuarielles/Mémoire/Code_Memoire_Hames_Thomas"


if(!require(devtools)){install.packages("devtools")}
if(!require(roxygen2)){devtools::install_github("klutometis/roxygen")}

devtools::create("ExtendedMort2DSmooth")

devtools::install_github("HamesThomas/Code_Memoire_Hames_Thomas")

findmypath=function(dir,file){
  path=system.file(dir,file,package="ExtendedMort2DSmooth")
  return(path)
}




devtools::install("ExtendedMort2DSmooth")
library("ExtendedMort2DSmooth")

Mort2Dsmooth_Binomial_logit()


devtools::install_github("HamesThomas/Code_Memoire_Hames_Thomas/ExtendedMort2DSmooth")
library("ExtendedMort2DSmooth")