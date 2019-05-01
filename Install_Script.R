# Script de création du package 
path =  "/Users/ThomasHames/Documents/Sciences_Actuarielles/Mémoire/Code_Memoire_Hames_Thomas"


if(!require(devtools)){install.packages("devtools")}
if(!require(roxygen2)){devtools::install_github("klutometis/roxygen")}


devtools::install_github("HamesThomas/Code_Memoire_Hames_Thomas/ExtendedMort2DSmooth")
library("ExtendedMort2DSmooth")

# Mort2DSmooth_Binomial_logit, Mort2DSmooth_Binomial_cloglog,  Mort2DSmooth_poisson_log ready to use ! 