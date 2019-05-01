rm(list = ls())
setwd("~/Documents/Sciences_Actuarielles/Mémoire/Memoire/DATA_R")

LifeTableMen <- read.table("LifeTableMen_Belgium.txt",
                             header = TRUE ,sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
#install.packages("gnm")
#install.packages("StMoMo")
#install.packages("stats")
library(gnm)
library(forecast)
library(StMoMo)
library(fanplot)
library(demography)
library(arima)
library(MortalitySmooth)
library(jpeg)
library(fields)
library(spatstat)
##StMoMo vignette
#kronecker product existe


#########################################
###        Traitement des data        ###
#########################################



LifeTableMen$qx = as.numeric(as.character(LifeTableMen$qx))
LifeTableMen$lx = as.numeric(as.character(LifeTableMen$lx))
LifeTableMen$dx = as.numeric(as.character(LifeTableMen$dx))
LifeTableMen$ex = as.numeric(as.character(LifeTableMen$ex))
LifeTableMen$mx = as.numeric(as.character(LifeTableMen$mx))
LifeTableMen$ax = as.numeric(as.character(LifeTableMen$ax))
LifeTableMen$Lx = as.numeric(as.character(LifeTableMen$Lx))
LifeTableMen$Year = as.numeric(as.character(LifeTableMen$Year))

Age = c(50:90)
Year = c(1960:2015)

n_y = 2015-1960 +1
n_a = length(Age)
n_c = n_y + n_a -1
n = n_a * n_y


## Création Matrice Décès
LT_M_1960 = LifeTableMen[LifeTableMen$Year >= 1960,]
LT_M_1960[LT_M_1960$Age == "110+",]$qx = 1
LT_M_1960$Age = as.numeric(as.character(LT_M_1960$Age))

LT_M_1960_50 = LT_M_1960[LT_M_1960$Age >= 50,]
LT_M_1960_50_90 = LT_M_1960_50[LT_M_1960_50$Age <= 90,]
LT_M_1960_50_90 = na.exclude(LT_M_1960_50_90)

n_a = length(LT_M_1960_50_90$Year[LT_M_1960_50_90$Year == 1960])
n_y = max(LT_M_1960_50_90$Year) - min(LT_M_1960_50_90$Year)+1
n = n_a * n_y

plot(log(LifeTableMen[LifeTableMen$Year == 2015,]$mx), type = "l", xlab = "age",ylab = "Logarithm of the observed force of mortality", main = "Belgian Population in 1990")
plot((LifeTableMen[LifeTableMen$Year == 1990,]$mx), type = "l", xlab = "age",ylab = "Observed force of mortality", main = "Belgian Population in 1990")





#####################################################################################
#####    Ecriture des modèles M3/M5/M6/M7/M8 avec leurs applications via glm    #####
#####################################################################################


InitialExpo = LT_M_1960_50_90$lx
Qx = LT_M_1960_50_90$dx/LT_M_1960_50_90$lx
CentralExpo = -LT_M_1960_50_90$lx*Qx / (log(1-Qx))
dx = LT_M_1960_50_90$dx
mu = -log(1-Qx)


## Création de la matrice de décès
dx_Mat = matrix(dx, nrow = n_a, ncol = n_y, byrow = FALSE )
Qx_Mat = matrix(Qx, nrow = n_a, ncol = n_y, byrow = FALSE )
mu_Mat = matrix(mu, nrow = n_a, ncol = n_y, byrow = FALSE )
## Création de la matrice d'exposition
InitialExposure_Mat = matrix(InitialExpo, nrow = n_a, ncol = n_y, byrow = FALSE)
CentralExposure_Mat = matrix(CentralExpo, nrow = n_a, ncol = n_y, byrow = FALSE)





      ##################
      ### Modèle APC ###    Done
      ##################


Model_Matrix_APC = function(n_a, n_y)
{  
  n = n_a * n_y
        #Création X_a
  Identity_n_a = matrix(0, ncol = n_a, nrow = n_a, byrow = TRUE)
  diag(Identity_n_a) = 1
  X_na = Identity_n_a 
  for (i in 2:n_y)
  {
    X_na = rbind(X_na,Identity_n_a)
  }
  
        #Création X_y
  Identity_n_y = matrix(0, ncol = n_y, nrow = n_y, byrow = TRUE)
  diag(Identity_n_y) = 1
  serie1_n_a = rep(1, n_a)
  X_ny = matrix(0 , ncol = n_y, nrow = n, byrow = TRUE )
  for (j in (1:n_y))
  {
    X_ny[(1:n_a) + n_a*(j-1) ,j] = 1
  }
 
  
  
  ### Rework de X_c! 
  n_c = n_a + n_y -1
  X_c = matrix(0 , ncol = n_c, nrow = n, byrow = TRUE )
  
  count = 1
  for(i in (0:(n_y-1)))
  {
    for(j in (0:(n_a-1)))
    {
      vect = rep(0,n_c)
      vect[n_a - j +i] = 1
      
      X_c[count,] = vect
      count = count +1 
    }
    
  }
  X = cbind(X_na,X_ny,X_c)
  
  return (X)
}


X_APC = Model_Matrix_APC(n_a, n_y)

GLM1_APC = glm(dx ~ -1 + X_APC + offset(log(CentralExpo)), family = poisson(link = "log"), epsilon = 1.4e-8)

GLM2_APC = glm(Qx ~ -1 + X_APC, family = binomial(link = "cloglog"), weight = InitialExpo, epsilon = 1.4e-8)

GLM3_APC = glm(Qx ~ -1 + X_APC, family = binomial(link = "logit"), weight = InitialExpo, epsilon = 1.3e-8)

APC = c(GLM1_APC$deviance, GLM2_APC$deviance, GLM3_APC$deviance)
APC

APC_aic = c(GLM1_APC$aic, GLM2_APC$aic, GLM3_APC$aic)
APC_aic


NumberConstraints_APC = dim(X_APC)[2] - qr(X_APC)$rank
NumberConstraints_APC

#Remplace les NA que R a mis pour les contraintes par des 0
GLM1_APC$coefficients[97] = 0
GLM1_APC$coefficients[192] = 0
GLM1_APC$coefficients[193] = 0

GLM2_APC$coefficients[97] = 0
GLM2_APC$coefficients[192] = 0
GLM2_APC$coefficients[193] = 0

GLM3_APC$coefficients[97] = 0
GLM3_APC$coefficients[192] = 0
GLM3_APC$coefficients[193] = 0

## Transformation des valeurs des coefficients pour que ça corresponde avec les contraintes posées en général
## alors que R pose lui même ses propres contraintes

#Création matrice H
H = c(rep(0,n_a), rep(1,n_y), rep(0,n_c),rep(0, n_a + n_y), rep(1,n_c), rep(0, n_a + n_y), c(1:n_c))
H = matrix(H, nrow = n_a + n_y + n_c, ncol = 3)
H = t(H)

H_R = c(rep(0,n_a+n_y), 1, rep(0,n_c -1), rep(0, n_a + n_y + n_c -2), 1,0,rep(0, n_a + n_y + n_c -1), 1)
H_R = matrix(H_R, nrow = n_a + n_y + n_c, ncol = 3)
H_R = t(H_R)

GLM1_APC$coefficients = solve(t(X_APC) %*% X_APC + t(H) %*% H) %*% (t(X_APC) %*% X_APC + t(H) %*% H_R) %*% GLM1_APC$coefficients
GLM2_APC$coefficients = solve(t(X_APC) %*% X_APC + t(H) %*% H) %*% (t(X_APC) %*% X_APC + t(H) %*% H_R) %*% GLM2_APC$coefficients
GLM3_APC$coefficients = solve(t(X_APC) %*% X_APC + t(H) %*% H) %*% (t(X_APC) %*% X_APC + t(H) %*% H_R) %*% GLM3_APC$coefficients



Compute_resids = function(fit_values, FamillyType){
  
  if(FamillyType == "Poisson"){
    mu_fit = -log(1-fit_values/LT_M_1960_50_90$lx)
    mu_fit_matrix = matrix(c(mu_fit), nrow = n_a, ncol = n_y, byrow = FALSE)
  }
  else if(FamillyType == "Binomial"){
    mu_fit = -log(1-fit_values)
    mu_fit_matrix = matrix(c(mu_fit), nrow = n_a, ncol = n_y, byrow = FALSE)
  }
  
  resids =  (mu_Mat - mu_fit_matrix)
  return(resids)
}
Resids_Plot = function(residuals_mu, blur, Main_name)
{
  picture <- as.matrix(blur(as.im(residuals_mu), sigma=blur))
  image.plot(picture, col=rainbow(500, s= 1, v = 1), midpoint = TRUE, main=paste("mortality rate residuals for model: ", Main_name), asp = 1, xlab = "Age", ylab = "Year", axes = F)
  axis(1, at=seq(0,1,0.25), labels=c(50,60,70,80,90))
  axis(2, at=seq(0,1, 0.2), labels = c(1960,1970,1980,1990,2000,2010))
  
  
  drape.plot(c(50:90), 1:ncol(picture),col=rainbow(100, s= 1, v = 1) ,picture, border=NA, theta= -35, phi=35, main=paste("standardized residuals for model: ", Main_name), xlab = "Age", ylab = "Year", zlab = "residuals value")
}


## then with Poisson assumption, we find
    ##for the log link
eta_APC_1 = (X_APC %*% GLM1_APC$coefficients)
mu_log_APC_1 = exp(eta_APC_1)
mu_log_APC_1 = matrix(mu_log_APC_1, nrow = n_a, ncol = n_y, byrow = FALSE)  #Mat format
persp(Age,Year, mu_log_APC_1, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "APC_log")
Resids_Plot(residuals_mu = matrix(c(rstandard(GLM1_APC)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "APC_poisson_log" )

## then with Binomial assumption, we find :
    ##For the cloglog link
eta_APC_2 = (X_APC %*% GLM2_APC$coefficients)
qx_cloglog_APC_2 = 1-exp(-exp(eta_APC_2))
qx_cloglog_APC_2 = matrix(qx_cloglog_APC_2, nrow = n_a, ncol = n_y, byrow = FALSE)  #Mat format
mu_cloglog_APC_2 = -log(1-qx_cloglog_APC_2)
persp(Age,Year, mu_cloglog_APC_2, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "APC_cloglog")
Resids_Plot(residuals_mu = matrix(c(rstandard(GLM2_APC)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "APC_Binomial_cloglog" )

    ##For the logit link
eta_APC_3 = (X_APC %*% GLM3_APC$coefficients)
qx_logit_APC_3 = (exp(eta_APC_3))/(1 + exp(eta_APC_3))
qx_logit_APC_3 = matrix((qx_logit_APC_3), nrow = n_a, ncol = n_y, byrow = FALSE)  #Mat format
mu_logit_APC_3 = -log(1-qx_logit_APC_3)
persp(Age,Year, mu_logit_APC_3, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "APC_logit")
Resids_Plot(residuals_mu = matrix(c(rstandard(GLM3_APC)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "APC_Binomial_logit" )




Deviance.FromBinomial = function(dx_obs, qx_fit, initial_ETR, devType){
  dx_fit = qx_fit * InitialExpo
  
  if(devType == "Poisson"){
    dev = 2 * (sum((dx_obs * log(dx_obs/dx_fit) - (dx_obs - dx_fit))))
  }
  else if(devType == "Binomial"){
    dev = 2 * (sum((dx_obs)* log(dx_obs/dx_fit) + (InitialExpo - dx_obs) * log((InitialExpo - dx_obs) / (InitialExpo-dx_fit))))
  }

  return (dev)
}

APC_dev_Poisson = c(GLM1_APC$deviance,
                    Deviance.FromBinomial(dx_obs = dx, qx_fit = GLM2_APC$fitted.values, initial_ETR = InitialExpo, devType = "Poisson"),
                    Deviance.FromBinomial(dx_obs = dx, qx_fit = GLM3_APC$fitted.values, initial_ETR = InitialExpo, devType = "Poisson"))




#Fonction permettant de plotter les coefs estimé pour chacun des 3 modèles
Plot_Superpose_Coef = function(TheRange, GLM1, GLM2, GLM3, GLM1_name, GLM2_name, GLM3_name ,main_name, y_name, pos_legend)
{

  min_value = min(GLM1$coefficients[TheRange], GLM2$coefficients[TheRange], GLM3$coefficients[TheRange])
  max_value = max(GLM1$coefficients[TheRange], GLM2$coefficients[TheRange], GLM3$coefficients[TheRange])
  for(i in (1:3)){
    if( i == 1 ){
      plot(GLM1$coefficients[TheRange], main = main_name, ylab = y_name, ylim = c(min_value, max_value))
    }
    else if (i== 2)
    {
      lines(GLM2$coefficients[TheRange], col = "red")
    }
    else if ( i == 3)
    {
      lines(GLM3$coefficients[TheRange], col = "green")
    }

  }
  legend(pos_legend, paste(c(GLM1_name, GLM2_name, GLM3_name)), lty=c(NA,1,1),pch=c(1, NA, NA) , col = c("black", "red", "green"))
  
}

Plot_Superpose_Coef(TheRange = c(1:n_a), GLM1_APC, GLM2_APC, GLM3_APC, GLM1_name = "APC_Poisson_log" , GLM2_name = "APC_Binomial_cloglog" , GLM3_name = "APC_Binomial_logit" , main_name = "APC", y_name = "Age_coef_Beta", pos_legend = "topleft")

Plot_Superpose_Coef(TheRange = c((n_a +1):(n_a + n_y)), GLM1_APC, GLM2_APC, GLM3_APC, GLM1_name = "APC_Poisson_log" , GLM2_name = "APC_Binomial_cloglog" , GLM3_name = "APC_Binomial_logit" , main_name = "APC", y_name = "Year_coef_Kappa", pos_legend = 'topright')

Plot_Superpose_Coef(TheRange = c((n_a + n_y + 1): (n_a + n_y + n_c)), GLM1_APC, GLM2_APC, GLM3_APC, GLM1_name = "APC_Poisson_log" , GLM2_name = "APC_Binomial_cloglog" , GLM3_name = "APC_Binomial_logit" , main_name = "APC", y_name = "Cohort_coef_Gamma", pos_legend = "topleft")
    




      ##################
      ### Modèle CBD ###    Done
      ##################

Model_Matrix_CBD = function(n_a, n_y)
{
  n = n_a * n_y
          #Création X1
  Identity_n_y = matrix(0, ncol = n_y, nrow = n_y, byrow = TRUE)
  diag(Identity_n_y) = 1
  serie1_n_a = rep(1, n_a)
  X_1 = matrix(0 , ncol = n_y, nrow = n, byrow = TRUE )
  
  for (j in 1:n_y)
  {
    X_1[(1:n_a) + n_a*(j-1) ,j] = 1
  }
  
        #Création X2
  Centred_Age_Vector = (50:90) - mean(50:90)
  X_2 = matrix(0 , ncol = n_y, nrow = n, byrow = TRUE )
  
  for (j in 1:n_y)
  {
    for( i in 1:length(Centred_Age_Vector))
    {
      X_2[i + (n_a) * (j-1),j] = Centred_Age_Vector[i]
    }
  }
  
  X = cbind(X_1,X_2)
  return (X)
}

X_CBD = Model_Matrix_CBD(n_a,n_y)

GLM1_CBD = glm(dx ~ -1 +  X_CBD + offset(log(CentralExpo)), family = poisson(link = "log"))
GLM2_CBD = glm(Qx ~ -1 + X_CBD, family = binomial(link = "cloglog"), weight = InitialExpo)
GLM3_CBD = glm(Qx ~ -1 + X_CBD, family = binomial(link = "logit"), weight = InitialExpo)

CBD = c(GLM1_CBD$deviance, GLM2_CBD$deviance, GLM3_CBD$deviance)
CBD

#Number of constraints needed to fit the model
NumberConstraints_CBD = dim(X_CBD)[2] - qr(X_CBD)$rank
NumberConstraints_CBD


## then with Poisson assumption, we find :
      ##For the log link
eta_CBD_1 = (X_CBD %*% GLM1_CBD$coefficients)
mu_log_CBD_1 = exp(eta_CBD_1)
mu_log_CBD_1 = matrix(mu_log_CBD_1, nrow = n_a, ncol = n_y, byrow = FALSE)  #Mat format
persp(Age,Year, mu_log_CBD_1, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "CBD_log")


## then with Binomial assumption, we find :
    ##For the cloglog link
eta_CBD_2 = (X_CBD %*% GLM2_CBD$coefficients)
qx_cloglog_CBD_2 = 1-exp(-exp(eta_CBD_2))
qx_cloglog_CBD_2 = matrix(qx_cloglog_CBD_2, nrow = n_a, ncol = n_y, byrow = FALSE)  #Mat format
mu_cloglog_CBD_2 = -log(1-qx_cloglog_CBD_2)
persp(Age,Year, mu_cloglog_CBD_2, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "CBD_cloglog")

    ##For the logit link
eta_CBD_3 = (X_CBD %*% GLM3_CBD$coefficients)
qx_logit_CBD_3 = (exp(eta_CBD_3))/(1 + exp(eta_CBD_3))
qx_logit_CBD_3 = matrix(qx_logit_CBD_3, nrow = n_a, ncol = n_y, byrow = FALSE)  #Mat format
mu_logit_CBD_3 = -log(1-qx_logit_CBD_3)
persp(Age,Year, mu_logit_CBD_3, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "CBD_logit")

CBD_dev_Poisson = c(GLM1_CBD$deviance,
                    Deviance.FromBinomial(dx_obs = dx, qx_fit = GLM2_CBD$fitted.values, initial_ETR = InitialExpo, devType = "Poisson"),
                    Deviance.FromBinomial(dx_obs = dx, qx_fit = GLM3_CBD$fitted.values, initial_ETR = InitialExpo, devType = "Poisson"))

CBD_aic = c(GLM1_CBD$aic, GLM2_CBD$aic, GLM3_CBD$aic)
CBD_aic




Resids_Plot(residuals_mu = matrix(c(rstandard(GLM1_CBD)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "CBD_poisson_log" )
Resids_Plot(residuals_mu = matrix(c(rstandard(GLM2_CBD)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "CBD_Binomial_cloglog" )
Resids_Plot(residuals_mu = matrix(c(rstandard(GLM3_CBD)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "CBD_Binomial_logit" )


Plot_Superpose_Coef(TheRange = c(1:n_y), GLM1_CBD, GLM2_CBD, GLM3_CBD, GLM1_name = "CBD_Poisson_log" , GLM2_name = "CBD_Binomial_cloglog" , GLM3_name = "CBD_Binomial_logit" , main_name = "CBD", y_name = "Year_coef_Kappa_1", pos_legend = "topright")

Plot_Superpose_Coef(TheRange = c((n_y +1):(2* n_y)), GLM1_CBD, GLM2_CBD, GLM3_CBD, GLM1_name = "CBD_Poisson_log" , GLM2_name = "CBD_Binomial_cloglog" , GLM3_name = "CBD_Binomial_logit" , main_name = "CBD", y_name = "Year_coef_Kappa_2", pos_legend = 'topleft')



      ###################################
      ### Modèle CBD  + Cohort effect ###   Done
      ###################################


Model_Matrix_CBD_cohort = function(n_a, n_y)
{
  n = n_a * n_y
  
        #Création X1
  Identity_n_y = matrix(0, ncol = n_y, nrow = n_y, byrow = TRUE)
  diag(Identity_n_y) = 1
  serie1_n_a = rep(1, n_a)
  X_1 = matrix(0 , ncol = n_y, nrow = n, byrow = TRUE )
  
  for (j in 1:n_y)
  {
    X_1[(1:n_a) + n_a*(j-1) ,j] = 1
  }
  
  #Création X2
  Centred_Age_Vector = (50:90) - mean(50:90)
  X_2 = matrix(0 , ncol = n_y, nrow = n, byrow = TRUE )
  
  for (j in 1:n_y)
  {
    for( i in 1:length(Centred_Age_Vector))
    {
      X_2[i + (n_a) * (j-1),j] = Centred_Age_Vector[i]
    }
  }
  

  ### Rework de X_c! 
  n_c = n_a + n_y -1
  X_c = matrix(0 , ncol = n_c, nrow = n_a * n_y, byrow = TRUE )
  
  count = 1
  for(i in (0:(n_y-1)))
  {
    for(j in (0:(n_a-1)))
    {
      vect = rep(0,n_c)
      vect[n_a - j +i] = 1
      
      X_c[count,] = vect
      count = count +1 
    }
    
  }
  
  
  X = cbind(X_1,X_2,X_c)
  return(X)
}

X_CBD_cohort = Model_Matrix_CBD_cohort(n_a,n_y)

GLM1_CBD_cohort = glm(dx ~ -1 +  X_CBD_cohort + offset(log(CentralExpo)), family = poisson(link = "log"))
GLM2_CBD_cohort = glm(Qx ~ -1 + X_CBD_cohort, family = binomial(link = "cloglog"), weight = InitialExpo)
GLM3_CBD_cohort = glm(Qx ~ -1 + X_CBD_cohort, family = binomial(link = "logit"), weight = InitialExpo)

CBD_cohort = c(GLM1_CBD_cohort$deviance, GLM2_CBD_cohort$deviance, GLM3_CBD_cohort$deviance)

NumberConstraints_CBD_cohort = dim(X_CBD_cohort)[2] - qr(X_CBD_cohort)$rank
NumberConstraints_CBD_cohort


#Remplace les NA que R a mis pour les contraintes par des 0
GLM1_CBD_cohort$coefficients[207] = 0
GLM1_CBD_cohort$coefficients[208] = 0

GLM2_CBD_cohort$coefficients[207] = 0
GLM2_CBD_cohort$coefficients[208] = 0

GLM3_CBD_cohort$coefficients[207] = 0
GLM3_CBD_cohort$coefficients[208] = 0


#Création matrice H
H = c(rep(0,2*n_y), rep(1,n_c), rep(0, 2* n_y), c(1:n_c))
H = matrix(H, nrow = 2*n_y + n_c, ncol = 2)
H = t(H)

H_R = c(rep(0,2*n_y + n_c -2), 1, 0, rep(0, 2*n_y + n_c -1), 1) 
H_R = matrix(H_R, nrow = 2*n_y + n_c, ncol = 2)
H_R = t(H_R)

GLM1_CBD_cohort$coefficients = solve(t(X_CBD_cohort) %*% X_CBD_cohort + t(H) %*% H) %*% (t(X_CBD_cohort) %*% X_CBD_cohort + t(H) %*% H_R) %*% GLM1_CBD_cohort$coefficients
GLM2_CBD_cohort$coefficients = solve(t(X_CBD_cohort) %*% X_CBD_cohort + t(H) %*% H) %*% (t(X_CBD_cohort) %*% X_CBD_cohort + t(H) %*% H_R) %*% GLM2_CBD_cohort$coefficients
GLM3_CBD_cohort$coefficients = solve(t(X_CBD_cohort) %*% X_CBD_cohort + t(H) %*% H) %*% (t(X_CBD_cohort) %*% X_CBD_cohort + t(H) %*% H_R) %*% GLM3_CBD_cohort$coefficients

## then with Poisson assumption, we find
    ##For the log link
eta_CBD_cohort_1 = (X_CBD_cohort %*% GLM1_CBD_cohort$coefficients)
mu_log_CBD_cohort_1 = exp(eta_CBD_cohort_1)
mu_log_CBD_cohort_1 = matrix(mu_log_CBD_cohort_1, nrow = n_a, ncol = n_y, byrow = FALSE)  #Mat format
persp(Age,Year, mu_log_CBD_cohort_1, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "CBD_cohort_log")

## then with Binomial assumption, we find :
    ##For the cloglog link
eta_CBD_cohort_2 = (X_CBD_cohort %*% GLM2_CBD_cohort$coefficients)
qx_cloglog_CBD_cohort_2 = 1-exp(-exp(eta_CBD_cohort_2))
qx_cloglog_CBD_cohort_2 = matrix(qx_cloglog_CBD_cohort_2, nrow = n_a, ncol = n_y, byrow = FALSE)  #Mat format
mu_cloglog_CBD_cohort_2 = -log(1-qx_cloglog_CBD_cohort_2)
persp(Age,Year, mu_cloglog_CBD_cohort_2, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "CBD_cohort_cloglog")


    ##For the logit link
eta_CBD_cohort_3 = (X_CBD_cohort %*% GLM3_CBD_cohort$coefficients)
qx_logit_CBD_cohort_3 = (exp(eta_CBD_cohort_3))/(1 + exp(eta_CBD_cohort_3))
qx_logit_CBD_cohort_3 = matrix(qx_logit_CBD_cohort_3, nrow = n_a, ncol = n_y, byrow = FALSE)  #Mat format
mu_logit_CBD_cohort_3 = -log(1-qx_logit_CBD_cohort_3)
persp(Age,Year, mu_logit_CBD_cohort_3, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "CBD_cohort_logit")


CBD_cohort_dev_Poisson = c(GLM1_CBD_cohort$deviance,
                    Deviance.FromBinomial(dx_obs = dx, qx_fit = GLM2_CBD_cohort$fitted.values, initial_ETR = InitialExpo, devType = "Poisson"),
                    Deviance.FromBinomial(dx_obs = dx, qx_fit = GLM3_CBD_cohort$fitted.values, initial_ETR = InitialExpo, devType = "Poisson"))

CBD_cohort_aic = c(GLM1_CBD_cohort$aic, GLM2_CBD_cohort$aic, GLM3_CBD_cohort$aic)
CBD_cohort_aic

Resids_Plot(residuals_mu = matrix(c(rstandard(GLM1_CBD_cohort)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "CBD_cohort_poisson_log" )
Resids_Plot(residuals_mu = matrix(c(rstandard(GLM2_CBD_cohort)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "CBD_cohort_Binomial_cloglog" )
Resids_Plot(residuals_mu = matrix(c(rstandard(GLM3_CBD_cohort)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "CBD_cohort_Binomial_logit" )



Plot_Superpose_Coef(TheRange = c(1:n_y), GLM1_CBD_cohort, GLM2_CBD_cohort, GLM3_CBD_cohort, GLM1_name = "CBD_cohort_Poisson_log" , GLM2_name = "CBD_cohort_Binomial_cloglog" , GLM3_name = "CBD_cohort_Binomial_logit" , main_name = "CBD_cohort", y_name = "Year_coef_Kappa_1", pos_legend = "topright")

Plot_Superpose_Coef(TheRange = c((n_y +1):(2* n_y)), GLM1_CBD_cohort, GLM2_CBD_cohort, GLM3_CBD_cohort, GLM1_name = "CBD_cohort_Poisson_log" , GLM2_name = "CBD_cohort_Binomial_cloglog" , GLM3_name = "CBD_cohort_Binomial_logit" , main_name = "CBD_cohort", y_name = "Year_coef_Kappa_2", pos_legend = 'topleft')

Plot_Superpose_Coef(TheRange = c((n_y + n_y + 1): (n_y + n_y + n_c)), GLM1_CBD_cohort, GLM2_CBD_cohort, GLM3_CBD_cohort, GLM1_name = "CBD_cohort_Poisson_log" , GLM2_name = "CBD_cohort_Binomial_cloglog" , GLM3_name = "CBD_cohort_Binomial_logit" , main_name = "CBD_cohort", y_name = "Cohort_coef_Gamma", pos_legend = "topleft")



    ######################################################
    ### Modèle CBD  + Cohort effect + quadratic effect ###    Done
    ######################################################


Model_Matrix_CBD_cohort_quad = function(n_a, n_y) 
{  
  #Création X1
  Identity_n_y = matrix(0, ncol = n_y, nrow = n_y, byrow = TRUE)
  diag(Identity_n_y) = 1
  serie1_n_a = rep(1, n_a)
  X_1 = matrix(0 , ncol = n_y, nrow = n_a * n_y, byrow = TRUE )
  
  for (j in 1:n_y)
  {
    X_1[(1:n_a) + n_a*(j-1) ,j] = 1
  }
  
  #Création X2
  Centred_Age_Vector = (50:90) - mean(50:90)
  X_2 = matrix(0 , ncol = n_y, nrow = n_a * n_y, byrow = TRUE )
  
  for (j in 1:n_y)
  {
    for( i in 1:length(Centred_Age_Vector))
    {
      X_2[i + (n_a) * (j-1),j] = Centred_Age_Vector[i]
    }
  }
  
  # quadratic effect
  
  Quadratic_Vector = ((50:90) - mean(50:90))^2 - (sd(50:90))^2
  X_q = matrix(0 , ncol = n_y, nrow = n_a * n_y, byrow = TRUE )
  
  for (j in 1:n_y)
  {
    for( i in 1:length(Quadratic_Vector))
    {
      X_q[i + (n_a) * (j-1),j] = Quadratic_Vector[i]
    }
  }
  
  
  # Cohort matrix
  ### Rework de X_c! 
  n_c = n_a + n_y -1
  X_c = matrix(0 , ncol = n_c, nrow = n_a * n_y, byrow = TRUE )
  
  count = 1
  for(i in (0:(n_y-1)))
  {
    for(j in (0:(n_a-1)))
    {
      vect = rep(0,n_c)
      vect[n_a - j +i] = 1
      
      X_c[count,] = vect
      count = count +1 
    }
    
  }
  
  X = cbind(X_1, X_2, X_q ,X_c)
}  

X_CBD_cohort_quad = Model_Matrix_CBD_cohort_quad(n_a, n_y)

GLM1_CBD_cohort_quad = glm(dx ~ -1 +  X_CBD_cohort_quad + offset(log(CentralExpo)), family = poisson(link = "log"), epsilon = 9e-6)
GLM2_CBD_cohort_quad = glm(Qx ~ -1 + X_CBD_cohort_quad, family = binomial(link = "cloglog"), weight = InitialExpo, epsilon = 9e-6)
GLM3_CBD_cohort_quad = glm(Qx ~ -1 + X_CBD_cohort_quad, family = binomial(link = "logit"), weight = InitialExpo, epsilon = 9e-6)

CBD_cohort_quad = c(GLM1_CBD_cohort_quad$deviance, GLM2_CBD_cohort_quad$deviance, GLM3_CBD_cohort_quad$deviance)

NumberConstraints_CBD_cohort_quad = dim(X_CBD_cohort_quad)[2] - qr(X_CBD_cohort_quad)$rank
NumberConstraints_CBD_cohort_quad

#Remplace les NA que R a mis pour les contraintes par des 0
GLM1_CBD_cohort_quad$coefficients[262] = 0
GLM1_CBD_cohort_quad$coefficients[263] = 0
GLM1_CBD_cohort_quad$coefficients[264] = 0

GLM2_CBD_cohort_quad$coefficients[262] = 0
GLM2_CBD_cohort_quad$coefficients[263] = 0
GLM2_CBD_cohort_quad$coefficients[264] = 0

GLM3_CBD_cohort_quad$coefficients[262] = 0
GLM3_CBD_cohort_quad$coefficients[263] = 0
GLM3_CBD_cohort_quad$coefficients[264] = 0


#Création matrice H
square = (1:n_c)^2
H = c(rep(0,3*n_y), rep(1,n_c), rep(0, 3* n_y), c(1:n_c), rep(0, 3*n_y), square )
H = matrix(H, nrow = 3*n_y + n_c, ncol = 3)
H = t(H)

H_R = c(rep(0,3*n_y + n_c -3),1, 0, 0, rep(0,3*n_y + n_c -2),1, 0, rep(0,3*n_y + n_c -1),1)
H_R = matrix(H_R, nrow = 3*n_y + n_c, ncol = 3)
H_R = t(H_R)

GLM1_CBD_cohort_quad$coefficients = solve(t(X_CBD_cohort_quad) %*% X_CBD_cohort_quad + t(H) %*% H) %*% (t(X_CBD_cohort_quad) %*% X_CBD_cohort_quad + t(H) %*% H_R) %*% GLM1_CBD_cohort_quad$coefficients
GLM2_CBD_cohort_quad$coefficients = solve(t(X_CBD_cohort_quad) %*% X_CBD_cohort_quad + t(H) %*% H) %*% (t(X_CBD_cohort_quad) %*% X_CBD_cohort_quad + t(H) %*% H_R) %*% GLM2_CBD_cohort_quad$coefficients
GLM3_CBD_cohort_quad$coefficients = solve(t(X_CBD_cohort_quad) %*% X_CBD_cohort_quad + t(H) %*% H) %*% (t(X_CBD_cohort_quad) %*% X_CBD_cohort_quad + t(H) %*% H_R) %*% GLM3_CBD_cohort_quad$coefficients



## then with Poisson assumption, we find :
    ##For the log link
eta_CBD_cohort_quad_1 = (X_CBD_cohort_quad %*% GLM1_CBD_cohort_quad$coefficients)
mu_log_CBD_cohort_quad_1 = exp(eta_CBD_cohort_quad_1)
mu_log_CBD_cohort_quad_1 = matrix(mu_log_CBD_cohort_quad_1, nrow = n_a, ncol = n_y, byrow = FALSE)  #Mat format
persp(Age,Year, mu_log_CBD_cohort_quad_1, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "CBD_cohort_quad_log")


## then with Binomial assumption, we find :
    ##For the cloglog link
eta_CBD_cohort_quad_2 = (X_CBD_cohort_quad %*% GLM2_CBD_cohort_quad$coefficients)
qx_cloglog_CBD_cohort_quad_2 = 1-exp(-exp(eta_CBD_cohort_quad_2))
qx_cloglog_CBD_cohort_quad_2 = matrix(qx_cloglog_CBD_cohort_quad_2, nrow = n_a, ncol = n_y, byrow = FALSE)  #Mat format
mu_cloglog_CBD_cohort_quad_2 = -log(1-qx_cloglog_CBD_cohort_quad_2)
persp(Age,Year, mu_cloglog_CBD_cohort_quad_2, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "CBD_cohort_quad_cloglog")

    ##For the logit link
eta_CBD_cohort_quad_3 = (X_CBD_cohort_quad %*% GLM3_CBD_cohort_quad$coefficients)
qx_logit_CBD_cohort_quad_3 = (exp(eta_CBD_cohort_quad_3))/(1 + exp(eta_CBD_cohort_quad_3))
qx_logit_CBD_cohort_quad_3 = matrix(qx_logit_CBD_cohort_quad_3, nrow = n_a, ncol = n_y, byrow = FALSE)  #Mat format
mu_logit_CBD_cohort_quad_3 = -log(1-qx_logit_CBD_cohort_quad_3)
persp(Age,Year, mu_logit_CBD_cohort_quad_3, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "CBD_cohort_quad_logit")


CBD_cohort_quad_dev_Poisson = c(GLM1_CBD_cohort_quad$deviance,
                               Deviance.FromBinomial(dx_obs = dx, qx_fit = GLM2_CBD_cohort_quad$fitted.values, initial_ETR = InitialExpo, devType = "Poisson"),
                               Deviance.FromBinomial(dx_obs = dx, qx_fit = GLM3_CBD_cohort_quad$fitted.values, initial_ETR = InitialExpo, devType = "Poisson"))

CBD_cohort_quad_aic = c(GLM1_CBD_cohort_quad$aic, GLM2_CBD_cohort_quad$aic, GLM3_CBD_cohort_quad$aic)
CBD_cohort_quad_aic

Resids_Plot(residuals_mu = matrix(c(rstandard(GLM1_CBD_cohort_quad)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "CBD_cohort_quad_poisson_log" )
Resids_Plot(residuals_mu = matrix(c(rstandard(GLM2_CBD_cohort_quad)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "CBD_cohort_quad_Binomial_cloglog" )
Resids_Plot(residuals_mu = matrix(c(rstandard(GLM3_CBD_cohort_quad)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "CBD_cohort_quad_Binomial_logit" )


Plot_Superpose_Coef(TheRange = c(1:n_y), GLM1_CBD_cohort_quad, GLM2_CBD_cohort_quad, GLM3_CBD_cohort_quad, GLM1_name = "CBD_cohort_quad_Poisson_log" , GLM2_name = "CBD_cohort_quad_Binomial_cloglog" , GLM3_name = "CBD_cohort_quad_Binomial_logit" , main_name = "CBD_cohort_quad", y_name = "Year_coef_Kappa_1", pos_legend = "topright")

Plot_Superpose_Coef(TheRange = c((n_y +1):(2* n_y)), GLM1_CBD_cohort_quad, GLM2_CBD_cohort_quad, GLM3_CBD_cohort_quad, GLM1_name = "CBD_cohort_quad_Poisson_log" , GLM2_name = "CBD_cohort_quad_Binomial_cloglog" , GLM3_name = "CBD_cohort_quad_Binomial_logit" , main_name = "CBD_cohort_quad", y_name = "Year_coef_Kappa_2", pos_legend = 'topleft')

Plot_Superpose_Coef(TheRange = c((n_y + n_y+ 1): (n_y + n_y + n_y)), GLM1_CBD_cohort_quad, GLM2_CBD_cohort_quad, GLM3_CBD_cohort_quad, GLM1_name = "CBD_cohort_quad_Poisson_log" , GLM2_name = "CBD_cohort_quad_Binomial_cloglog" , GLM3_name = "CBD_cohort_quad_Binomial_logit" , main_name = "CBD_cohort_quad", y_name = "Year_coef_Kappa_3", pos_legend = "topleft")

Plot_Superpose_Coef(TheRange = c((n_y + n_y + n_y + 1): (n_y + n_y +  n_y + n_c)), GLM1_CBD_cohort_quad, GLM2_CBD_cohort_quad, GLM3_CBD_cohort_quad, GLM1_name = "CBD_cohort_quad_Poisson_log" , GLM2_name = "CBD_cohort_quad_Binomial_cloglog" , GLM3_name = "CBD_cohort_quad_Binomial_logit" , main_name = "CBD_cohort_quad", y_name = "Cohort_coef_Gamma", pos_legend = "topleft")


    #############################################
    ### Modèle CBD  + modulated cohort effect ###   Done
    #############################################


Model_Matrix_CBD_modulated_cohort = function(n_a, n_y, delta) 
{
  #Création X1
  Identity_n_y = matrix(0, ncol = n_y, nrow = n_y, byrow = TRUE)
  diag(Identity_n_y) = 1
  serie1_n_a = rep(1, n_a)
  X_1 = matrix(0 , ncol = n_y, nrow = n_a * n_y, byrow = TRUE )
  
  for (j in 1:n_y)
  {
    X_1[(1:n_a) + n_a*(j-1) ,j] = 1
  }
  
  #Création X2
  Centred_Age_Vector = (50:90) - mean(50:90)
  X_2 = matrix(0 , ncol = n_y, nrow = n_a * n_y, byrow = TRUE )
  
  for (j in 1:n_y)
  {
    for( i in 1:length(Centred_Age_Vector))
    {
      X_2[i + (n_a) * (j-1),j] = Centred_Age_Vector[i]
    }
  }
  
  
  
  ### Rework de X_c! 
  n_c = n_a + n_y -1
  X_c = matrix(0 , ncol = n_c, nrow = n_a * n_y, byrow = TRUE )
  
  count = 1
  for(i in (0:(n_y-1)))
  {
    for(j in (0:(n_a-1)))
    {
      vect = rep(0,n_c)
      vect[n_a - j +i] = 1
      
      X_c[count,] = vect
      count = count +1 
    }
    
  }
  
  # Création de v(delta)
  
  v_delta = rep(delta* serie1_n_a - (50:90), n_y)

  X = cbind(X_1, X_2 ,v_delta*X_c)
  
  return (X)
  
}  
delta = 80
X_CBD_modulated_cohort= Model_Matrix_CBD_modulated_cohort(n_a, n_y, delta)

GLM1_CBD_modulated_cohort = glm(dx ~ -1 +  X_CBD_modulated_cohort + offset(log(CentralExpo)), family = poisson(link = "log"))
GLM2_CBD_modulated_cohort = glm(Qx ~ -1 + X_CBD_modulated_cohort, family = binomial(link = "cloglog"), weight = InitialExpo)
GLM3_CBD_modulated_cohort = glm(Qx ~ -1 + X_CBD_modulated_cohort, family = binomial(link = "logit"), weight = InitialExpo)

CBD_modulated_cohort_fixed_delta = c(GLM1_CBD_modulated_cohort$deviance, GLM2_CBD_modulated_cohort$deviance, GLM3_CBD_modulated_cohort$deviance)
CBD_modulated_cohort_fixed_delta


CBD_modulated_cohort_fixed_delta_dev_Poisson = c(GLM1_CBD_modulated_cohort$deviance,
                                          Deviance.FromBinomial(dx_obs = dx, qx_fit = GLM2_CBD_modulated_cohort$fitted.values, initial_ETR = InitialExpo, devType = "Poisson"),
                                          Deviance.FromBinomial(dx_obs = dx, qx_fit = GLM3_CBD_modulated_cohort$fitted.values, initial_ETR = InitialExpo, devType = "Poisson"))



NumberConstraints_CBD_modulated_cohort = dim(X_CBD_modulated_cohort)[2] - qr(X_CBD_modulated_cohort)$rank
NumberConstraints_CBD_modulated_cohort

Plotting_modulated_cohort = function(range, arg )  
{
  CBD_modulated_cohort_vector = c()
  serie1_n_a = rep(1, n_a)
  
  for (i in (range[1]:range[length(range)]))
  {
    delta = i
    v_delta = rep(delta* serie1_n_a - (50:90), n_a)
    X_CBD_modulated_cohort= Model_Matrix_CBD_modulated_cohort(n_a, n_y, delta)
    if (arg == "log")
    {
      GLM_CBD_modulated_cohort = glm(dx ~ -1 +  X_CBD_modulated_cohort + offset(log(CentralExpo)), family = poisson(link = "log"))
    }
    else if (arg == "cloglog")
    {
      GLM_CBD_modulated_cohort = glm(Qx ~ -1 + X_CBD_modulated_cohort, family = binomial(link = "cloglog"), weight = InitialExpo)
    }
    else if (arg == "logit")
    {
      GLM_CBD_modulated_cohort = glm(Qx ~ -1 + X_CBD_modulated_cohort, family = binomial(link = "logit"), weight = InitialExpo)
    }
    CBD_modulated_cohort_vector = c(CBD_modulated_cohort_vector, GLM_CBD_modulated_cohort$deviance)
  }
  min = min(CBD_modulated_cohort_vector)
  arg_min = which(CBD_modulated_cohort_vector == min) + range[1]
  plot(c(range[1]:range[length(range)]), CBD_modulated_cohort_vector, type = "l", ylab = "Deviance", xlab = "Constant delta", main = paste("Deviance of glm model with ",arg, " link"))
  return (arg_min)
}
delta = Plotting_modulated_cohort(c(-200:0), "log")
delta = Plotting_modulated_cohort(c(0:50), "log")
delta_log = Plotting_modulated_cohort(c(50:90), "log")
#delta = Plotting_modulated_cohort(c(90:140), "log")
#delta = Plotting_modulated_cohort(c(140:250), "log")
#delta = Plotting_modulated_cohort(c(250:400), "log")
delta = Plotting_modulated_cohort(c(-200:90), "log")
#delta = Plotting_modulated_cohort(c(50:55), "log")

delta_cloglog = Plotting_modulated_cohort(c(50:90), "cloglog")
delta_logit = Plotting_modulated_cohort(c(50:90), "logit")


#Remplace les NA que R a mis pour les contraintes par des 0
GLM1_CBD_modulated_cohort$coefficients[208] = 0

GLM2_CBD_modulated_cohort$coefficients[208] = 0

GLM3_CBD_modulated_cohort$coefficients[208] = 0

#Calcul des paramètres cohérents
weight = c()
for(i in 1:n_a)
{
  weight[i] = i
}
for (i in 1:(n_y -n_a))
{
  weight[i + n_a] = n_a
}
for(i in 1:(n_a-1))
{
  weight[i + n_y] = n_a -i
}

Kappa_1_t_CBD_modulated_cohort_1 = GLM3_CBD_modulated_cohort$coefficients[1:n_y]
Kappa_2_t_CBD_modulated_cohort_1 = GLM3_CBD_modulated_cohort$coefficients[(n_y +1):(2*n_y)]
Gamma_t_CBD_modulated_cohort_1 = GLM3_CBD_modulated_cohort$coefficients[(2*n_y +1):(length(GLM1_CBD_modulated_cohort$coefficients))]

Mean_Gamma = (weight %*% Gamma_t_CBD_modulated_cohort_1) / sum(weight)

Gamma_t_CBD_modulated_cohort_1 = Gamma_t_CBD_modulated_cohort_1 - Mean_Gamma * rep(1,n_c)
Kappa_2_t_CBD_modulated_cohort_1 = Kappa_2_t_CBD_modulated_cohort_1 - Mean_Gamma * rep(1,n_y)
Kappa_1_t_CBD_modulated_cohort_1 = Kappa_1_t_CBD_modulated_cohort_1 + Mean_Gamma * (delta - mean(c(50:90)) )* rep(1,n_y)
teta_CBD_modulated_cohort_1 = c(Kappa_1_t_CBD_modulated_cohort_1, Kappa_2_t_CBD_modulated_cohort_1, Gamma_t_CBD_modulated_cohort_1)



## then with Poisson assumption, we find :
    ##For the log link
eta_CBD_modulated_cohort_1 = (X_CBD_modulated_cohort %*% teta_CBD_modulated_cohort_1)
mu_log_CBD_modulated_cohort_1 = exp(eta_CBD_modulated_cohort_1)
mu_log_CBD_modulated_cohort_1 = matrix(mu_log_CBD_modulated_cohort_1, nrow = n_a, ncol = n_y, byrow = FALSE)  #Mat format
persp(Age,Year, mu_log_CBD_modulated_cohort_1, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "CBD_modulated_cohort_log")

## then with Binomial assumption, we find :
    ##For the cloglog link
eta_CBD_modulated_cohort_2 = (X_CBD_modulated_cohort %*% GLM2_CBD_modulated_cohort$coefficients)
qx_cloglog_CBD_modulated_cohort_2 = 1-exp(-exp(eta_CBD_modulated_cohort_2))
qx_cloglog_CBD_modulated_cohort_2 = matrix(qx_cloglog_CBD_modulated_cohort_2, nrow = n_a, ncol = n_y, byrow = FALSE)  #Mat format
mu_cloglog_CBD_modulated_cohort_2 = -log(1-qx_cloglog_CBD_modulated_cohort_2)
persp(Age,Year, mu_cloglog_CBD_modulated_cohort_2, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "CBD_modulated_cohort_cloglog")

    ##For the logit link
eta_CBD_modulated_cohort_3 = (X_CBD_modulated_cohort %*% GLM3_CBD_modulated_cohort$coefficients)
qx_logit_CBD_modulated_cohort_3 = (exp(eta_CBD_modulated_cohort_3))/(1 + exp(eta_CBD_modulated_cohort_3))
qx_logit_CBD_modulated_cohort_3 = matrix(qx_logit_CBD_modulated_cohort_3, nrow = n_a, ncol = n_y, byrow = FALSE)  #Mat format
mu_logit_CBD_modulated_cohort_3 = -log(1-qx_logit_CBD_modulated_cohort_3)
persp(Age,Year, mu_logit_CBD_modulated_cohort_3, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "CBD_modulated_cohort_logit")

CBD_modulated_cohort_aic = c(GLM1_CBD_modulated_cohort$aic, GLM2_CBD_modulated_cohort$aic, GLM3_CBD_modulated_cohort$aic)
CBD_modulated_cohort_aic

Resids_Plot(residuals_mu = matrix(c(rstandard(GLM1_CBD_modulated_cohort)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "CBD_modulated_cohort_poisson_log" )
Resids_Plot(residuals_mu = matrix(c(rstandard(GLM2_CBD_modulated_cohort)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "CBD_modulated_cohort_Binomial_cloglog" )
Resids_Plot(residuals_mu = matrix(c(rstandard(GLM3_CBD_modulated_cohort)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "CBD_modulated_cohort_Binomial_logit" )


      ###########################################################
      ### Modèle CBD  + modulated cohort effect Optimal Delta ###   Done
      ###########################################################


OptimDelta_GLM1 = function(delta)
{
  X_CBD_modulated_cohort = Model_Matrix_CBD_modulated_cohort(n_a, n_y, delta[1])
  glm(dx ~ -1 +  X_CBD_modulated_cohort + offset(log(CentralExpo)), family = poisson(link = "log"))$deviance
}
OptimDelta_GLM2 = function(delta)
{
  X_CBD_modulated_cohort = Model_Matrix_CBD_modulated_cohort(n_a, n_y, delta[1])
  glm(Qx ~ -1 + X_CBD_modulated_cohort, family = binomial(link = "cloglog"), weight = InitialExpo)$deviance
}
OptimDelta_GLM3 = function(delta)
{
  X_CBD_modulated_cohort = Model_Matrix_CBD_modulated_cohort(n_a, n_y, delta[1])
  glm(Qx ~ -1 + X_CBD_modulated_cohort, family = binomial(link = "logit"), weight = InitialExpo)$deviance
}

resultatOptim_GLM1 = nlm(OptimDelta_GLM1, p = c(10), hessian = TRUE)
resultatOptim_GLM2 = nlm(OptimDelta_GLM2, p = c(10), hessian = TRUE)
resultatOptim_GLM3 = nlm(OptimDelta_GLM3, p = c(10), hessian = TRUE)


CBD_modulated_cohort_optim = c(resultatOptim_GLM1$minimum, resultatOptim_GLM2$minimum, resultatOptim_GLM3$minimum)
CBD_modulated_cohort_optim_delta = c(resultatOptim_GLM1$estimate, resultatOptim_GLM2$estimate, resultatOptim_GLM3$estimate)



## Calcul des nouveaux paramètres sur base des deltas optimaux
X_CBD_modulated_cohort= Model_Matrix_CBD_modulated_cohort(n_a, n_y, CBD_modulated_cohort_optim_delta[1])
GLM1_CBD_modulated_cohort_optim = glm(dx ~ -1 +  X_CBD_modulated_cohort + offset(log(CentralExpo)), family = poisson(link = "log"))

#Remplace les NA que R a mis pour les contraintes par des 0
GLM1_CBD_modulated_cohort_optim$coefficients[208] = 0

## then with Poisson assumption, we find :
    ##For the log link
eta_CBD_modulated_cohort_optim_1 = (X_CBD_modulated_cohort %*% GLM1_CBD_modulated_cohort_optim$coefficients)
mu_log_CBD_modulated_cohort_optim_1 = exp(eta_CBD_modulated_cohort_optim_1)
mu_log_CBD_modulated_cohort_optim_1 = matrix(mu_log_CBD_modulated_cohort_optim_1, nrow = n_a, ncol = n_y, byrow = FALSE)  #Mat format
persp(Age,Year, mu_log_CBD_modulated_cohort_optim_1, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "CBD_modulated_cohort_optim_log")



## Calcul des nouveaux paramètres sur base des deltas optimaux
X_CBD_modulated_cohort= Model_Matrix_CBD_modulated_cohort(n_a, n_y, CBD_modulated_cohort_optim_delta[2])
GLM2_CBD_modulated_cohort_optim = glm(Qx ~ -1 + X_CBD_modulated_cohort, family = binomial(link = "cloglog"), weight = InitialExpo)

#Remplace les NA que R a mis pour les contraintes par des 0
GLM2_CBD_modulated_cohort_optim$coefficients[208] = 0

## then with Binomial assumption, we find :
    ##For the cloglog link
eta_CBD_modulated_cohort_optim_2 = (X_CBD_modulated_cohort %*% GLM2_CBD_modulated_cohort_optim$coefficients)
qx_cloglog_CBD_modulated_cohort_optim_2 = 1-exp(-exp(eta_CBD_modulated_cohort_optim_2))
qx_cloglog_CBD_modulated_cohort_optim_2 = matrix(qx_cloglog_CBD_modulated_cohort_optim_2, nrow = n_a, ncol = n_y, byrow = FALSE)  #Mat format
mu_cloglog_CBD_modulated_cohort_optim_2 = -log(1-qx_cloglog_CBD_modulated_cohort_optim_2)
persp(Age,Year, mu_cloglog_CBD_modulated_cohort_optim_2, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "CBD_modulated_cohort_optim_cloglog")




## Calcul des nouveaux paramètres sur base des deltas optimaux
X_CBD_modulated_cohort= Model_Matrix_CBD_modulated_cohort(n_a, n_y, CBD_modulated_cohort_optim_delta[3])
GLM3_CBD_modulated_cohort_optim = glm(Qx ~ -1 + X_CBD_modulated_cohort, family = binomial(link = "logit"), weight = InitialExpo)

#Remplace les NA que R a mis pour les contraintes par des 0
GLM3_CBD_modulated_cohort_optim$coefficients[208] = 0

  ##For the logit link
eta_CBD_modulated_cohort_optim_3 = (X_CBD_modulated_cohort %*% GLM3_CBD_modulated_cohort_optim$coefficients)
qx_logit_CBD_modulated_cohort_optim_3 = (exp(eta_CBD_modulated_cohort_optim_3))/(1 + exp(eta_CBD_modulated_cohort_optim_3))
qx_logit_CBD_modulated_cohort_optim_3 = matrix(qx_logit_CBD_modulated_cohort_optim_3, nrow = n_a, ncol = n_y, byrow = FALSE)  #Mat format
mu_logit_CBD_modulated_cohort_optim_3 = -log(1-qx_logit_CBD_modulated_cohort_optim_3)
persp(Age,Year, mu_logit_CBD_modulated_cohort_optim_3, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "CBD_modulated_cohort_optim_logit")




CBD_modulated_cohort_optim_delta_dev_Poisson = c(GLM1_CBD_modulated_cohort_optim$deviance,
                                                 Deviance.FromBinomial(dx_obs = dx, qx_fit = GLM2_CBD_modulated_cohort_optim$fitted.values, initial_ETR = InitialExpo, devType = "Poisson"),
                                                 Deviance.FromBinomial(dx_obs = dx, qx_fit = GLM3_CBD_modulated_cohort_optim$fitted.values, initial_ETR = InitialExpo, devType = "Poisson"))

CBD_modulated_cohort_optim_aic = c(GLM1_CBD_modulated_cohort_optim$aic, GLM2_CBD_modulated_cohort_optim$aic, GLM3_CBD_modulated_cohort_optim$aic)
CBD_modulated_cohort_optim_aic


Resids_Plot(residuals_mu = matrix(c(rstandard(GLM1_CBD_modulated_cohort_optim)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "CBD_modulated_cohort_optim_poisson_log" )
Resids_Plot(residuals_mu = matrix(c(rstandard(GLM2_CBD_modulated_cohort_optim)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "CBD_modulated_cohort_optim_Binomial_cloglog" )
Resids_Plot(residuals_mu = matrix(c(rstandard(GLM3_CBD_modulated_cohort_optim)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "CBD_modulated_cohort_optim_Binomial_logit" )


#####################################################
#####     Creating plot functions for parameters#####   Done
#####################################################

Changings_Parameters = function(glm , delta)
{
  
  Model_Matrix_CBD_modulated_cohort = function(n_a, n_y, delta) 
  {
    #Création X1
    Identity_n_y = matrix(0, ncol = n_y, nrow = n_y, byrow = TRUE)
    diag(Identity_n_y) = 1
    serie1_n_a = rep(1, n_a)
    X_1 = matrix(0 , ncol = n_y, nrow = n_a * n_y, byrow = TRUE )
    
    for (j in 1:n_y)
    {
      X_1[(1:n_a) + n_a*(j-1) ,j] = 1
    }
    
    #Création X2
    Centred_Age_Vector = (50:90) - mean(50:90)
    X_2 = matrix(0 , ncol = n_y, nrow = n_a * n_y, byrow = TRUE )
    
    for (j in 1:n_y)
    {
      for( i in 1:length(Centred_Age_Vector))
      {
        X_2[i + (n_a) * (j-1),j] = Centred_Age_Vector[i]
      }
    }
    
    
    ### Rework de X_c! 
    n_c = n_a + n_y -1
    X_c = matrix(0 , ncol = n_c, nrow = n_a * n_y, byrow = TRUE )
    
    count = 1
    for(i in (0:(n_y-1)))
    {
      for(j in (0:(n_a-1)))
      {
        vect = rep(0,n_c)
        vect[n_a - j +i] = 1
        
        X_c[count,] = vect
        count = count +1 
      }
      
    }
    
    # Création de v(delta)
    
    v_delta = rep(delta* serie1_n_a - (50:90), n_y)
    
    X = cbind(X_1, X_2 ,v_delta*X_c)
    
    return (X)
    
  }  
  delta = 50
  X_CBD_modulated_cohort= Model_Matrix_CBD_modulated_cohort(n_a, n_y, delta)
  GLM1_CBD_modulated_cohort_1 = glm(dx ~ -1 +  X_CBD_modulated_cohort + offset(log(CentralExpo)), family = poisson(link = "log"))
  
  
  Model_Matrix_CBD_modulated_cohort = function(n_a, n_y, delta) 
  {
    #Création X1
    Identity_n_y = matrix(0, ncol = n_y, nrow = n_y, byrow = TRUE)
    diag(Identity_n_y) = 1
    serie1_n_a = rep(1, n_a)
    X_1 = matrix(0 , ncol = n_y, nrow = n_a * n_y, byrow = TRUE )
    
    for (j in 1:n_y)
    {
      X_1[(1:n_a) + n_a*(j-1) ,j] = 1
    }
    
    #Création X2
    Centred_Age_Vector = (50:90) - mean(50:90)
    X_2 = matrix(0 , ncol = n_y, nrow = n_a * n_y, byrow = TRUE )
    
    for (j in 1:n_y)
    {
      for( i in 1:length(Centred_Age_Vector))
      {
        X_2[i + (n_a) * (j-1),j] = Centred_Age_Vector[i]
      }
    }
    
    
    # Cohort matrix
    #n_c = n_a + n_y 
    #X_c = matrix(0 , ncol = n_c, nrow = n, byrow = TRUE )
    
    #  for (age in 1:n_a)
    # {
    #  for( j in 1:n_c)
    # {
    #  if((n_a + j)  <= 97){
    #   X_c[ j + n_y * (age-1) , n_a + (j-1) - (age-1)] = 1 
    # }
    
    #  }
    #}
    
    ### Rework de X_c! 
    n_c = n_a + n_y -1
    X_c = matrix(0 , ncol = n_c, nrow = n_a * n_y, byrow = TRUE )
    
    count = 1
    for(i in (0:(n_y-1)))
    {
      for(j in (0:(n_a-1)))
      {
        vect = rep(0,n_c)
        vect[n_a - j +i] = 1
        
        X_c[count,] = vect
        count = count +1 
      }
      
    }
    
    # Création de v(delta)
    
    v_delta = rep(delta* serie1_n_a - (50:90), n_y)
    
    X = cbind(X_1, X_2 ,v_delta*X_c)
    
    return (X)
    
  }  
  delta = 90
  X_CBD_modulated_cohort= Model_Matrix_CBD_modulated_cohort(n_a, n_y, delta)
  GLM1_CBD_modulated_cohort_2 = glm(dx ~ -1 +  X_CBD_modulated_cohort + offset(log(CentralExpo)), family = poisson(link = "log"))
  
  Model_Matrix_CBD_modulated_cohort = function(n_a, n_y, delta) 
  {
    #Création X1
    Identity_n_y = matrix(0, ncol = n_y, nrow = n_y, byrow = TRUE)
    diag(Identity_n_y) = 1
    serie1_n_a = rep(1, n_a)
    X_1 = matrix(0 , ncol = n_y, nrow = n_a * n_y, byrow = TRUE )
    
    for (j in 1:n_y)
    {
      X_1[(1:n_a) + n_a*(j-1) ,j] = 1
    }
    
    #Création X2
    Centred_Age_Vector = (50:90) - mean(50:90)
    X_2 = matrix(0 , ncol = n_y, nrow = n_a * n_y, byrow = TRUE )
    
    for (j in 1:n_y)
    {
      for( i in 1:length(Centred_Age_Vector))
      {
        X_2[i + (n_a) * (j-1),j] = Centred_Age_Vector[i]
      }
    }
    
    
    # Cohort matrix
    #n_c = n_a + n_y 
    #X_c = matrix(0 , ncol = n_c, nrow = n, byrow = TRUE )
    
    #  for (age in 1:n_a)
    # {
    #  for( j in 1:n_c)
    # {
    #  if((n_a + j)  <= 97){
    #   X_c[ j + n_y * (age-1) , n_a + (j-1) - (age-1)] = 1 
    # }
    
    #  }
    #}
    
    ### Rework de X_c! 
    n_c = n_a + n_y -1
    X_c = matrix(0 , ncol = n_c, nrow = n_a * n_y, byrow = TRUE )
    
    count = 1
    for(i in (0:(n_y-1)))
    {
      for(j in (0:(n_a-1)))
      {
        vect = rep(0,n_c)
        vect[n_a - j +i] = 1
        
        X_c[count,] = vect
        count = count +1 
      }
      
    }
    
    # Création de v(delta)
    
    v_delta = rep(delta* serie1_n_a - (50:90), n_y)
    
    X = cbind(X_1, X_2 ,v_delta*X_c)
    
    return (X)
    
  }  
  delta = resultatOptim_GLM1$estimate
  X_CBD_modulated_cohort= Model_Matrix_CBD_modulated_cohort(n_a, n_y, delta)
  GLM1_CBD_modulated_cohort_3 = glm(dx ~ -1 +  X_CBD_modulated_cohort + offset(log(CentralExpo)), family = poisson(link = "log"))
  #Remplace les NA que R a mis pour les contraintes par des 0
  GLM1_CBD_modulated_cohort_1$coefficients[207] = 0
  GLM1_CBD_modulated_cohort_1$coefficients[208] = 0
  GLM1_CBD_modulated_cohort_1$coefficients[113] = 0
  
  GLM1_CBD_modulated_cohort_2$coefficients[207] = 0
  GLM1_CBD_modulated_cohort_2$coefficients[208] = 0
  GLM1_CBD_modulated_cohort_2$coefficients[113] = 0
  
  GLM1_CBD_modulated_cohort_3$coefficients[207] = 0  
  GLM1_CBD_modulated_cohort_3$coefficients[208] = 0
  GLM1_CBD_modulated_cohort_3$coefficients[113] = 0
  
  plot(c(1960:2015), GLM1_CBD_modulated_cohort_3$coefficients[1:n_y], xlab = "year", ylab = "Kappa_1",ylim = c(min(GLM1_CBD_modulated_cohort_1$coefficients[1:n_y],
                                                                                                                   GLM1_CBD_modulated_cohort_2$coefficients[1:n_y],
                                                                                                                   GLM1_CBD_modulated_cohort_3$coefficients[1:n_y]) ,
                                                                                                               max(GLM1_CBD_modulated_cohort_1$coefficients[1:n_y],
                                                                                                                   GLM1_CBD_modulated_cohort_2$coefficients[1:n_y],
                                                                                                                   GLM1_CBD_modulated_cohort_3$coefficients[1:n_y])), type = "l")
  lines(c(1960:2015), GLM1_CBD_modulated_cohort_2$coefficients[1:n_y], col = "red")
  lines(c(1960:2015), GLM1_CBD_modulated_cohort_1$coefficients[1:n_y], col = "green")
  legend("topright", c(paste("delta = ", resultatOptim_GLM1$estimate),"delta = 90", "delta = 50"), col = c("black", "red", "green"), lty = 1)
  
  
  plot(c(1960:2015), GLM1_CBD_modulated_cohort_3$coefficients[(n_y +1) :(2*n_y)],xlab = "year", ylab = "Kappa_2", ylim = c(min(GLM1_CBD_modulated_cohort_1$coefficients[(n_y +1) :(2*n_y)],
                                                                                                                               GLM1_CBD_modulated_cohort_2$coefficients[(n_y +1) :(2*n_y)],
                                                                                                                               GLM1_CBD_modulated_cohort_3$coefficients[(n_y +1) :(2*n_y)]) ,
                                                                                                                           max(GLM1_CBD_modulated_cohort_1$coefficients[(n_y +1) :(2*n_y)],
                                                                                                                               GLM1_CBD_modulated_cohort_2$coefficients[(n_y +1) :(2*n_y)],
                                                                                                                               GLM1_CBD_modulated_cohort_3$coefficients[(n_y +1) :(2*n_y)])), type = "l")
  lines(c(1960:2015), GLM1_CBD_modulated_cohort_2$coefficients[(n_y +1) :(2*n_y)], col = "red")
  lines(c(1960:2015), GLM1_CBD_modulated_cohort_1$coefficients[(n_y +1) :(2*n_y)], col = "green")
  legend("right", c(paste("delta = ", resultatOptim_GLM1$estimate),"delta = 90", "delta = 50"), col = c("black", "red", "green"), lty = 1)
  
  
  
  plot(c((2015-n_c +2):2015), GLM1_CBD_modulated_cohort_3$coefficients[(2*n_y +2) :(2*n_y + n_c)],xlab = "year of birth", ylab = "Gamma" ,ylim = c(min(GLM1_CBD_modulated_cohort_1$coefficients[(2*n_y +2) :(2*n_y + n_c)],
                                                                                                                                                       GLM1_CBD_modulated_cohort_2$coefficients[(2*n_y +2) :(2*n_y + n_c)],
                                                                                                                                                       GLM1_CBD_modulated_cohort_3$coefficients[(2*n_y +2) :(2*n_y + n_c)]) ,
                                                                                                                                                   max(GLM1_CBD_modulated_cohort_1$coefficients[(2*n_y +2) :(2*n_y + n_c)],
                                                                                                                                                       GLM1_CBD_modulated_cohort_2$coefficients[(2*n_y +2) :(2*n_y + n_c)],
                                                                                                                                                       GLM1_CBD_modulated_cohort_3$coefficients[(2*n_y +2) :(2*n_y + n_c)])), type = "l")
  lines(c((2015-n_c +1 ):2015), GLM1_CBD_modulated_cohort_2$coefficients[(2*n_y +1) :(2*n_y + n_c)], col = "red")
  lines(c((2015-n_c +2 ):2015), GLM1_CBD_modulated_cohort_1$coefficients[(2*n_y +2) :(2*n_y + n_c)], col = "green")
  legend("topright", c(paste("delta = ", resultatOptim_GLM1$estimate),"delta = 90", "delta = 50"), col = c("black", "red", "green"), lty = 1)
  
}






      #####################################################
      #####     Recupatilatif des Modèles via glm     #####   Done
      #####################################################

tab_GLM = data.frame(APC , CBD, CBD_cohort, CBD_cohort_quad, CBD_modulated_cohort_fixed_delta, CBD_modulated_cohort_optim)
rownames(tab_GLM) = c("GLM1_log", "GLM2_cloglog", "GLM3_logit")

tab_optim_delta = data.frame(CBD_modulated_cohort_optim,CBD_modulated_cohort_optim_delta)
rownames(tab_optim_delta) = c("GLM1_log", "GLM2_cloglog", "GLM3_logit")
colnames(tab_optim_delta) = c("CBD_modulated_cohort_optim_deviance", "Local_optimal_delta")

tab_constraints = data.frame(NumberConstraints_APC,NumberConstraints_CBD, NumberConstraints_CBD_cohort, NumberConstraints_CBD_cohort_quad ,NumberConstraints_CBD_modulated_cohort)
rownames(tab_constraints) = c("Number of constraints: ")
colnames(tab_constraints) = c("APC", "CBD" ,"CBD_c", "CBD_c_q", "CBD_m_c")

tab_aic = data.frame(APC_aic, CBD_aic, CBD_cohort_aic, CBD_cohort_quad_aic, CBD_modulated_cohort_aic, CBD_modulated_cohort_optim_aic)
rownames(tab_constraints) = c("AIC criterion: ")
colnames(tab_constraints) = c("APC", "CBD" ,"CBD_c", "CBD_c_q", "CBD_m_c")

tab_GLM
tab_optim_delta
tab_constraints
tab_aic








#################################################################################
######    Ecriture des modèles M1/M2/M4 avec leurs applications via glm    ######
#################################################################################


      ##############################
      #### Modèle de Lee-Carter ####    Done
      ##############################

# As it has been described before



# The vector of age suffices to the dx vector
X_suff_d = rep(c(50:90), n_y)
# The vector of year suffices to the dx vector: y kronenker product 1_n_a
Y_suff_d = c() 
years = c(1960:2015)
for (i in years)
{
  Y_suff_d = c(Y_suff_d , rep(i,n_a))
}

### Converting variables to factors

Age_F = factor(X_suff_d)
Year_F = factor(Y_suff_d)

set.seed(1)
GNM1_Lee_Carter = gnm(dx ~ -1 + Age_F + Mult(Age_F, Year_F) + offset(log(CentralExpo)), family = poisson(link = "log"))
set.seed(1)
GNM2_Lee_Carter = gnm(Qx ~ -1 + Age_F + Mult(Age_F, Year_F), family = binomial(link = "cloglog"), weight = InitialExpo)
set.seed(1)
GNM3_Lee_Carter = gnm(Qx ~ -1 + Age_F + Mult(Age_F, Year_F), family = binomial(link = "logit"), weight = InitialExpo)

Lee_Carter = c(GNM1_Lee_Carter$deviance, GNM2_Lee_Carter$deviance, GNM3_Lee_Carter$deviance)


Lee_Carter_dev_Poisson = c(GNM1_Lee_Carter$deviance,
                           Deviance.FromBinomial(dx_obs = dx, qx_fit = GNM2_Lee_Carter$fitted.values, initial_ETR = InitialExpo, devType = "Poisson"),
                           Deviance.FromBinomial(dx_obs = dx, qx_fit = GNM3_Lee_Carter$fitted.values, initial_ETR = InitialExpo, devType = "Poisson"))

## with the Poisson assumption, we find
    # log link
eta_LC_1 = c()
for( i in ((2*n_a +1) : (length(GNM1_Lee_Carter$coefficients))))
{
  eta_LC_1 = c(eta_LC_1, GNM1_Lee_Carter$coefficients[1:n_a] + 
    GNM1_Lee_Carter$coefficients[(n_a+1): (n_a +n_a)] * GNM1_Lee_Carter$coefficients[i])
}
eta_LC_1 = matrix(eta_LC_1, nrow = n_a, ncol = n_y, byrow = FALSE)
mu_log_LC_1 = exp(eta_LC_1)
persp(Age,Year, mu_log_LC_1, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "Lee_Carter_log")

## with the Binomial assumption, we find
    # cloglog link
eta_LC_2 = c()
for( i in ((2*n_a +1) : (length(GNM2_Lee_Carter$coefficients))))
{
  eta_LC_2 = c(eta_LC_2, GNM2_Lee_Carter$coefficients[1:n_a] + 
                 GNM2_Lee_Carter$coefficients[(n_a+1): (n_a +n_a)] * GNM2_Lee_Carter$coefficients[i])
}
eta_LC_2 = matrix(eta_LC_2, nrow = n_a, ncol = n_y, byrow = FALSE)
qx_cloglog_LC_2 = 1-exp(-exp(eta_LC_2))
mu_cloglog_LC_2 = -log(1-qx_cloglog_LC_2)
persp(Age,Year, mu_cloglog_LC_2, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "Lee_Carter_cloglog")

    # logit link
eta_LC_3 = c()
for( i in ((2*n_a +1) : (length(GNM3_Lee_Carter$coefficients))))
{
  eta_LC_3 = c(eta_LC_3, GNM3_Lee_Carter$coefficients[1:n_a] + 
                 GNM3_Lee_Carter$coefficients[(n_a+1): (n_a +n_a)] * GNM3_Lee_Carter$coefficients[i])
}
eta_LC_3 = matrix(eta_LC_3, nrow = n_a, ncol = n_y, byrow = FALSE)
qx_logit_LC_3 = (exp(eta_LC_3))/(1 + exp(eta_LC_3))
mu_logit_LC_3 = -log(1-qx_logit_LC_3)
persp(Age,Year, mu_logit_LC_3, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "Lee_Carter_logit")


Resids_Plot(residuals_mu = matrix(c(rstandard(GNM1_Lee_Carter)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "Lee_Carter_poisson_log" )
Resids_Plot(residuals_mu = matrix(c(rstandard(GNM2_Lee_Carter)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "Lee_Carter_Binomial_cloglog" )
Resids_Plot(residuals_mu = matrix(c(rstandard(GNM3_Lee_Carter)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "Lee_Carter_Binomial_logit" )




Plot_Superpose_Coef(TheRange = c(1:n_a), GNM1_Lee_Carter, GNM2_Lee_Carter, GNM3_Lee_Carter, GLM1_name = "Lee_Carter_Poisson_log" , GLM2_name = "Lee_Carter_Binomial_cloglog" , GLM3_name = "Lee_Carter_Binomial_logit" , main_name = "Lee_Carter", y_name = "Age_coef_Beta_1", pos_legend = "topleft")

Plot_Superpose_Coef(TheRange = c((n_a +1):( 2* n_a)), GNM1_Lee_Carter, GNM2_Lee_Carter, GNM3_Lee_Carter, GLM1_name = "Lee_Carter_Poisson_log" , GLM2_name = "Lee_Carter_Binomial_cloglog" , GLM3_name = "Lee_Carter_Binomial_logit" , main_name = "Lee_Carter", y_name = "Age_coef_Beta_2", pos_legend = 'topright')

Plot_Superpose_Coef(TheRange = c((2* n_a + 1): (2* n_a + n_y)), GNM1_Lee_Carter, GNM2_Lee_Carter, GNM3_Lee_Carter, GLM1_name = "Lee_Carter_Poisson_log" , GLM2_name = "Lee_Carter_Binomial_cloglog" , GLM3_name = "Lee_Carter_Binomial_logit" , main_name = "Lee_Carter", y_name = "Year_coef_Kappa_2", pos_legend = "topright")

Plot_Superpose_Coef(TheRange = c((n_y + n_y + n_y + 1): (n_y + n_y +  n_y + n_c)), GLM1_Lee_Carter, GLM2_Lee_Carter, GLM3_Lee_Carter, GLM1_name = "Lee_Carter_Poisson_log" , GLM2_name = "Lee_Carter_Binomial_cloglog" , GLM3_name = "Lee_Carter_Binomial_logit" , main_name = "Lee_Carter", y_name = "Cohort_coef_Gamma", pos_legend = "topleft")



      ####################################
      #### Modèle de Renshaw-Haberman ####    Done    
      ####################################

### matrice X_c! 
n_c = n_a + n_y -1
X_c = matrix(0 , ncol = n_c, nrow = n, byrow = TRUE )

count = 1
for(i in (0:(n_y-1)))
{
  for(j in (0:(n_a-1)))
  {
    vect = rep(0,n_c)
    vect[n_a - j +i] = 1
    
    X_c[count,] = vect
    count = count +1 
  }
  
}


z_vector = c() 
for (i in (1:n))
{
  z_vector = c(z_vector , which(X_c[i,] == 1))
}
Cohort_F = factor(z_vector)


### Les GNM suivant convergent si j'augmente le nombre d'iterations :-)
#La pramétrisation étant aléatoire, cela complique grandement les forecasts... Faire run jusqu'à avoir un bon set de paramètres à forecaster
#Remarquons que l'algorithme ne s'arrête car uniquement arrivé a itermax
conv_count_log = 0
conv_count_cloglog = 0
conv_count_logit = 0
Dev_log = c()
Dev_cloglog = c()
Dev_logit = c()
iter_log = c()
iter_cloglog = c()
iter_logit = c()
for(i in (1:10))
{
  print(i)
  GNM1_Renshaw_Haberman = gnm(dx ~ -1 + Age_F + Mult(Age_F, Year_F) + Mult(Age_F, Cohort_F) + offset(log(CentralExpo)), family = poisson(link = "log"), iterMax = 1500)

  if((is.null(GNM1_Renshaw_Haberman) == FALSE))
  {
    if((GNM1_Renshaw_Haberman$converged == TRUE) )
    {
      conv_count_log = conv_count_log + 1
      Dev_log[i] = GNM1_Renshaw_Haberman$deviance
      iter_log[i] = GNM1_Renshaw_Haberman$iter
    }

  }else
  {
    Dev_log[i] = "NA"
    iter_log[i] = "NULL"
  }
  print(Dev_log)
  print(iter_log)

}
for(i in (1:10))
{
  print(i)
  GNM2_Renshaw_Haberman = gnm(Qx ~ -1 + Age_F + Mult(Age_F, Year_F) + Mult(Age_F, Cohort_F), family = binomial(link = "cloglog"), weight = InitialExpo, iterMax = 1500)

  if((is.null(GNM2_Renshaw_Haberman) == FALSE))
  {
    if((GNM2_Renshaw_Haberman$converged == TRUE) )
    {
      conv_count_cloglog = conv_count_cloglog + 1
      Dev_cloglog[i] = GNM2_Renshaw_Haberman$deviance
      iter_cloglog[i] = GNM2_Renshaw_Haberman$iter
    }
    
  }else
  {
    Dev_cloglog[i] = "NA"
    iter_cloglog[i] = "NULL"
  }
  print(Dev_cloglog)
  print(iter_cloglog)

  
}
for(i in (1:10))
{
  print(i)
  GNM3_Renshaw_Haberman = gnm(Qx ~ -1 + Age_F + Mult(Age_F, Year_F) + Mult(Age_F, Cohort_F), family = binomial(link = "logit"), weight = InitialExpo, iterMax = 1500)
  
 
  if((is.null(GNM3_Renshaw_Haberman) == FALSE))
  {
    if((GNM3_Renshaw_Haberman$converged == TRUE) )
    {
      conv_count_logit = conv_count_logit + 1
      Dev_logit[i] = GNM3_Renshaw_Haberman$deviance
      iter_logit[i] = GNM3_Renshaw_Haberman$iter
    }
  }else
  {
    Dev_logit[i] = "NA"
    iter_logit[i] = "NULL"
  }
  print(Dev_logit)
  print(iter_logit)
  
}


set.seed(10)
GNM1_Renshaw_Haberman = gnm(dx ~ -1 + Age_F + Mult(Age_F, Year_F) + Mult(Age_F, Cohort_F) + offset(log(CentralExpo)), family = poisson(link = "log"), iterMax = 1500)
set.seed(10)
GNM2_Renshaw_Haberman = gnm(Qx ~ -1 + Age_F + Mult(Age_F, Year_F) + Mult(Age_F, Cohort_F), family = binomial(link = "cloglog"), weight = InitialExpo, iterMax = 1500)
set.seed(10)
GNM3_Renshaw_Haberman = gnm(Qx ~ -1 + Age_F + Mult(Age_F, Year_F) + Mult(Age_F, Cohort_F), family = binomial(link = "logit"), weight = InitialExpo, iterMax = 1500)

Renshaw_Haberman = c(GNM1_Renshaw_Haberman$deviance, GNM2_Renshaw_Haberman$deviance, GNM3_Renshaw_Haberman$deviance)
Renshaw_Haberman

Renshaw_Haberman_dev_Poisson = c(GNM1_Renshaw_Haberman$deviance,
                           Deviance.FromBinomial(dx_obs = dx, qx_fit = GNM2_Renshaw_Haberman$fitted.values, initial_ETR = InitialExpo, devType = "Poisson"),
                           Deviance.FromBinomial(dx_obs = dx, qx_fit = GNM3_Renshaw_Haberman$fitted.values, initial_ETR = InitialExpo, devType = "Poisson"))

GNM1_Renshaw_Haberman$iter
GNM2_Renshaw_Haberman$iter
GNM3_Renshaw_Haberman$iter

## with the Poisson assumption, we find
   # log link
Beta_X_1_RH_log_R = GNM1_Renshaw_Haberman$coefficients[1:n_a]
Beta_X_2_RH_log_R = GNM1_Renshaw_Haberman$coefficients[(n_a +1):(2*n_a)]
Kappa_t_1_RH_log_R = GNM1_Renshaw_Haberman$coefficients[(2*n_a +1):((2*n_a) +n_y)]
Beta_X_3_RH_log_R = GNM1_Renshaw_Haberman$coefficients[((2*n_a) +n_y+1):((2*n_a) +n_y +n_a)]
Gamma_t_x_RH_log_R = GNM1_Renshaw_Haberman$coefficients[((2*n_a) +n_y +n_a +1):(length(GNM1_Renshaw_Haberman$coefficients))]

###Calcul des paramètres sur base des contraintes choisies à partir de ceux estimé par R

Beta_X_2_RH_log_mean = sum(Beta_X_2_RH_log_R)/n_a
Beta_X_3_RH_log_mean = sum(Beta_X_3_RH_log_R)/n_a
Kappa_t_1_RH_log_mean = sum(Kappa_t_1_RH_log_R)/n_y
Gamma_t_x_RH_log_mean = sum(Gamma_t_x_RH_log_R)/n_c

#Then
Beta_X_2_RH_log = Beta_X_2_RH_log_R/(n_a * Beta_X_2_RH_log_mean)
Kappa_t_1_RH_log = n_a*Beta_X_2_RH_log_mean * (Kappa_t_1_RH_log_R - Kappa_t_1_RH_log_mean * (rep(1,n_y)))
Beta_X_3_RH_log = Beta_X_3_RH_log_R/(n_a * Beta_X_3_RH_log_mean)
Gamma_t_x_RH_log = (n_a * Beta_X_3_RH_log_mean) * (Gamma_t_x_RH_log_R - Gamma_t_x_RH_log_mean * rep(1,n_c))
Beta_X_1_RH_log = Beta_X_1_RH_log_R + Kappa_t_1_RH_log_mean * Beta_X_2_RH_log_R + Gamma_t_x_RH_log_mean * Beta_X_3_RH_log_R


eta_RH_1 = matrix( nrow = n_a, ncol = n_y, byrow = FALSE)
for(j in 1:n_y)
{
  for(i in 1:n_a)
  {
    eta_RH_1[i,j] = Beta_X_1_RH_log[i] + Beta_X_2_RH_log[i] * Kappa_t_1_RH_log[j] + Beta_X_3_RH_log[i] * Gamma_t_x_RH_log[n_a - (i-1) + (j-1)] 
  }
}
eta_RH_1 = matrix(eta_RH_1, nrow = n_a, ncol = n_y, byrow = FALSE)
mu_log_RH_1 = exp(eta_RH_1)
persp(Age,Year, mu_log_RH_1, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "RH_log")


## with the Binomial assumption, we find
   # cloglog link
Beta_X_1_RH_cloglog_R = GNM2_Renshaw_Haberman$coefficients[1:n_a]
Beta_X_2_RH_cloglog_R = GNM2_Renshaw_Haberman$coefficients[(n_a +1):(2*n_a)]
Kappa_t_1_RH_cloglog_R = GNM2_Renshaw_Haberman$coefficients[(2*n_a +1):((2*n_a) +n_y)]
Beta_X_3_RH_cloglog_R = GNM2_Renshaw_Haberman$coefficients[((2*n_a) +n_y+1):((2*n_a) +n_y +n_a)]
Gamma_t_x_RH_cloglog_R = GNM2_Renshaw_Haberman$coefficients[((2*n_a) +n_y +n_a +1):(length(GNM1_Renshaw_Haberman$coefficients))]

###Calcul des paramètres sur base des contraintes choisies à partir de ceux estimé par R

Beta_X_2_RH_cloglog_mean = sum(Beta_X_2_RH_cloglog_R)/n_a
Beta_X_3_RH_cloglog_mean = sum(Beta_X_3_RH_cloglog_R)/n_a
Kappa_t_1_RH_cloglog_mean = sum(Kappa_t_1_RH_cloglog_R)/n_y
Gamma_t_x_RH_cloglog_mean = sum(Gamma_t_x_RH_cloglog_R)/n_c

#Then
Beta_X_2_RH_cloglog = Beta_X_2_RH_cloglog_R/(n_a * Beta_X_2_RH_cloglog_mean)
Kappa_t_1_RH_cloglog = n_a*Beta_X_2_RH_cloglog_mean * (Kappa_t_1_RH_cloglog_R - Kappa_t_1_RH_cloglog_mean * (rep(1,n_y)))
Beta_X_3_RH_cloglog = Beta_X_3_RH_cloglog_R/(n_a * Beta_X_3_RH_cloglog_mean)
Gamma_t_x_RH_cloglog = (n_a * Beta_X_3_RH_cloglog_mean) * (Gamma_t_x_RH_cloglog_R - Gamma_t_x_RH_cloglog_mean * rep(1,n_c))
Beta_X_1_RH_cloglog = Beta_X_1_RH_cloglog_R + Kappa_t_1_RH_cloglog_mean * Beta_X_2_RH_cloglog_R + Gamma_t_x_RH_cloglog_mean * Beta_X_3_RH_cloglog_R


eta_RH_2 = matrix( nrow = n_a, ncol = n_y, byrow = FALSE)
for(j in 1:n_y)
{
  for(i in 1:n_a)
  {
    eta_RH_2[i,j] = Beta_X_1_RH_cloglog[i] + Beta_X_2_RH_cloglog[i] * Kappa_t_1_RH_cloglog[j] + Beta_X_3_RH_cloglog[i] * Gamma_t_x_RH_cloglog[n_a - (i-1) + (j-1)] 
  }
}
qx_cloglog_RH_2 = 1-exp(-exp(eta_RH_2))
mu_cloglog_RH_2 = -log(1-qx_cloglog_RH_2)
persp(Age,Year, mu_cloglog_RH_2, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "RH_cloglog")



    # logit link
Beta_X_1_RH_logit_R = GNM3_Renshaw_Haberman$coefficients[1:n_a]
Beta_X_2_RH_logit_R = GNM3_Renshaw_Haberman$coefficients[(n_a +1):(2*n_a)]
Kappa_t_1_RH_logit_R = GNM3_Renshaw_Haberman$coefficients[(2*n_a +1):((2*n_a) +n_y)]
Beta_X_3_RH_logit_R = GNM3_Renshaw_Haberman$coefficients[((2*n_a) +n_y+1):((2*n_a) +n_y +n_a)]
Gamma_t_x_RH_logit_R = GNM3_Renshaw_Haberman$coefficients[((2*n_a) +n_y +n_a +1):(length(GNM1_Renshaw_Haberman$coefficients))]

###Calcul des paramètres sur base des contraintes choisies à partir de ceux estimé par R

Beta_X_2_RH_logit_mean = sum(Beta_X_2_RH_logit_R)/n_a
Beta_X_3_RH_logit_mean = sum(Beta_X_3_RH_logit_R)/n_a
Kappa_t_1_RH_logit_mean = sum(Kappa_t_1_RH_logit_R)/n_y
Gamma_t_x_RH_logit_mean = sum(Gamma_t_x_RH_logit_R)/n_c

#Then
Beta_X_2_RH_logit = Beta_X_2_RH_logit_R/(n_a * Beta_X_2_RH_logit_mean)
Kappa_t_1_RH_logit = n_a*Beta_X_2_RH_logit_mean * (Kappa_t_1_RH_logit_R - Kappa_t_1_RH_logit_mean * (rep(1,n_y)))
Beta_X_3_RH_logit = Beta_X_3_RH_logit_R/(n_a * Beta_X_3_RH_logit_mean)
Gamma_t_x_RH_logit = (n_a * Beta_X_3_RH_logit_mean) * (Gamma_t_x_RH_logit_R - Gamma_t_x_RH_logit_mean * rep(1,n_c))
Beta_X_1_RH_logit = Beta_X_1_RH_logit_R + Kappa_t_1_RH_logit_mean * Beta_X_2_RH_logit_R + Gamma_t_x_RH_logit_mean * Beta_X_3_RH_logit_R


eta_RH_3 = matrix( nrow = n_a, ncol = n_y, byrow = FALSE)
for(j in 1:n_y)
{
  for(i in 1:n_a)
  {
    eta_RH_3[i,j] = Beta_X_1_RH_logit[i] + Beta_X_2_RH_logit[i] * Kappa_t_1_RH_logit[j] + Beta_X_3_RH_logit[i] * Gamma_t_x_RH_logit[n_a - (i-1) + (j-1)] 
  }
}
qx_logit_RH_3 = (exp(eta_RH_3))/(1 + exp(eta_RH_3))
mu_logit_RH_3 = -log(1-qx_logit_RH_3)
persp(Age,Year, mu_logit_RH_3, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "RH_logit")



Resids_Plot(residuals_mu = matrix(c(rstandard(GNM1_Renshaw_Haberman)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "Renshaw_Haberman_poisson_log" )
Resids_Plot(residuals_mu = matrix(c(rstandard(GNM2_Renshaw_Haberman)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "Renshaw_Haberman_Binomial_cloglog" )
Resids_Plot(residuals_mu = matrix(c(rstandard(GNM3_Renshaw_Haberman)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "Renshaw_Haberman_Binomial_logit" )



Plot_Superpose_Coef(TheRange = c(1:n_a), GNM1_Renshaw_Haberman, GNM2_Renshaw_Haberman, GNM3_Renshaw_Haberman, GLM1_name = "Renshaw_Haberman_Poisson_log" , GLM2_name = "Renshaw_Haberman_Binomial_cloglog" , GLM3_name = "Renshaw_Haberman_Binomial_logit" , main_name = "Renshaw_Haberman", y_name = "Age_coef_Beta_1", pos_legend = "topleft")

Plot_Superpose_Coef(TheRange = c((n_a +1):( 2* n_a)), GNM1_Renshaw_Haberman, GNM2_Renshaw_Haberman, GNM3_Renshaw_Haberman, GLM1_name = "Renshaw_Haberman_Poisson_log" , GLM2_name = "Renshaw_Haberman_Binomial_cloglog" , GLM3_name = "Renshaw_Haberman_Binomial_logit" , main_name = "Renshaw_Haberman", y_name = "Age_coef_Beta_2", pos_legend = 'topright')

Plot_Superpose_Coef(TheRange = c((2* n_a + 1): (2* n_a + n_y)), GNM1_Renshaw_Haberman, GNM2_Renshaw_Haberman, GNM3_Renshaw_Haberman, GLM1_name = "Renshaw_Haberman_Poisson_log" , GLM2_name = "Renshaw_Haberman_Binomial_cloglog" , GLM3_name = "Renshaw_Haberman_Binomial_logit" , main_name = "Renshaw_Haberman", y_name = "Year_coef_Kappa_2", pos_legend = "topright")

Plot_Superpose_Coef(TheRange = c((2*n_a + n_y + 1): (2*n_a + n_y + n_a)), GNM1_Renshaw_Haberman, GNM2_Renshaw_Haberman, GNM3_Renshaw_Haberman, GLM1_name = "Renshaw_Haberman_Poisson_log" , GLM2_name = "Renshaw_Haberman_Binomial_cloglog" , GLM3_name = "Renshaw_Haberman_Binomial_logit" , main_name = "Renshaw_Haberman", y_name = "Age_coef_Beta_3", pos_legend = "bottomleft")

Plot_Superpose_Coef(TheRange = c((3*n_a + n_y + 1): (3*n_a + n_y + n_c)), GNM1_Renshaw_Haberman, GNM2_Renshaw_Haberman, GNM3_Renshaw_Haberman, GLM1_name = "Renshaw_Haberman_Poisson_log" , GLM2_name = "Renshaw_Haberman_Binomial_cloglog" , GLM3_name = "Renshaw_Haberman_Binomial_logit" , main_name = "Renshaw_Haberman", y_name = "Cohort_coef_Gamma_3", pos_legend = "topleft")


# Comparison of the kappa coefficients
plot(Kappa_t_1_RH_cloglog, type = "l")
legend("bottomright", "sum over kappa = -3.2e-15")
plot(Kappa_t_1_RH_cloglog_R, type = "l")
legend("topright", paste("sum over kappa = ", round(sum(Kappa_t_1_RH_cloglog_R), 5)))

    #=>        ## Marche pas alors que censé être plus optimal, je dois tester le fit en fixant simplement le facteur age
               ## Même si du coup c'est déjà pas mal que j'obtienne une réponse 


        #####################################################
        #####     Recupatilatif des Modèles via gnm     #####   Done
        #####################################################


tab_GNM = data.frame(Lee_Carter, Renshaw_Haberman)
rownames(tab_GNM) = c("GNM1_log", "GNM2_cloglog", "GNM3_logit")
tab_GNM







#####################################################
#### Modèle de Renshaw-Haberman avec start values ####    Done    
#####################################################


### Going to go trough the start input of the gnm function taking the estimated parameters of the APC model as start parameters
#Marche pas alors que censé être plus optimal, je dois tester le fit en fixant simplement le facteur age
a = c()
b = c()
c = c()
for( i in 1:10)
{  GNM1_Renshaw_Haberman_start = gnm(dx ~ -1 + Age_F + Mult(Age_F, Year_F) + Mult(Age_F, Cohort_F) + offset(log(CentralExpo)), 
                                     family = poisson(link = "log") ,start = Start)
GNM2_Renshaw_Haberman_start = gnm(Qx ~ -1 + Age_F + Mult(Age_F, Year_F) + Mult(Age_F, Cohort_F), 
                                  family = binomial(link = "cloglog"), weight = InitialExpo, start = Start)
GNM3_Renshaw_Haberman_start = gnm(Qx ~ -1 + Age_F + Mult(Age_F, Year_F) + Mult(Age_F, Cohort_F),
                                  family = binomial(link = "logit"), weight = InitialExpo, start = Start)

Renshaw_Haberman_start = c(GNM1_Renshaw_Haberman_start$deviance, GNM2_Renshaw_Haberman_start$deviance, GNM3_Renshaw_Haberman_start$deviance)
Renshaw_Haberman_start

a = c(a,GNM1_Renshaw_Haberman_start$iter)
b = c(b,GNM2_Renshaw_Haberman_start$iter)
c = c(c,GNM3_Renshaw_Haberman_start$iter)
}



#### For Poisson log

Alpha = GLM1_APC$coefficients[1:n_a]
Kappa = GLM1_APC$coefficients[(n_a+1) : (n_a + n_y) ]
Gamma = GLM1_APC$coefficients[(n_a + n_y +1): dim(X_APC)[2]]

Beta2 = rep(1, n_a)
Beta3 = Beta2

Beta2 = Beta2/n_a
Beta3 = Beta3/n_a

# only consider kappa and beta3 and alpha 
Start = c(NA* Alpha, Beta2, Kappa, Beta3,  NA*Gamma)

GNM1_Renshaw_Haberman_start = gnm(dx ~ -1 + Age_F + Mult(Age_F, Year_F) + Mult(Age_F, Cohort_F) + offset(log(CentralExpo)), 
                                  family = poisson(link = "log"), start = Start, iterMax = 1500) #8 est optimal en nombre d'it
GNM1_Renshaw_Haberman_start$deviance

#### For Binomial cloglog

Alpha = GLM2_APC$coefficients[1:n_a]
Kappa = GLM2_APC$coefficients[(n_a+1) : (n_a + n_y) ]
Gamma = GLM2_APC$coefficients[(n_a + n_y +1): dim(X_APC)[2]]

Beta2 = rep(1, n_a)
Beta3 = Beta2

Beta2 = Beta2/n_a
Beta3 = Beta3/n_a

# only consider kappa and beta3 and alpha 
Start = c( NA*Alpha, Beta2, Kappa, Beta3,  NA*Gamma)

GNM2_Renshaw_Haberman_start = gnm(Qx ~ -1 + Age_F + Mult(Age_F, Year_F) + Mult(Age_F, Cohort_F), 
                                  family = binomial(link = "cloglog"), weight = InitialExpo, start = Start, iterMax = 1500)



#### For Binomial logit

Alpha = GLM3_APC$coefficients[1:n_a]
Kappa = GLM3_APC$coefficients[(n_a+1) : (n_a + n_y) ]
Gamma = GLM3_APC$coefficients[(n_a + n_y +1): dim(X_APC)[2]]

Beta2 = rep(1, n_a)
Beta3 = Beta2

Beta2 = Beta2/n_a
Beta3 = Beta3/n_a

# only consider kappa and beta2/3
Start = c( NA*Alpha,Beta2, Kappa, Beta3,  NA*Gamma)


GNM3_Renshaw_Haberman_start = gnm(Qx ~ -1 + Age_F + Mult(Age_F, Year_F) + Mult(Age_F, Cohort_F),
                                  family = binomial(link = "logit"), weight = InitialExpo, start = Start, iterMax = 1500)




Renshaw_Haberman_start_dev_Poisson = c(GNM1_Renshaw_Haberman_start$deviance,
                           Deviance.FromBinomial(dx_obs = dx, qx_fit = GNM2_Renshaw_Haberman_start$fitted.values, initial_ETR = InitialExpo, devType = "Poisson"),
                           Deviance.FromBinomial(dx_obs = dx, qx_fit = GNM3_Renshaw_Haberman_start$fitted.values, initial_ETR = InitialExpo, devType = "Poisson"))



## with the Poisson assumption, we find
# log link
Beta_X_1_RH_log_R = GNM1_Renshaw_Haberman_start$coefficients[1:n_a]
Beta_X_2_RH_log_R = GNM1_Renshaw_Haberman_start$coefficients[(n_a +1):(2*n_a)]
Kappa_t_1_RH_log_R = GNM1_Renshaw_Haberman_start$coefficients[(2*n_a +1):((2*n_a) +n_y)]
Beta_X_3_RH_log_R = GNM1_Renshaw_Haberman_start$coefficients[((2*n_a) +n_y+1):((2*n_a) +n_y +n_a)]
Gamma_t_x_RH_log_R = GNM1_Renshaw_Haberman_start$coefficients[((2*n_a) +n_y +n_a +1):(length(GNM1_Renshaw_Haberman$coefficients))]

###Calcul des paramètres sur base des contraintes choisies à partir de ceux estimé par R

Beta_X_2_RH_log_mean = sum(Beta_X_2_RH_log_R)/n_a
Beta_X_3_RH_log_mean = sum(Beta_X_3_RH_log_R)/n_a
Kappa_t_1_RH_log_mean = sum(Kappa_t_1_RH_log_R)/n_y
Gamma_t_x_RH_log_mean = sum(Gamma_t_x_RH_log_R)/n_c

#Then
Beta_X_2_RH_log = Beta_X_2_RH_log_R/(n_a * Beta_X_2_RH_log_mean)
Kappa_t_1_RH_log = n_a*Beta_X_2_RH_log_mean * (Kappa_t_1_RH_log_R - Kappa_t_1_RH_log_mean * (rep(1,n_y)))
Beta_X_3_RH_log = Beta_X_3_RH_log_R/(n_a * Beta_X_3_RH_log_mean)
Gamma_t_x_RH_log = (n_a * Beta_X_3_RH_log_mean) * (Gamma_t_x_RH_log_R - Gamma_t_x_RH_log_mean * rep(1,n_c))
Beta_X_1_RH_log = Beta_X_1_RH_log_R + Kappa_t_1_RH_log_mean * Beta_X_2_RH_log_R + Gamma_t_x_RH_log_mean * Beta_X_3_RH_log_R


eta_RH_1 = matrix( nrow = n_a, ncol = n_y, byrow = FALSE)
for(j in 1:n_y)
{
  for(i in 1:n_a)
  {
    eta_RH_1[i,j] = Beta_X_1_RH_log[i] + Beta_X_2_RH_log[i] * Kappa_t_1_RH_log[j] + Beta_X_3_RH_log[i] * Gamma_t_x_RH_log[n_a - (i-1) + (j-1)] 
  }
}
eta_RH_1 = matrix(eta_RH_1, nrow = n_a, ncol = n_y, byrow = FALSE)
mu_log_RH_1 = exp(eta_RH_1)
persp(Age,Year, mu_log_RH_1, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "RH_log")


## with the Binomial assumption, we find
# cloglog link
Beta_X_1_RH_cloglog_R = GNM2_Renshaw_Haberman_start$coefficients[1:n_a]
Beta_X_2_RH_cloglog_R = GNM2_Renshaw_Haberman_start$coefficients[(n_a +1):(2*n_a)]
Kappa_t_1_RH_cloglog_R = GNM2_Renshaw_Haberman_start$coefficients[(2*n_a +1):((2*n_a) +n_y)]
Beta_X_3_RH_cloglog_R = GNM2_Renshaw_Haberman_start$coefficients[((2*n_a) +n_y+1):((2*n_a) +n_y +n_a)]
Gamma_t_x_RH_cloglog_R = GNM2_Renshaw_Haberman_start$coefficients[((2*n_a) +n_y +n_a +1):(length(GNM1_Renshaw_Haberman$coefficients))]

###Calcul des paramètres sur base des contraintes choisies à partir de ceux estimé par R

Beta_X_2_RH_cloglog_mean = sum(Beta_X_2_RH_cloglog_R)/n_a
Beta_X_3_RH_cloglog_mean = sum(Beta_X_3_RH_cloglog_R)/n_a
Kappa_t_1_RH_cloglog_mean = sum(Kappa_t_1_RH_cloglog_R)/n_y
Gamma_t_x_RH_cloglog_mean = sum(Gamma_t_x_RH_cloglog_R)/n_c

#Then
Beta_X_2_RH_cloglog = Beta_X_2_RH_cloglog_R/(n_a * Beta_X_2_RH_cloglog_mean)
Kappa_t_1_RH_cloglog = n_a*Beta_X_2_RH_cloglog_mean * (Kappa_t_1_RH_cloglog_R - Kappa_t_1_RH_cloglog_mean * (rep(1,n_y)))
Beta_X_3_RH_cloglog = Beta_X_3_RH_cloglog_R/(n_a * Beta_X_3_RH_cloglog_mean)
Gamma_t_x_RH_cloglog = (n_a * Beta_X_3_RH_cloglog_mean) * (Gamma_t_x_RH_cloglog_R - Gamma_t_x_RH_cloglog_mean * rep(1,n_c))
Beta_X_1_RH_cloglog = Beta_X_1_RH_cloglog_R + Kappa_t_1_RH_cloglog_mean * Beta_X_2_RH_cloglog_R + Gamma_t_x_RH_cloglog_mean * Beta_X_3_RH_cloglog_R


eta_RH_2 = matrix( nrow = n_a, ncol = n_y, byrow = FALSE)
for(j in 1:n_y)
{
  for(i in 1:n_a)
  {
    eta_RH_2[i,j] = Beta_X_1_RH_cloglog[i] + Beta_X_2_RH_cloglog[i] * Kappa_t_1_RH_cloglog[j] + Beta_X_3_RH_cloglog[i] * Gamma_t_x_RH_cloglog[n_a - (i-1) + (j-1)] 
  }
}
qx_cloglog_RH_2 = 1-exp(-exp(eta_RH_2))
mu_cloglog_RH_2 = -log(1-qx_cloglog_RH_2)
persp(Age,Year, mu_cloglog_RH_2, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "RH_cloglog")



# logit link
Beta_X_1_RH_logit_R = GNM3_Renshaw_Haberman_start$coefficients[1:n_a]
Beta_X_2_RH_logit_R = GNM3_Renshaw_Haberman_start$coefficients[(n_a +1):(2*n_a)]
Kappa_t_1_RH_logit_R = GNM3_Renshaw_Haberman_start$coefficients[(2*n_a +1):((2*n_a) +n_y)]
Beta_X_3_RH_logit_R = GNM3_Renshaw_Haberman_start$coefficients[((2*n_a) +n_y+1):((2*n_a) +n_y +n_a)]
Gamma_t_x_RH_logit_R = GNM3_Renshaw_Haberman_start$coefficients[((2*n_a) +n_y +n_a +1):(length(GNM1_Renshaw_Haberman$coefficients))]

###Calcul des paramètres sur base des contraintes choisies à partir de ceux estimé par R

Beta_X_2_RH_logit_mean = sum(Beta_X_2_RH_logit_R)/n_a
Beta_X_3_RH_logit_mean = sum(Beta_X_3_RH_logit_R)/n_a
Kappa_t_1_RH_logit_mean = sum(Kappa_t_1_RH_logit_R)/n_y
Gamma_t_x_RH_logit_mean = sum(Gamma_t_x_RH_logit_R)/n_c

#Then
Beta_X_2_RH_logit = Beta_X_2_RH_logit_R/(n_a * Beta_X_2_RH_logit_mean)
Kappa_t_1_RH_logit = n_a*Beta_X_2_RH_logit_mean * (Kappa_t_1_RH_logit_R - Kappa_t_1_RH_logit_mean * (rep(1,n_y)))
Beta_X_3_RH_logit = Beta_X_3_RH_logit_R/(n_a * Beta_X_3_RH_logit_mean)
Gamma_t_x_RH_logit = (n_a * Beta_X_3_RH_logit_mean) * (Gamma_t_x_RH_logit_R - Gamma_t_x_RH_logit_mean * rep(1,n_c))
Beta_X_1_RH_logit = Beta_X_1_RH_logit_R + Kappa_t_1_RH_logit_mean * Beta_X_2_RH_logit_R + Gamma_t_x_RH_logit_mean * Beta_X_3_RH_logit_R


eta_RH_3 = matrix( nrow = n_a, ncol = n_y, byrow = FALSE)
for(j in 1:n_y)
{
  for(i in 1:n_a)
  {
    eta_RH_3[i,j] = Beta_X_1_RH_logit[i] + Beta_X_2_RH_logit[i] * Kappa_t_1_RH_logit[j] + Beta_X_3_RH_logit[i] * Gamma_t_x_RH_logit[n_a - (i-1) + (j-1)] 
  }
}
qx_logit_RH_3 = (exp(eta_RH_3))/(1 + exp(eta_RH_3))
mu_logit_RH_3 = -log(1-qx_logit_RH_3)
persp(Age,Year, mu_logit_RH_3, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "RH_logit")

















#################################################################################
######                 Ecriture du modèle M4 2D-Splines                   #######
#################################################################################

Age = c(50:90)
Year = c(1960:2015)
n_a = length(Age)
n_y = length(Year)
## Création de la matrice de décès
dx_Mat = matrix(dx, nrow = n_a, ncol = n_y, byrow = FALSE )
Qx_Mat = matrix(Qx, nrow = n_a, ncol = n_y, byrow = FALSE )

## Création de la matrice d'exposition
InitialExposure_Mat = matrix(InitialExpo, nrow = n_a, ncol = n_y, byrow = FALSE)
CentralExposure_Mat = matrix(CentralExpo, nrow = n_a, ncol = n_y, byrow = FALSE)

## Utilisation de la fonction de Mort2DSmooth: après réécriture
P_Spline_2D_Poisson_log = Mort2Dsmooth_poisson(x = Age, y = Year, Z = dx_Mat, offset = log(CentralExposure_Mat), overdispersion = FALSE)
P_Spline_2D_Poisson_log$aic
P_Spline_2D_Poisson_log_over = Mort2Dsmooth_poisson(x = Age, y = Year, Z = dx_Mat, offset = log(CentralExposure_Mat), overdispersion = TRUE)
P_Spline_2D_Poisson_log_over$aic


P_Spline_2D_Binomial_cloglog = Mort2Dsmooth_Binomial_cloglog(x = Age, y = Year, Z = Qx_Mat, W = InitialExposure_Mat, overdispersion = FALSE)
P_Spline_2D_Binomial_cloglog$aic   
P_Spline_2D_Binomial_cloglog_over = Mort2Dsmooth_Binomial_cloglog(x = Age, y = Year, Z = Qx_Mat, W = InitialExposure_Mat, overdispersion = TRUE)
P_Spline_2D_Binomial_cloglog_over$aic  


P_Spline_2D_Binomial_logit = Mort2Dsmooth_Binomial_logit(x = Age, y = Year, Z = Qx_Mat, W = InitialExposure_Mat, overdispersion = FALSE)
P_Spline_2D_Binomial_logit$aic  
P_Spline_2D_Binomial_logit_over = Mort2Dsmooth_Binomial_logit(x = Age, y = Year, Z = Qx_Mat, W = InitialExposure_Mat, overdispersion = TRUE)
P_Spline_2D_Binomial_logit_over$aic 



P_Spline_2D = c(P_Spline_2D_Poisson_log$dev ,P_Spline_2D_Binomial_cloglog$dev,P_Spline_2D_Binomial_logit$dev)
P_Spline_2D_aic = c(P_Spline_2D_Poisson_log$aic ,P_Spline_2D_Binomial_cloglog$aic,P_Spline_2D_Binomial_logit$aic)
P_Spline_2D_psi2 = c(P_Spline_2D_Poisson_log$psi2 ,P_Spline_2D_Binomial_cloglog$psi2,P_Spline_2D_Binomial_logit$psi2)
## then with Poisson assumption, we find :
##For the log link
#eta_2D_B_Spline_1 = (GAM1_2D_B_Spline$model %*% GAM1_2D_B_Spline$coefficients)
mu_log_2D_P_Spline_1 = exp(P_Spline_2D_Poisson_log$logmortality)
mu_log_2D_P_Spline_1 = matrix(mu_log_2D_P_Spline_1, nrow = n_a, ncol = n_y, byrow = FALSE)  #Mat format
persp(Age,Year, mu_log_2D_P_Spline_1, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "2D_P_Spline_log")

## then with Binomial assumption, we find :
##For the cloglog link
mu_cloglog_2D_P_Spline_2 = -log(1-P_Spline_2D_Binomial_cloglog$fitted.values)
persp(Age,Year, mu_cloglog_2D_P_Spline_2, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "2D_P_Spline_cloglog")

##For the logit link
mu_logit_2D_P_Spline_3 = -log(1-P_Spline_2D_Binomial_logit$fitted.values)
persp(Age,Year, mu_logit_2D_P_Spline_3, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "2D_P_Spline_logit")





Is.it.too.smooth = function(Qx_smooth_mat, Qx_notsmooth_mat, Exposure_mat)
{
  p = 2 # Because we smooth based on 2 paramaters to chose the number of internal knots for each axes ages an years
  m = n_a 
  for(i in 1:n_y){
    StatsTestH = sum( (Exposure_mat[,i]) * ((c(Qx_smooth_mat[,i]) - Qx_notsmooth_mat[,i])^2)/(c(Qx_smooth_mat[,i]) * (i-c(Qx_smooth_mat[,i]) )))
    
    Chisq_995 = qchisq(0.99, m - p - 1)

    if(StatsTestH > Chisq_995){
      print(paste("The Year ", Year[i], " is oversmoothed"))
    }
    
  }
}
Is.it.too.smooth(Qx_smooth_mat = P_Spline_2D_Poisson_log$fitted.values/LT_M_1960_50_90$lx, Qx_notsmooth_mat = Qx_Mat, Exposure_mat = InitialExposure_Mat)
Is.it.too.smooth(Qx_smooth_mat = P_Spline_2D_Binomial_logit$fitted.values, Qx_notsmooth_mat = Qx_Mat, Exposure_mat = InitialExposure_Mat)
Is.it.too.smooth(Qx_smooth_mat = P_Spline_2D_Binomial_cloglog$fitted.values, Qx_notsmooth_mat = Qx_Mat, Exposure_mat = InitialExposure_Mat)



resid_P_Spline_log = (Compute_resids(fit_values = P_Spline_2D_Poisson_log$fitted.values, FamillyType = "Poisson") - mean(Compute_resids(fit_values = P_Spline_2D_Poisson_log$fitted.values, FamillyType = "Poisson"))) / sd(Compute_resids(fit_values = P_Spline_2D_Poisson_log$fitted.values, FamillyType = "Poisson"))
resid_P_Spline_cloglog = (Compute_resids(fit_values = P_Spline_2D_Binomial_cloglog$fitted.values, FamillyType = "Binomial") - mean(Compute_resids(fit_values = P_Spline_2D_Binomial_cloglog$fitted.values, FamillyType = "Binomial"))) / sd(Compute_resids(fit_values = P_Spline_2D_Binomial_cloglog$fitted.values, FamillyType = "Binomial"))
resid_P_Spline_logit = (Compute_resids(fit_values = P_Spline_2D_Binomial_logit$fitted.values, FamillyType = "Binomial") - mean(Compute_resids(fit_values = P_Spline_2D_Binomial_logit$fitted.values, FamillyType = "Binomial"))) / sd(Compute_resids(fit_values = P_Spline_2D_Binomial_logit$fitted.values, FamillyType = "Binomial"))


Resids_Plot(residuals_mu = matrix(c(resid_P_Spline_log), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "2D_P_Spline_poisson_log" )
Resids_Plot(residuals_mu = matrix(c(resid_P_Spline_cloglog), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "2D_P_Spline_Binomial_cloglog" )
Resids_Plot(residuals_mu = matrix(c(resid_P_Spline_logit), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "2D_P_Spline_Binomial_logit" )



# Plots with overdispersion
## then with Poisson assumption, we find :
##For the log link
#eta_2D_B_Spline_1 = (GAM1_2D_B_Spline$model %*% GAM1_2D_B_Spline$coefficients)
mu_log_over_2D_P_Spline_1 = exp(P_Spline_2D_Poisson_log_over$logmortality)
mu_log_over_2D_P_Spline_1 = matrix(mu_log_over_2D_P_Spline_1, nrow = n_a, ncol = n_y, byrow = FALSE)  #Mat format
persp(Age,Year, mu_log_over_2D_P_Spline_1, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "2D_P_Spline_log_over")

## then with Binomial assumption, we find :
##For the cloglog_over link
mu_cloglog_over_2D_P_Spline_2 = -log(1-P_Spline_2D_Binomial_cloglog_over$fitted.values)
persp(Age,Year, mu_cloglog_over_2D_P_Spline_2, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "2D_P_Spline_cloglog_over")

##For the logit link
mu_logit_over_2D_P_Spline_3 = -log(1-P_Spline_2D_Binomial_logit_over$fitted.values)
persp(Age,Year, mu_logit_over_2D_P_Spline_3, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "2D_P_Spline_logit_over")


resid_P_Spline_log_over = (Compute_resids(fit_values = P_Spline_2D_Poisson_log_over$fitted.values, FamillyType = "Poisson") - mean(Compute_resids(fit_values = P_Spline_2D_Poisson_log_over$fitted.values, FamillyType = "Poisson"))) / sd(Compute_resids(fit_values = P_Spline_2D_Poisson_log_over$fitted.values, FamillyType = "Poisson"))
resid_P_Spline_cloglog_over = (Compute_resids(fit_values = P_Spline_2D_Binomial_cloglog_over$fitted.values, FamillyType = "Binomial") - mean(Compute_resids(fit_values = P_Spline_2D_Binomial_cloglog_over$fitted.values, FamillyType = "Binomial"))) / sd(Compute_resids(fit_values = P_Spline_2D_Binomial_cloglog_over$fitted.values, FamillyType = "Binomial"))
resid_P_Spline_logit_over = (Compute_resids(fit_values = P_Spline_2D_Binomial_logit_over$fitted.values, FamillyType = "Binomial") - mean(Compute_resids(fit_values = P_Spline_2D_Binomial_logit_over$fitted.values, FamillyType = "Binomial"))) / sd(Compute_resids(fit_values = P_Spline_2D_Binomial_logit_over$fitted.values, FamillyType = "Binomial"))


Resids_Plot(residuals_mu = matrix(c(resid_P_Spline_log_over), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "2D_P_Spline_poisson_over_log_over" )
Resids_Plot(residuals_mu = matrix(c(resid_P_Spline_cloglog_over), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "2D_P_Spline_Binomial_over_cloglog_over" )
Resids_Plot(residuals_mu = matrix(c(resid_P_Spline_logit_over), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "2D_P_Spline_Binomial_over_logit_over" )




########################################################################################################
######                 Ecriture des modèles glm AP avec un lissage B-Spline fait à priori        #######   
########################################################################################################


B_x = MortSmooth_bbase(unique(LT_M_1960_50_90$Age), min(LT_M_1960_50_90$Age), max(LT_M_1960_50_90$Age) , ndx = 11, deg = 3)
B_y = MortSmooth_bbase(unique(LT_M_1960_50_90$Year), min(LT_M_1960_50_90$Year), max(LT_M_1960_50_90$Year) , ndx = 14, deg = 3)


Stock = c()
X_model = c()
for( i in 1:n_y)
{
  Stock = c()
  for(j in 1:dim(B_y)[2])
  {
    Stock = cbind(Stock,B_y[i,j] * B_x)
  }
  X_model = rbind(X_model,Stock)
}

##### Results through only the cubic B-spline regression matrix function 

GLM1_2D_B_Spline = glm(dx ~ -1 + X_model + offset(log(CentralExpo)), family = poisson(link = "log"))
GLM2_2D_B_Spline = glm(Qx ~ -1 + X_model, family = binomial(link = "cloglog"), weight = InitialExpo)
GLM3_2D_B_Spline = glm(Qx ~ -1 + X_model, family = binomial(link = "logit"), weight = InitialExpo)

GLM_Splines = c(GLM1_2D_B_Spline$deviance, GLM2_2D_B_Spline$deviance, GLM3_2D_B_Spline$deviance)
GLM_Splines
GLM_Splines_aic = c(GLM1_2D_B_Spline$aic, GLM2_2D_B_Spline$aic, GLM3_2D_B_Spline$aic)
GLM_Splines_aic

## then with Poisson assumption, we find :
##For the log link
#eta_2D_B_Spline_1 = (GAM1_2D_B_Spline$model %*% GAM1_2D_B_Spline$coefficients)
eta_2D_B_Spline_1 = (X_model %*% GLM1_2D_B_Spline$coefficients)
mu_log_2D_B_Spline_1 = exp(eta_2D_B_Spline_1)
mu_log_2D_B_Spline_1 = matrix(mu_log_2D_B_Spline_1, nrow = n_a, ncol = n_y, byrow = FALSE)  #Mat format
persp(Age,Year, mu_log_2D_B_Spline_1, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "2D_B_Spline_log")

## then with Binomial assumption, we find :
##For the cloglog link
eta_2D_B_Spline_2 = (X_model %*% GLM2_2D_B_Spline$coefficients)
qx_cloglog_2D_B_Spline_2 = 1-exp(-exp(eta_2D_B_Spline_2))
qx_cloglog_2D_B_Spline_2 = matrix(qx_cloglog_2D_B_Spline_2, nrow = n_a, ncol = n_y, byrow = FALSE)  #Mat format
mu_cloglog_2D_B_Spline_2 = -log(1-qx_cloglog_2D_B_Spline_2)
persp(Age,Year, mu_cloglog_2D_B_Spline_2, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "2D_B_Spline_cloglog")

##For the logit link
eta_2D_B_Spline_3 = (X_model %*% GLM3_2D_B_Spline$coefficients)
qx_logit_2D_B_Spline_3 = (exp(eta_2D_B_Spline_3))/(1 + exp(eta_2D_B_Spline_3))
qx_logit_2D_B_Spline_3 = matrix(qx_logit_2D_B_Spline_3, nrow = n_a, ncol = n_y, byrow = FALSE)  #Mat format
mu_logit_2D_B_Spline_3 = -log(1-qx_logit_2D_B_Spline_3)
persp(Age,Year, mu_logit_2D_B_Spline_3, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "2D_B_Spline_logit")




Resids_Plot(residuals_mu = matrix(c(rstandard(GLM1_2D_B_Spline)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "2D_B_Spline_poisson_log" )
Resids_Plot(residuals_mu = matrix(c(rstandard(GLM2_2D_B_Spline)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "2D_B_Spline_Binomial_cloglog" )
Resids_Plot(residuals_mu = matrix(c(rstandard(GLM3_2D_B_Spline)), nrow = n_a, ncol = n_y, byrow = FALSE), blur = 4, Main_name = "2D_B_Spline_Binomial_logit" )



