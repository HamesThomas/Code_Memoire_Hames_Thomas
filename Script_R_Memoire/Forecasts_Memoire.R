

#######################################################################################################
######                Projection of the mortality and the Kappa_t parameter                     #######   
#######################################################################################################
## Explanation ARIMA: https://www.datascience.com/blog/introduction-to-forecasting-with-arima-in-r-learn-data-science-tutorials
install.packages("StMoMo")
library(forecast)
library(StMoMo)
library(fanplot)
library(demography)
install.packages("stats")
library(arima)
##StMoMo vignette

#function forecast forecasts the mutlivariate random walk with drift in equation 19 using the approach described in RH
#Where function Arima is used to estimate and forecast the ARIMA processes of Eq. 20 and 21

Number_Years_Forecasted = 51 


###############################################################
######          Forecast CBD cohort quad effect         #######   Done mais pas l'hypothèse des kappa-t ils sont forcés en arima
###############################################################

##For the GLM3
Kappa_1_t_CBD_cohort_quad_3 = GLM3_CBD_cohort_quad$coefficients[1:n_y]
Kappa_1_t_Forecasted_CBD_cohort_quad_3 = forecast(Arima(Kappa_1_t_CBD_cohort_quad_3, c(1,0,1), include.drift = TRUE), h = Number_Years_Forecasted)
Kappa_1_t_Forecasted_CBD_cohort_quad_3
plot(Kappa_1_t_Forecasted_CBD_cohort_quad_3, main = "CBD quad cohort logit forecast kappa_1 ARIMA(1,0,1), with drift", xlab = "Time years", ylab = "Period index forecasted (Kappa_1)")

Kappa_2_t_CBD_cohort_quad_3 = GLM3_CBD_cohort_quad$coefficients[(n_y +1):(2*n_y)]
Kappa_2_t_Forecasted_CBD_cohort_quad_3 = forecast(Arima(Kappa_2_t_CBD_cohort_quad_3, c(1,0,1)), h = Number_Years_Forecasted)
Kappa_2_t_Forecasted_CBD_cohort_quad_3
plot(Kappa_2_t_Forecasted_CBD_cohort_quad_3,  main = "CBD quad cohort logit forecast kappa_2 ARIMA(1,0,1)", xlab = "Time years", ylab = "Period index forecasted (Kappa_2)")


Kappa_3_t_CBD_cohort_quad_3 = GLM3_CBD_cohort_quad$coefficients[(2*n_y +1):(3*n_y)]
Kappa_3_t_Forecasted_CBD_cohort_quad_3 = forecast(Arima((Kappa_3_t_CBD_cohort_quad_3),c(1,0,1), include.drift = TRUE), h = Number_Years_Forecasted)
Kappa_3_t_Forecasted_CBD_cohort_quad_3
plot(Kappa_3_t_Forecasted_CBD_cohort_quad_3,  main = "CBD quad cohort logit forecast kappa_3 ARIMA(1,0,1), with drift", xlab = "Time years", ylab = "Period index forecasted (Kappa_3)")


Gamma_t_CBD_cohort_quad_3 = GLM3_CBD_cohort_quad$coefficients[(3*n_y +1):(length(GLM3_CBD_cohort_quad$coefficients))]
Gamma_t_Forecasted_CBD_cohort_quad_3 = forecast(Arima(Gamma_t_CBD_cohort_quad_3, c(0,0,0), include.mean = TRUE), h = Number_Years_Forecasted + n_a -1)
Gamma_t_Forecasted_CBD_cohort_quad_3
plot(Gamma_t_Forecasted_CBD_cohort_quad_3,  main = "CBD quad cohort logit forecast gamma ARIMA(0,0,0), with mean", xlab = "Time years", ylab = "Period index forecasted (Gamma)")



X_CBD_cohort_quad_forecast = Model_Matrix_CBD_cohort_quad(n_a, Number_Years_Forecasted)
Teta_CBD_cohort_quad_forecast_3 = c(Kappa_1_t_Forecasted_CBD_cohort_quad_3$mean, Kappa_2_t_Forecasted_CBD_cohort_quad_3$mean,Kappa_3_t_Forecasted_CBD_cohort_quad_3$mean , Gamma_t_Forecasted_CBD_cohort_quad_3$mean)
eta_CBD_cohort_quad_forecast_3 = (X_CBD_cohort_quad_forecast %*% Teta_CBD_cohort_quad_forecast_3)
mu_logit_CBD_cohort_quad_forecast_3 = exp(eta_CBD_cohort_quad_forecast_3)
mu_logit_CBD_cohort_quad_forecast_3 = matrix(mu_logit_CBD_cohort_quad_forecast_3, nrow = n_a, ncol = Number_Years_Forecasted, byrow = FALSE)
Years = (2016:(2016 + Number_Years_Forecasted-1))

### Plots based only on the forecast
#Mortality rate before 2015 and projected Mortality rate after 2015
persp(Age,Year, mu_logit_CBD_cohort_quad_3, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "CBD_cohort_quad_logit")
persp(Age,Years, mu_logit_CBD_cohort_quad_forecast_3, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "CBD_cohort_quad_logit_forecast")
mu_logit_CBD_cohort_quad_total_1 = cbind(mu_logit_CBD_cohort_quad_3, mu_logit_CBD_cohort_quad_forecast_3)
persp(Age,c(Year,Years), mu_logit_CBD_cohort_quad_total_1, theta=-115,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "CBD_cohort_quad_logit_total")




# Plots based on the trends of the forecasted mortality rates
ratio_forecast_mu_logit_CBD_cohort_quad = matrix(nrow = n_a, ncol = Number_Years_Forecasted-1)
for(i in (2:Number_Years_Forecasted))
{
  ratio_forecast_mu_logit_CBD_cohort_quad[,i-1] = mu_logit_CBD_cohort_quad_forecast_3[,i-1]/mu_logit_CBD_cohort_quad_forecast_3[,1]
}

extrapol_mu_logit_CBD_cohort_quad = rowMeans(mu_logit_CBD_cohort_quad_3[, (n_y-2):n_y]) * ratio_forecast_mu_logit_CBD_cohort_quad
Years = (2016:(2016 + Number_Years_Forecasted-2))
persp(Age,Year, mu_logit_CBD_cohort_quad_3, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "CBD_cohort_quad_logit")
persp(Age,Years, extrapol_mu_logit_CBD_cohort_quad, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "CBD_cohort_quad_logit_forecast")
mu_logit_CBD_cohort_quad_total_3 = cbind(mu_logit_CBD_cohort_quad_3, extrapol_mu_logit_CBD_cohort_quad)
persp(Age,c(Year,Years), mu_logit_CBD_cohort_quad_total_3, theta=-115,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "CBD_cohort_quad_logit_total")





##################################################
######          Forecast Lee-Carter        #######    Done
##################################################


##For the GNM3
Kappa_t_LC_3 = GNM3_Lee_Carter$coefficients[(2*n_a+1):(2*n_a + n_y) ]
Kappa_t_Forecasted_LC_3 = forecast(Arima(Kappa_t_LC_3,c(1,0,1), include.drift = TRUE), h = Number_Years_Forecasted)
Kappa_t_Forecasted_LC_3
plot(Kappa_t_Forecasted_LC_3,  main = "LC logit forecast kappa ARIMA(1,0,1), with drift", xlab = "Time years", ylab = "Period index forecasted (Kappa)")

eta_LC_forecast_3 = c()
for( i in (1: (length(Kappa_t_Forecasted_LC_3$mean))))
{
  eta_LC_forecast_3 = c(eta_LC_forecast_3, GNM3_Lee_Carter$coefficients[1:n_a] + 
                          GNM3_Lee_Carter$coefficients[(n_a+1): (n_a +n_a)] * Kappa_t_Forecasted_LC_3$mean[i])
}
eta_LC_forecast_3 = matrix(eta_LC_forecast_3, nrow = n_a, ncol = Number_Years_Forecasted, byrow = FALSE)
mu_logit_LC_forecast_3 = exp(eta_LC_forecast_3)

Years = (2016:(2016 + Number_Years_Forecasted-1))

#Mortality rate before 2015 and projected Mortality rate after 2015
persp(Age,Year, mu_logit_LC_3, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "LC_logit")
persp(Age,Years, mu_logit_LC_forecast_3, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "LC_logit_forecast")
mu_logit_LC_total_3 = cbind(mu_logit_LC_3, mu_logit_LC_forecast_3)
persp(Age,c(Year,Years), mu_logit_LC_total_3, theta=-115,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "LC_logit_total")


# Plots based on the trends of the forecasted mortality rates
ratio_forecast_mu_logit_LC = matrix(nrow = n_a, ncol = Number_Years_Forecasted-1)
for(i in (2:Number_Years_Forecasted))
{
  ratio_forecast_mu_logit_LC[,i-1] = mu_logit_LC_forecast_3[,i-1]/mu_logit_LC_forecast_3[,1]
}

extrapol_mu_logit_LC = rowMeans(mu_logit_LC_3[, (n_y-2):n_y]) * ratio_forecast_mu_logit_LC
Years = (2016:(2016 + Number_Years_Forecasted-2))
persp(Age,Year, mu_logit_LC_3, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "LC_logit")
persp(Age,Years, extrapol_mu_logit_LC, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "LC_logit_forecast")
mu_logit_LC_total_3 = cbind(mu_logit_LC_3, extrapol_mu_logit_LC)
persp(Age,c(Year,Years), mu_logit_LC_total_3, theta=-115,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "LC_logit_total")




########################################################
######          Forecast Renshaw-Haberman        #######    Done
########################################################


##For the GNM2
par(mfrow = c(1,2) )


Kappa_t_1_RH_cloglog
#Kappa_t_Forecasted_RH_2 = forecast(Arima(Kappa_t_RH_2, c(2,1,0)), h = Number_Years_Forecasted)
Kappa_t_Forecasted_RH_2 = forecast(Arima(Kappa_t_1_RH_cloglog, c(1,0,1), include.drift = TRUE), h = Number_Years_Forecasted)

Kappa_t_Forecasted_RH_2
plot(Kappa_t_Forecasted_RH_2,main = "RH cloglog forecast kappa ARIMA(1,0,1), with drift", xlab = "Time years", ylab = "Period index forecasted (Kappa)")

Gamma_t_x_RH_cloglog
#Gamma_t_Forecasted_RH_2 = forecast(Arima(Gamma_t_RH_2, c(2,2,1)), h = Number_Years_Forecasted +n_a -1)
Gamma_t_Forecasted_RH_2 = forecast(Arima(Gamma_t_x_RH_cloglog, c(1,0,1), include.drift = TRUE), h = Number_Years_Forecasted +n_a -1)
Gamma_t_Forecasted_RH_2
plot(Gamma_t_Forecasted_RH_2, main = "RH cloglog forecast gamma ARIMA(1,0,1), with drift", xlab = "Time years", ylab = "Cohort index forecasted (Gamma)")



## with the Poisson assumption, we find
# cloglog link

eta_RH_forecast_2 = matrix( nrow = n_a, ncol = Number_Years_Forecasted, byrow = FALSE)
for(j in 1:Number_Years_Forecasted)
{
  for(i in 1:n_a)
  {
    eta_RH_forecast_2[i,j] = Beta_X_1_RH_cloglog[i] + Beta_X_2_RH_cloglog[i] * Kappa_t_Forecasted_RH_2$mean[j] + Beta_X_3_RH_cloglog[i] * Gamma_t_Forecasted_RH_2$mean[n_a - (i-1) + (j-1)] 
  }
}
qx_cloglog_RH_forecast_2 = 1-exp(-exp(eta_RH_forecast_2))
mu_cloglog_RH_forecast_2 = -log(1-qx_cloglog_RH_forecast_2)
#mu_cloglogit_RH_forecast_3 = -cloglog(1-qx_cloglogit_RH_forecast_3)

Years = (2016:(2016 + Number_Years_Forecasted-1))

#Mortality rate before 2015 and projected Mortality rate after 2015
#par(mfrow = c(1,2) )
persp(Age,Year, mu_cloglog_RH_2, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "RH_cloglog")
persp(Age,Years, mu_cloglog_RH_forecast_2, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "RH_cloglog_forecast")




# Plots based on the trends of the forecasted mortality rates
ratio_forecast_mu_cloglog_RH = matrix(nrow = n_a, ncol = Number_Years_Forecasted-1)
for(i in (2:Number_Years_Forecasted))
{
  ratio_forecast_mu_cloglog_RH[,i-1] = mu_cloglog_RH_forecast_2[,i-1]/mu_cloglog_RH_forecast_2[,1]
}

extrapol_mu_cloglog_RH = rowMeans(mu_cloglog_RH_2[, (n_y-2):n_y]) * ratio_forecast_mu_cloglog_RH
Years = (2016:(2016 + Number_Years_Forecasted-2))
persp(Age,Year, mu_cloglog_RH_2, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "RH_cloglog")
persp(Age,Years, extrapol_mu_cloglog_RH, theta=-35,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "RH_cloglog_forecast")
mu_cloglog_RH_total_2 = cbind(mu_cloglog_RH_2, extrapol_mu_cloglog_RH)
persp(Age,c(Year,Years), mu_cloglog_RH_total_2, theta=-115,phi=15, shade=0.5, ticktype="detailed", xlab = "Age", ylab = "Year", zlab = "force of mortality", main = "RH_cloglog_total")





###########################################################################################
##### Computation of the last table showing the evolution of the different forecasted mu ##
###########################################################################################


extrapol_mu_logit_APC
extrapol_mu_logit_CBD_cohort_quad
extrapol_mu_logit_LC
extrapol_mu_cloglog_RH


summary_forecast = matrix(c(mu_logit_LC_3[30,n_y], 1-exp(-mu_logit_LC_3[30,n_y]),
                             extrapol_mu_logit_LC[30,25], extrapol_mu_logit_LC[30,25]/mu_logit_LC_3[30,n_y], 1-exp(-extrapol_mu_logit_LC[30,25]), 
                             extrapol_mu_logit_LC[30,50], extrapol_mu_logit_LC[30,50]/mu_logit_LC_3[30,n_y], 1-exp(-extrapol_mu_logit_LC[30,50]),

                            mu_logit_CBD_cohort_quad_3[30,n_y], 1-exp(-mu_logit_CBD_cohort_quad_3[30,n_y]),
                             extrapol_mu_logit_CBD_cohort_quad[30,25], extrapol_mu_logit_CBD_cohort_quad[30,25]/mu_logit_CBD_cohort_quad_3[30,n_y], 1-exp(-extrapol_mu_logit_CBD_cohort_quad[30,25]), 
                             extrapol_mu_logit_CBD_cohort_quad[30,50], extrapol_mu_logit_CBD_cohort_quad[30,50]/mu_logit_CBD_cohort_quad_3[30,n_y], 1-exp(-extrapol_mu_logit_CBD_cohort_quad[30,50]),
                             
                            mu_cloglog_RH_2[30,n_y], 1-exp(-mu_cloglog_RH_2[30,n_y]),
                             extrapol_mu_cloglog_RH[30,25], extrapol_mu_cloglog_RH[30,25]/mu_cloglog_RH_2[30,n_y], 1-exp(-extrapol_mu_cloglog_RH[30,25]), 
                             extrapol_mu_cloglog_RH[30,50], extrapol_mu_cloglog_RH[30,50]/mu_cloglog_RH_2[30,n_y], 1-exp(-extrapol_mu_cloglog_RH[30,50])),
                             
         ncol = 8, nrow = 3, byrow = TRUE)

rownames(summary_forecast) = c("LC","CBD_quad_C","RH")
colnames(summary_forecast) = c("mu_2015", "q_2015", "mu_2040", "var_2040/2015", "q_2040", "mu_2065", "var_2065/2015", "q_2065" )

summary_forecast



