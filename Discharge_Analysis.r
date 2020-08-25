##############
# Analysis of river discharge data
# Authors: Alexander Arzt, Sven Kohler, Martin Mächler, Florian Schottmann
# Date: 22.06.2020
##############

# Step 0:
# source functions file for analysis and load libraries
source("Functions.r")

library(sfsmisc)
library(ggplot2)
library(forecast)
library(doParallel)

# Step 1:
# Data Preparation
# 
# Step 1.1
# Preparation of discharge values
# prepare an RDS-file of the following structure: dataframe of 2 columns ("time" - as date, "discharge" - as num)
#
# time | discharge
# ----------------
# 1999-01-01 | 23.4
# 1999-01-02 | 27.1
# 1999-01-03 | 21.0
# ...
# 
# Step 1.2
# Define coordinates of power turbine (longitude/latitude)
#
longitude <- 37.3
latitude <- -0.5


# Step 2: 
# load discharge values
#
riv <- list(d=readRDS("load_your_file_here.rds"), longitude=longitude, latitude=latitude)
check_riv_data(riv)

# Step 2.1 Compatibility of discharge timeline
# 
# The discharge values have to be adjusted accordingly (see notes setion in readme file)
relevant_discharge <- adjust_discharge_to_meteo(riv)

# RUN THIS LINE if the discharge timeline has missing values, the following operation fetches the longest coherent time series.
relevant_discharge <- longest_coh_ts(relevant_discharge)

# Step 3:
# load meteo data. make sure extracted weather files are in your working directory
# 
# load high level data
meteo <- read_meteo()
#
# get meteo data for river-specific geographical region. Define grid granularity. Define lags to include
meteo_rel <- get_ts_all(meteo, 
                        lon=c(36.7, 38), lat=c(0.6,-0.6), range = TRUE, #                 specify lon/lat range to include
                        start_time = as.character(relevant_discharge$time[1]), #          start time of meteo data
                        end_time = as.character(relevant_discharge$time #                 end time of meteo data
                                                [length(relevant_discharge$time)]),
                        pool=TRUE, pool_grid = c(3,3), pool_fun=mean, #                   pooling functionalities
                        make_matrix=T, #                                                  should matrices be computed for further computations
                        lag=TRUE, lags=1:90, lag_pool=T, lag_pool_span = 5) #             lagging functionalities

# Step 4:
# Model fitting
#

# Step 4.1 Linear Model Fit including all variables
#
fit_lm <- river_reg_all(meteo_data = meteo_rel$master_matrix_lagpooled_pooled,
                     dis_data = relevant_discharge, 
                     pooled_lags = TRUE, 
                     max_lag_used=90) # make sure max_lag_used matches the highest lag from the lags argument in get_ts_all()
# get a summary output of the lm fit
summary(fit_lm$reg_fit)
# plot of the lm fit
plot_reg(fit_lm)
# Residual analysis of the lm fit. Collection of plots (Tukey-Anscombe, Normal QQ, Residuals vs. time, ACF, PACF)
plot_res_reg(fit_lm)

# Step 4.2 Determining error structure
#
# look out in the ACF/PACF plots of the residuals for the typical structures of an ARMA process
#
# In absence of a clear ARMA structure or as complement to an existing one: consult auto.arima
#
resid_struct <- auto.arima(resid(fit_lm$reg_fit), 
                           max.p = 10, 
                           max.q = 10, 
                           max.order = 15) # can expand model comlexity when deemed necessary

# if several candidate structures should be tested, run the following variable selection for the candidate models and
# compare the AIC values. If those are similar, compare the predicitve performance of the models in terms of e.g. mean squared error

# define (candidate) ARIMA order
arima_order <- c(1, 0, 6)

# Step 4.3 Variable selection
#
# reduce variable space. This is a computaionally expensive task. To that end the code is parallelized and can run on multiple cores.
# We propose to use a high performance cluster. If not available, let the task run on all cores of your local computer.
#
forw.select <- step.arima.forw(ts=ts(log10(fit_lm$dis_data$discharge)[(fit_lm$max_lag+1):length(fit_lm$dis_data$discharge)]), 
                               order = arima_order, 
                               xreg =meteo_rel$master_matrix_lagpooled_pooled, 
                               verbose = FALSE, 
                               numCores = 48) # Adjust depending on the number of cores you want to use

back.select <- step.arima.backw(ts=ts(log10(fit_lm$dis_data$discharge)[(fit_lm$max_lag+1):length(fit_lm$dis_data$discharge)], start = fit_lm$dis_data$time[1]), 
                                order = arima_order, 
                                xreg =forw.select, 
                                verbose = FALSE, 
                                numCores = 48) # Adjust depending on the number of cores you want to use


# Step 4.4 Final model fit with ARIMA
#
fit_ARIMA <- arima_reg(dis_data=relevant_discharge, meteo_data=back.select, max_lag_used=90, order= arima_order)

# Analysing model
plot_arima(fit_ARIMA)

# Step 4.5 Prediction with meteo data forecast
# 
pred <- arima_pred(model=fit_ARIMA, n.ahead=365, new_meteo=meteo_rel_pred_new) 
#
# meteo_rel_pred_new needs to have same number of rows as n.ahead argument and the  same covariates as back.select (matrix resulting from backward selection)
# It needs to contain the forecasts of the respective covariates.

# plot prediction
#
plot_pred(pred=pred, model=fit_ARIMA, discharge=relevant_discharge, order= arima_order)


# Step 5:
# Duration curve
#
# For comparison: extract true values in the year that was predicted
true_discharge <- relevant_discharge[(nrow(relevant_discharge)-(365)+1):nrow(relevant_discharge),]

# Convert discharge values into duration curve (see notes section of the readme file)
dur_true <- sort(true_discharge$discharge, decreasing=T)

# Plot the duration curve and add the true values
plott_dur(pred, main="Duration curve incl. 95% prediction interval")
lines(dur_true)

# Step 6:
# Power curve
#
# For comparison: compute true power values
true_power <- convert_discharge_to_power(dur_true)

# Plot the power curve and add the true values
plott_power(pred, main="Electricity production")
lines(true_power)

### END












