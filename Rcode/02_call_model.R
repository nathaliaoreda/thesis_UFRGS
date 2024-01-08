##KEY FEATURES FOR IMPROVING MACROECONOMIC FORECASTING -EVIDENCE FOR BRAZIL
##October 8th, 2023
#Nathalia Scholles Oreda

rm(list = ls())

# Import packages
library(glmnet)
library(randomForest)
library(dplyr)
library(mboost)
library(e1071)

# Variables of interest

UNRATE <- "SEADE12_TDTGSP12"
IPCA <- "PRECOS12_IPCA12"
SPREAD <- "JPM366_EMBI366"

####################
# Get data
####################
setwd("C:/Users/natha/OneDrive/Desktop/TCC/R code")
source("00_functions.R")
load("data/df.Rda")
df= df%>%as.matrix()
df=as.matrix(scale(df))
head(df)

####################
# Select parameters
####################

variable <-  c(UNRATE,IPCA, SPREAD)
hyperparameter_list <- c("bic", "aic","cv")
window <- c("rolling", "expanding")
horizon_list <- c(1, 3, 6, 9, 12)
lag_orders <- c(1:12)
factor_orders <- c(3, 6, 10)
regularization_list <- c("boost","ridge","lasso","en")

####### run window ##########

# Set seed
set.seed(1234)

model_name = "rf_rich"
model_function = rf_rich
model_list = list()
forecasts_list = list()
start_window=155# "2009-01-01"
end_window=nrow(df) #279 "2019-05-01"
nwindows = end_window - start_window +1 #125
actual_values=df[start_window:end_window]

window = c("expanding")
hyperparameter_list <- c("cv")
regularization_list <- c(NA)

for(v in variable){
  for(h in horizon_list){
    for(hyper in hyperparameter_list){
      for (regu in regularization_list){
        # The train set should be h steps behind the test set
        window_sup <- start_window - h - 1
        for (i in 1:nwindows){
          reoptimize_hyperparameters=FALSE
          if (window=="expanding"){
            Df=df[1:(window_sup+i),]
          }
          if (window=="rolling"){
            Df=df[i:(window_sup+i),]
          }
          # Check if it's time to re-optimize hyperparameters (the first iteration and every year)
          if (i %% 12 == 1) {
            reoptimize_hyperparameters=TRUE
          }
          model = model_function(Df,horizon=h,variable=v,lag_orders,type=hyper,reoptimize_hyperparameters,best,regu=regu)
          model_list[[i]] <- list(pred=model$pred)
          best <- model$best
        }
        #All forecasts should have same test size
        forecast=head(unlist(lapply(model_list,function(x)x$pred)),nwindows)
        if (is.na(regu)==TRUE){
          forecasts_list[[length(forecasts_list)+1]] <- list(variable=v,horizon=h,model=paste(model_name,hyper,sep="_"),pred=forecast)
          }
        if (is.na(regu)==FALSE){
            forecasts_list[[length(forecasts_list)+1]] <- list(variable=v,horizon=h,model=paste(model_name,hyper,regu,sep="_"),pred=forecast)  
        }
      }
    }
    cat(h,"\n")
  }
  cat(v,"\n")
}

saveRDS(forecasts_list,file = paste("forecasts/",model_name,".rda",sep = ""))
#saveRDS(actual_values, file= paste("forecasts/","actual_values",".rda",sep = ""))


    