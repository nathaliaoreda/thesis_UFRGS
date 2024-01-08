rm(list = ls())
##########################
# Load necessary packages
##########################
library(dplyr)
library(ipeadatar)
library(readxl)
library(lubridate)
library(urca)
library(tidyr)
library(tseries)
####################
# Data 
# Code adapted from (Lindenmeyer and Torrent, 2023)
####################
##############################
# 1 - Source file
##############################
setwd("C:/Users/natha/OneDrive/Desktop/TCC/R code")
source("00_functions.R")
dataset <- read_excel("data/dataset.xlsx",
                      col_types = c(
                        "text", "text", "text",
                        "date", "date", "numeric", "text"
                      )
)
# Removed CE12_CUTIND12 - 80 NA values. Included IPCA and EMBI+

UNRATE <- "SEADE12_TDTGSP12"
#INDPRO <- "PAN12_QIIGG12"
IPCA <- "PRECOS12_IPCA12"
SPREAD <- "JPM366_EMBI366"
##############################
# 2 - Get data
##############################

## Acquiring the metadata from each CODE
metadados <- metadata(dataset$codigo)
metadados_spread <- metadata("JPM366_EMBI366")
## Acquiring the values from each CODE
data <- ipeadata(metadados$code)
data_spread <- ipeadata(metadados_spread$code)

## Converting data to wide format
df <- data %>%
  pivot_wider(names_from = "code") %>%
  select(-c(uname, tcode))
df_sorted <- df[order(as.Date(df$date, format = "%Y/%m/%d")), ]
df_spread <- data_spread %>%
  pivot_wider(names_from = "code") %>%
  select(-c(uname, tcode))
# Transaform SPREAD from daily data to monthly data using mean
df_spread <- df_spread %>%
  mutate(year_month = format(date, "%Y-%m")) %>%  # Create a new column with year-month format
  group_by(year_month) %>%                        # Group by year-month
  summarise(JPM366_EMBI366 = mean(JPM366_EMBI366)) %>% # Calculate the average for each month
  mutate(date = as.Date(paste(year_month, "-01", sep = ""))) %>%  # Create a date column for the beginning of each month
  select(date, JPM366_EMBI366) 

# Merge results in the same df
df <- left_join(df_sorted, df_spread, by = "date")

## Filtering by dates
df_filtered <- subset(df, date >= as.Date("1996-01-01") & date < as.Date("2019-06-01"))
df_filtered <- df_filtered  %>% select_if(~ !any(is.na(.)))

## We want to know how to make the data set stationary using ADF Test ##

non_stationary <-  df

## Identifying the columns
tipo <- c(4)
for (i in 2:ncol(non_stationary)) {
  print(i)
  if (sum(non_stationary[, i] < 0, na.rm = T) > 0) { ## if there is any negative value
    
    if (sum(non_stationary[, i] == 0, na.rm = T) > 0) { ## and if there is any zero value
      
      tipo <- append(tipo, 3) # negative and zero
    } else {
      if (sum(non_stationary[, i] < 0, na.rm = T) == dim(non_stationary[, i])[1]) {
        tipo <- append(tipo,5) # only negative values
      } else {
        tipo <- append(tipo, 1) }# some negative values
    }
  } else { ## no negative values
    
    if (sum(non_stationary[, i] == 0, na.rm = T) > 0) { ## no negative values but a zero
      tipo <- append(tipo, 2)
    } else {
      tipo <- append(tipo, 0) ## no negative values and no zero
    }
  }
}


## 0 - no neg nor zero, 1 - any neg values, 2 - only pos but zero, 3 - neg values and zero.
test <- c(4)
for (i in 2:(ncol(non_stationary))) {
  print(i)
  X <- na.exclude(as.matrix(non_stationary[, i]))
  
  k <- 0
  
  j <- 0
  
  status <- "non-stationary"
  
  while (status == "non-stationary") {
    adf_test <- summary(ur.df(X, "none", lags = 12))
    if (k == 2) {
      dale <- colnames(non_stationary)[i]
    }
    # adf.test(X)[4] <= 0.05
    
    if (adf_test@teststat[1] <= adf_test@cval[1, 2] && kpss.test(X, null = "T")$p.value >= 0.05) {
      status <- "stationary"
      
      if (j == 1) {
        k <- 0.5
      }
    } else {
      if (j == 0) {
        X <- preparacao(X, i)
        
        j <- 1
      } else {
        k <- k + 1
        
        if (tipo[i] == 0 | tipo[i] == 2 | tipo[i] == 3 | tipo[i]==1) {
          X <- diff(X)
        }
        if (tipo[i] == 5) {
          X <- cresc_discreto(X)
        }
        
        j <- j + 1
      }
    }
    
    
    #  print(k)
  }
  
  test <- append(test, k)
}


## Assessing transformation value
transformation <- c(-1)

for (i in 2:ncol(non_stationary)) {
  if (tipo[i] == 0 && test[i] == 0) {
    transformation[i] <- 4
  } else if (tipo[i] == 0 && test[i] == 1) {
    transformation[i] <- 5
  } else if (tipo[i] == 0 && test[i] == 0.5) {
    transformation[i] <- 2
  } else if (tipo[i] == 0 && test[i] == 2) {
    transformation[i] <- 6
  } else if (tipo[i] == 1 && test[i] == 0) {
    transformation[i] <- 1
  } else if (tipo[i] == 1 && test[i] == 1) {
    transformation[i] <- 2
  } else if (tipo[i] == 1 && test[i] == 2) {
    transformation[i] <- 3
  } else if (tipo[i] == 2 && test[i] == 0) {
    transformation[i] <- 1
  } else if (tipo[i] == 3 && test[i] == 0) {
    transformation[i] <- 1
  } else if (tipo[i] == 3 & test[i] == 1) {
    transformation[i] <- 2
  } else if (tipo[i] == 2 & test[i] == 1) {
    transformation[i] <- 2
  } else if (tipo[i] == 5 && test[i] == 1) {
    transformation[i] <- 7}
}


month <- t(as.matrix(month(df_filtered$date)))
year <- t(as.matrix(year(df_filtered$date)))
reference <- available_subjects()
metadados <- metadados %>% inner_join(reference)
nome <- names(df_filtered[2:ncol(df_filtered)])
metadados <- metadados[metadados$code %in% nome,]
dataset <- list(
  "data" = as.matrix(df_filtered[2:ncol(df_filtered)]),
  "transformation" = t(as.matrix(transformation[2:ncol(df_filtered)])),
  month = month, year = year, "metadados" = metadados#,
)
dataset$names <- names(as.data.frame(dataset$data))
names(dataset$data) <- NULL


# Now do the transformation
transformed_dataset <- data.frame(matrix(NA, nrow = nrow(dataset$data), ncol = ncol(dataset$data)))
for (i in (1:ncol(dataset$data))){
  transformed_dataset[,i] <- transform_singlestep(dataset$data[,i],dataset$transform[i])
}
names(transformed_dataset) <- dataset$names

# Remove NA values
df <- transformed_dataset[-c(1:2),] %>% select_if(~ !any(is.na(.)))
head(df)
# Save df
save(df, file = paste("data/","df.rda",sep = ""))

