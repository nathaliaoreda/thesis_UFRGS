rm(list = ls())

library(tidyverse)
library(forecast)
library(MCS)
library(xtable)
library(tidyr)
library(ggplot2)
library(ggsci)
library(ggpubr)

setwd("C:/Users/natha/OneDrive/Desktop/TCC/R code")
actual_values <- readRDS("forecasts/actual_values.rda")
# Variables of interest

UNRATE <- "SEADE12_TDTGSP12"
IPCA <- "PRECOS12_IPCA12"
SPREAD <- "JPM366_EMBI366"

variable <-  c(UNRATE, IPCA, SPREAD)
test_matrix <- as.matrix(actual_values)
horizon_list <- c(1, 3, 6, 9, 12)
rmspe_data <- data.frame(variable = character(),
                         horizon = character(),
                         model = character(),
                         RMSPE = numeric(),
                         Diebold_Mariano_Test = numeric(),
                         stringsAsFactors = FALSE)


# ###########################
# # Run this code to save AR-BIC results###
# # This is the benchmark used for comparison
# ###########################
# model_list <- readRDS("forecasts/ar_bic_rw.rda")
# for (v in variable) {
#   for (h in horizon_list) {
#     filtered_models <- model_list[sapply(model_list, function(x) x$horizon == h & x$variable == v & x$model == "ar_bic_rw_bic")]
#     model_matrix <- sapply(filtered_models, function(x) x$pred)
#     pred_error_bic <- (test_matrix - model_matrix)
#     rmspe_bic <- sqrt(mean(pred_error_bic[is.finite(pred_error_bic)]^2, na.rm = TRUE))
#     rmspe_data <- rbind(rmspe_data, data.frame(variable = v, horizon = h, model="ar_bic", RMSPE = rmspe_bic))
#   }
# }
# 
# # Print the resulting data frame
# print(rmspe_data)
# print(pred_error_bic)
# saveRDS(rmspe_data,file = paste("forecasts/","ar_bic_results_rw",".rda",sep = ""))
# saveRDS(pred_error_bic,file = paste("forecasts/","pred_error_bic_rw",".rda",sep = ""))


############################
##### For other forecasts###
############################

# Data poor
ar <- readRDS("forecasts/ar_rw.rda")
shrink_poor <- readRDS("forecasts/shrink_poor_rw.rda")
shrink_poor_en <- readRDS("forecasts/shrink_poor_en_rw.rda")
rf_poor <- readRDS("forecasts/rf_poor_rw.rda")
bols_poor <- readRDS("forecasts/bols_poor_rw.rda")
bbs_poor <- readRDS("forecasts/bbs_poor_rw.rda")
svr_linear_poor <- readRDS("forecasts/svr_linear_poor_rw.rda")
svr_rbf_poor <- readRDS("forecasts/svr_rbf_poor_rw.rda")
# Data rich
ardi <- readRDS("forecasts/ardi_rw.rda")
shrink_rich <- readRDS("forecasts/shrink_rich_rw.rda")
shrink_rich_en <- readRDS("forecasts/shrink_rich_en_rw.rda")
rf_rich <- readRDS("forecasts/rf_rich_rw.rda")
bols_rich <- readRDS("forecasts/bols_rich_rw.rda")
bbs_rich <- readRDS("forecasts/bbs_rich_rw.rda")
svr_linear_rich <- readRDS("forecasts/svr_linear_rich_rw.rda")
svr_rbf_rich <- readRDS("forecasts/svr_rbf_rich_rw.rda")

# Bs
b1 <- readRDS("forecasts/b1_rw.rda")
b1_en <- readRDS("forecasts/b1_en_rw.rda")
b2 <- readRDS("forecasts/b2_rw.rda")
b2_en <- readRDS("forecasts/b2_en_rw.rda")
b3 <- readRDS("forecasts/b3_rw.rda")
b3_en <- readRDS("forecasts/b3_en_rw.rda")

model_list <- c(ar,shrink_poor,shrink_poor_en,rf_poor,bols_poor,bbs_poor,svr_linear_poor,svr_rbf_poor,
                ardi,shrink_rich,shrink_rich_en,rf_rich,bols_rich,bbs_rich,svr_linear_rich,svr_rbf_rich,
                b1,b1_en,b2,b2_en,b3,b3_en)

models_list <- c("ar_rw_aic","ar_rw_cv","shrink_poor_rw_cv_ridge","shrink_poor_rw_cv_lasso","shrink_poor_en_2_rw_cv_en","rf_poor_rw_cv","bols_poor_rw_cv","bbs_poor_rw_cv","svr_linear_poor_rw_cv","svr_rbf_poor_rw_cv",
                 "ardi_rw_bic","ardi_rw_aic","ardi_rw_cv","shrink_rich_rw_cv_ridge","shrink_rich_rw_cv_lasso","shrink_rich_en_2_rw_cv_en","rf_rich_rw_cv","bols_rich_rw_cv","bbs_rich_rw_cv","svr_linear_rich_rw_cv","svr_rbf_rich_rw_cv",
                 "b1_rw_cv_boost","b1_rw_cv_ridge","b1_en_rw_cv_en","b1_rw_cv_lasso","b2_rw_cv_boost","b2_rw_cv_ridge","b2_en_rw_cv_en","b2_rw_cv_lasso","b3_rw_cv_boost","b3_rw_cv_ridge","b3_en_rw_cv_en","b3_rw_cv_lasso")

ar_bic_rmspe <- readRDS("forecasts/ar_bic_results_rw.rda")
pred_error_bic_rw <- readRDS("forecasts/pred_error_bic_rw.rda")
matrix_list <- list()
CSFE <- list()
CSFE_list <- list()
squared_error_list <- list()

for (v in variable) {
  for (h in horizon_list) {
    for (model in models_list){
      model_list_filtered <- model_list[sapply(model_list, function(x) x$horizon == h & x$variable == v & x$model == model)]
      model_matrix <- sapply(model_list_filtered, function(x) x$pred)
      pred_error <- (test_matrix - model_matrix)
      matrix_list[[model]] <- pred_error[is.finite(pred_error)]^2
      CSFE[[model]] <- cumsum((pred_error_bic_rw[is.finite(pred_error_bic_rw)]^2) - (pred_error[is.finite(pred_error)]^2))
      rmspe <- (sqrt(mean(pred_error[is.finite(pred_error)]^2, na.rm = TRUE)) 
                / ar_bic_rmspe$RMSPE[ar_bic_rmspe$variable == v & ar_bic_rmspe$horizon == h])
      # Diebold-Mariano test
      dm_test <- forecast::dm.test(pred_error, pred_error_bic_rw,h=h,varestimator = "bartlett")
      rmspe_data <- rbind(rmspe_data, data.frame(variable = v, horizon = h,model=model,RMSPE = rmspe,Diebold_Mariano_Test = format(dm_test$p.value, scientific = FALSE)))
    }
    # Store the squared error values in results_list
    result_matrix <- do.call(cbind, matrix_list)
    CSFE_matrix <- do.call(cbind, CSFE)
    squared_error_list[[paste(v,h,sep="_")]] <- result_matrix
    CSFE_list[[paste(v,h,sep="_")]] <- CSFE_matrix
  }
}

# Print the resulting data frame
results<-list(rmspe_data=rmspe_data,squared_error=squared_error_list,CSFE=CSFE_list)
str(results)


### MCS ###
# Horizon 1
SPREAD_1 <- results$squared_error[["JPM366_EMBI366_1"]]
columns_to_exclude <- c("shrink_poor_rw_cv_lasso","b2_rw_cv_lasso", "b3_rw_cv_lasso")
SPREAD_1 <- SPREAD_1[, -which(colnames(SPREAD_1) %in% columns_to_exclude)]
IPCA_1 <- results$squared_error[["PRECOS12_IPCA12_1"]]
UNRATE_1 <- results$squared_error[["SEADE12_TDTGSP12_1"]]
MCS_SPREAD_1 <- MCSprocedure(SPREAD_1)
MCS_IPCA_1 <- MCSprocedure(IPCA_1)
MCS_UNRATE_1 <- MCSprocedure(UNRATE_1)

# Horizon 3
SPREAD_3 <- results$squared_error[["JPM366_EMBI366_3"]]
columns_to_exclude <- c("shrink_poor_rw_cv_lasso","b2_rw_cv_lasso", "b3_rw_cv_lasso")
SPREAD_3 <- SPREAD_3[, -which(colnames(SPREAD_3) %in% columns_to_exclude)]
IPCA_3 <- results$squared_error[["PRECOS12_IPCA12_3"]]
UNRATE_3 <- results$squared_error[["SEADE12_TDTGSP12_3"]]
MCS_SPREAD_3 <- MCSprocedure(SPREAD_3)
MCS_IPCA_3 <- MCSprocedure(IPCA_3)
MCS_UNRATE_3 <- MCSprocedure(UNRATE_3)

# Horizon 6
SPREAD_6 <- results$squared_error[["JPM366_EMBI366_6"]]
columns_to_exclude <- c("shrink_poor_rw_cv_lasso","b2_rw_cv_lasso", "b3_rw_cv_lasso","ardi_rw_aic")
SPREAD_6 <- SPREAD_6[, -which(colnames(SPREAD_6) %in% columns_to_exclude)]
IPCA_6 <- results$squared_error[["PRECOS12_IPCA12_6"]]
columns_to_exclude <- c("b2_rw_cv_lasso", "b3_rw_cv_lasso","ardi_rw_aic")
IPCA_6 <- IPCA_6[, -which(colnames(IPCA_6) %in% columns_to_exclude)]
UNRATE_6 <- results$squared_error[["SEADE12_TDTGSP12_6"]]
columns_to_exclude <- c("shrink_poor_rw_cv_lasso","b2_rw_cv_lasso", "b3_rw_cv_lasso","ardi_rw_aic")
UNRATE_6 <- UNRATE_6[, -which(colnames(UNRATE_6) %in% columns_to_exclude)]
MCS_SPREAD_6 <- MCSprocedure(SPREAD_6)
MCS_IPCA_6 <- MCSprocedure(IPCA_6)
MCS_UNRATE_6 <- MCSprocedure(UNRATE_6)

# Horizon 9
SPREAD_9 <- results$squared_error[["JPM366_EMBI366_9"]]
columns_to_exclude <- c("shrink_poor_rw_cv_lasso","b2_rw_cv_lasso", "b3_rw_cv_lasso","ardi_rw_aic")
SPREAD_9 <- SPREAD_9[, -which(colnames(SPREAD_9) %in% columns_to_exclude)]
IPCA_9 <- results$squared_error[["PRECOS12_IPCA12_9"]]
columns_to_exclude <- c("b2_rw_cv_lasso", "b3_rw_cv_lasso","ardi_rw_aic")
IPCA_9 <- IPCA_9[, -which(colnames(IPCA_9) %in% columns_to_exclude)]
UNRATE_9 <- results$squared_error[["SEADE12_TDTGSP12_9"]]
columns_to_exclude <- c("shrink_poor_rw_cv_lasso","b2_rw_cv_lasso", "b3_rw_cv_lasso","ardi_rw_aic")
UNRATE_9 <- UNRATE_9[, -which(colnames(UNRATE_9) %in% columns_to_exclude)]
MCS_SPREAD_9 <- MCSprocedure(SPREAD_9)
MCS_IPCA_9 <- MCSprocedure(IPCA_9)
MCS_UNRATE_9 <- MCSprocedure(UNRATE_9)

# Horizon 12
SPREAD_12 <- results$squared_error[["JPM366_EMBI366_12"]]
columns_to_exclude <- c("shrink_poor_rw_cv_lasso","b2_rw_cv_lasso", "b3_rw_cv_lasso","ardi_rw_aic")
SPREAD_12 <- SPREAD_12[, -which(colnames(SPREAD_12) %in% columns_to_exclude)]
IPCA_12 <- results$squared_error[["PRECOS12_IPCA12_12"]]
columns_to_exclude <- c("b2_rw_cv_lasso", "b3_rw_cv_lasso","ardi_rw_aic")
IPCA_12 <- IPCA_12[, -which(colnames(IPCA_12) %in% columns_to_exclude)]
UNRATE_12 <- results$squared_error[["SEADE12_TDTGSP12_12"]]
columns_to_exclude <- c("shrink_poor_rw_cv_lasso","b3_rw_cv_lasso","ardi_rw_aic")
UNRATE_12<- UNRATE_12[, -which(colnames(UNRATE_12) %in% columns_to_exclude)]
MCS_SPREAD_12 <- MCSprocedure(SPREAD_12)
MCS_IPCA_12 <- MCSprocedure(IPCA_12)
MCS_UNRATE_12 <- MCSprocedure(UNRATE_12)

############################
##### Store results########
############################

get_stars <- function(p_values) {
  significance <- ifelse(is.na(p_values), "",
                  ifelse(p_values < 0.01, "***",
                  ifelse(p_values < 0.05, "**",
                  ifelse(p_values < 0.1, "*", ""))))
  return(significance)
}

### RMSPE ###
rmspe<-results$rmspe_data
rmspe<-bind_rows(ar_bic_rmspe,rmspe_data)%>%arrange(variable,horizon)
rmspe$stars<-get_stars(rmspe$Diebold_Mariano_Test)
rmspe$RMSPE<-round(rmspe$RMSPE,4)
rmspe_SPREAD <- rmspe %>%
  filter(variable == SPREAD)%>%
  unite("values", RMSPE, stars, sep = "", remove = FALSE)%>%
  select(model, horizon, values) %>%
  pivot_wider(names_from = horizon, names_prefix = "h=",values_from = values, values_fill = NULL)
rmspe_IPCA <- rmspe %>%
  filter(variable == IPCA)%>%
  unite("values", RMSPE, stars, sep = "", remove = FALSE)%>%
  select(model, horizon, values) %>%
  pivot_wider(names_from = horizon, names_prefix = "h=",values_from = values, values_fill = NULL)
rmspe_UNRATE <- rmspe %>%
  filter(variable == UNRATE)%>%
  unite("values", RMSPE, stars, sep = "", remove = FALSE)%>%
  select(model, horizon, values) %>%
  pivot_wider(names_from = horizon, names_prefix = "h=",values_from = values, values_fill = NULL)

print(xtable(rmspe_SPREAD,type='latex'))
print(xtable(rmspe_IPCA,type='latex'))
print(xtable(rmspe_UNRATE,type='latex'))

### PLOTS ###

### SPREAD ###
SPREAD_1_CSFE <- results$CSFE[["JPM366_EMBI366_1"]]
SPREAD_3_CSFE <- results$CSFE[["JPM366_EMBI366_3"]]
SPREAD_6_CSFE <- results$CSFE[["JPM366_EMBI366_6"]]
SPREAD_9_CSFE <- results$CSFE[["JPM366_EMBI366_9"]]
SPREAD_12_CSFE <- results$CSFE[["JPM366_EMBI366_12"]]

# Get best models according to RANK_M of MCS
SPREAD_1_CSFE_df <- SPREAD_1_CSFE %>%
  as.data.frame()%>%
  mutate(dates = seq(as.Date("2009-01-01"), by = "month", length.out = 125))%>%
  select(b1_rw_cv_lasso,b3_rw_cv_boost,svr_rbf_rich_rw_cv,b2_rw_cv_boost,shrink_poor_rw_cv_ridge,ar_rw_cv,b2_rw_cv_ridge,bbs_poor_rw_cv,b3_rw_cv_ridge,rf_poor_rw_cv,dates) %>%
  gather(key = "Model", value = "value", -dates)
SPREAD_3_CSFE_df <- SPREAD_3_CSFE %>%
  as.data.frame()%>%
  mutate(dates = seq(as.Date("2009-01-01"), by = "month", length.out = 125))%>%
  select(svr_linear_poor_rw_cv, b1_rw_cv_lasso,bbs_poor_rw_cv, bols_poor_rw_cv, svr_rbf_rich_rw_cv, bbs_rich_rw_cv, svr_linear_rich_rw_cv,shrink_rich_rw_cv_lasso, svr_rbf_poor_rw_cv,b2_rw_cv_ridge,dates) %>%
  gather(key = "Model", value = "value", -dates)
SPREAD_6_CSFE_df <- SPREAD_6_CSFE %>%
  as.data.frame()%>%
  mutate(dates = seq(as.Date("2009-01-01"), by = "month", length.out = 125))%>%
  select(svr_linear_poor_rw_cv, b1_rw_cv_lasso, bols_poor_rw_cv, bols_rich_rw_cv,ar_rw_aic,svr_rbf_rich_rw_cv,bbs_rich_rw_cv,bbs_poor_rw_cv, shrink_rich_rw_cv_lasso,b3_rw_cv_ridge,dates) %>%
  gather(key = "Model", value = "value", -dates)
SPREAD_9_CSFE_df <- SPREAD_9_CSFE %>%
  as.data.frame()%>%
  mutate(dates = seq(as.Date("2009-01-01"), by = "month", length.out = 125))%>%
  select(ar_rw_aic,b1_rw_cv_lasso,bols_poor_rw_cv, b1_rw_cv_boost, bbs_rich_rw_cv,b2_rw_cv_ridge,svr_rbf_rich_rw_cv,shrink_rich_rw_cv_lasso,shrink_poor_en_2_rw_cv_en ,b3_rw_cv_ridge,dates) %>%
  rename(shrink_poor_rw_cv_en = shrink_poor_en_2_rw_cv_en) %>%
  gather(key = "Model", value = "value", -dates)
SPREAD_12_CSFE_df <- SPREAD_12_CSFE %>%
  as.data.frame()%>%
  mutate(dates = seq(as.Date("2009-01-01"), by = "month", length.out = 125))%>%
  select(bols_poor_rw_cv, shrink_poor_rw_cv_ridge, b1_rw_cv_lasso, bbs_rich_rw_cv, b1_rw_cv_boost, svr_rbf_rich_rw_cv,b2_rw_cv_boost,b3_rw_cv_ridge, shrink_poor_en_2_rw_cv_en,shrink_rich_rw_cv_lasso,dates) %>%
  rename(shrink_poor_rw_cv_en = shrink_poor_en_2_rw_cv_en) %>%
  gather(key = "Model", value = "value", -dates)


SPREAD_1_plot <- ggplot(SPREAD_1_CSFE_df, aes(x=dates,y= value)) + 
  geom_line(aes(color = Model),linetype = "solid", size=0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_npg()+
  labs(x = "Time", title = "h=1")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 15))
SPREAD_3_plot <- ggplot(SPREAD_3_CSFE_df, aes(x=dates,y= value)) + 
  geom_line(aes(color = Model),linetype = "solid",size=0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_npg()+
  labs(x = "Time", title = "h=3")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 15))
SPREAD_6_plot <- ggplot(SPREAD_6_CSFE_df, aes(x=dates,y= value)) + 
  geom_line(aes(color = Model),linetype = "solid",size=0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_npg()+
  labs(x = "Time", title = "h=6")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 15))
SPREAD_9_plot <- ggplot(SPREAD_9_CSFE_df, aes(x=dates,y= value)) + 
  geom_line(aes(color = Model),linetype = "solid",size=0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_npg()+
  labs(x = "Time", title = "h=9")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 15))
SPREAD_12_plot <- ggplot(SPREAD_9_CSFE_df, aes(x=dates,y= value)) + 
  geom_line(aes(color = Model),linetype = "solid",size=0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_npg()+
  labs(x = "Time", title = "h=12")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 15))

plot_spread=ggarrange(SPREAD_1_plot , SPREAD_3_plot , SPREAD_6_plot , SPREAD_9_plot, SPREAD_12_plot, 
                      ncol = 2, nrow = 3)


### IPCA ###
IPCA_1_CSFE <- results$CSFE[["PRECOS12_IPCA12_1"]]
IPCA_3_CSFE <- results$CSFE[["PRECOS12_IPCA12_3"]]
IPCA_6_CSFE <- results$CSFE[["PRECOS12_IPCA12_6"]]
IPCA_9_CSFE <- results$CSFE[["PRECOS12_IPCA12_9"]]
IPCA_12_CSFE <- results$CSFE[["PRECOS12_IPCA12_12"]]

# Get best models according to RANK_M of MCS
IPCA_1_CSFE_df <- IPCA_1_CSFE %>%
  as.data.frame()%>%
  mutate(dates = seq(as.Date("2009-01-01"), by = "month", length.out = 125))%>%
  select(b1_rw_cv_lasso,b3_rw_cv_lasso,b2_rw_cv_lasso,b2_en_rw_cv_en,b3_rw_cv_ridge,b2_rw_cv_ridge,svr_rbf_rich_rw_cv,b1_en_rw_cv_en, shrink_poor_en_2_rw_cv_en ,shrink_poor_rw_cv_lasso, dates) %>%
  rename(b2_rw_cv_en = b2_en_rw_cv_en) %>%
  rename(b1_rw_cv_en = b1_en_rw_cv_en) %>%
  rename(shrink_poor_rw_cv_en = shrink_poor_en_2_rw_cv_en) %>%
  gather(key = "Model", value = "value", -dates)
IPCA_3_CSFE_df <- IPCA_3_CSFE %>%
  as.data.frame()%>%
  mutate(dates = seq(as.Date("2009-01-01"), by = "month", length.out = 125))%>%
  select(b3_rw_cv_boost,b1_rw_cv_boost,b2_rw_cv_boost,b3_rw_cv_ridge,b3_en_rw_cv_en,b1_en_rw_cv_en,b1_rw_cv_lasso,b2_rw_cv_ridge,ar_rw_cv, svr_rbf_rich_rw_cv,dates) %>%
  rename(b3_rw_cv_en = b3_en_rw_cv_en) %>%
  rename(b1_rw_cv_en = b1_en_rw_cv_en) %>%
  gather(key = "Model", value = "value", -dates)
IPCA_6_CSFE_df <- IPCA_6_CSFE %>%
  as.data.frame()%>%
  mutate(dates = seq(as.Date("2009-01-01"), by = "month", length.out = 125))%>%
  select(b3_en_rw_cv_en,b1_rw_cv_lasso,b1_en_rw_cv_en,b2_rw_cv_boost,svr_rbf_rich_rw_cv,b3_rw_cv_boost,bbs_rich_rw_cv,ar_rw_cv,b1_rw_cv_boost,bols_rich_rw_cv,dates) %>%
  rename(b3_rw_cv_en = b3_en_rw_cv_en) %>%
  rename(b1_rw_cv_en = b1_en_rw_cv_en) %>%
  gather(key = "Model", value = "value", -dates)
IPCA_9_CSFE_df <- IPCA_9_CSFE %>%
  as.data.frame()%>%
  mutate(dates = seq(as.Date("2009-01-01"), by = "month", length.out = 125))%>%
  select(b1_rw_cv_lasso,b3_en_rw_cv_en,b1_en_rw_cv_en,b2_rw_cv_boost,svr_rbf_rich_rw_cv,b3_rw_cv_boost,ar_rw_cv,bbs_rich_rw_cv,svr_rbf_poor_rw_cv,bols_rich_rw_cv,dates) %>%
  rename(b3_rw_cv_en = b3_en_rw_cv_en) %>%
  rename(b1_rw_cv_en = b1_en_rw_cv_en) %>%
  gather(key = "Model", value = "value", -dates)
IPCA_12_CSFE_df <- IPCA_12_CSFE %>%
  as.data.frame()%>%
  mutate(dates = seq(as.Date("2009-01-01"), by = "month", length.out = 125))%>%
  select(b1_rw_cv_lasso,b1_en_rw_cv_en,b3_en_rw_cv_en,b2_rw_cv_boost,b3_rw_cv_boost,svr_rbf_rich_rw_cv,bols_rich_rw_cv,shrink_rich_rw_cv_ridge,svr_rbf_poor_rw_cv,bbs_rich_rw_cv,dates) %>%
  rename(b3_rw_cv_en = b3_en_rw_cv_en) %>%
  rename(b1_rw_cv_en = b1_en_rw_cv_en) %>%
  gather(key = "Model", value = "value", -dates)


IPCA_1_plot <- ggplot(IPCA_1_CSFE_df, aes(x=dates,y= value)) + 
  geom_line(aes(color = Model),linetype = "solid", size=0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_npg()+
  labs(x = "Time", title = "h=1")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 15))
IPCA_3_plot <- ggplot(IPCA_3_CSFE_df, aes(x=dates,y= value)) + 
  geom_line(aes(color = Model),linetype = "solid",size=0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_npg()+
  labs(x = "Time", title = "h=3")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 15))
IPCA_6_plot <- ggplot(IPCA_6_CSFE_df, aes(x=dates,y= value)) + 
  geom_line(aes(color = Model),linetype = "solid",size=0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_npg()+
  labs(x = "Time", title = "h=6")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 15))
IPCA_9_plot <- ggplot(IPCA_9_CSFE_df, aes(x=dates,y= value)) + 
  geom_line(aes(color = Model),linetype = "solid",size=0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_npg()+
  labs(x = "Time", title = "h=9")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 15))
IPCA_12_plot <- ggplot(IPCA_9_CSFE_df, aes(x=dates,y= value)) + 
  geom_line(aes(color = Model),linetype = "solid",size=0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_npg()+
  labs(x = "Time", title = "h=12")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 15))

plot_IPCA=ggarrange(IPCA_1_plot , IPCA_3_plot , IPCA_6_plot , IPCA_9_plot, IPCA_12_plot, 
                    ncol = 2, nrow = 3)

### UNRATE ###

UNRATE_1_CSFE <- results$CSFE[["SEADE12_TDTGSP12_1"]]
UNRATE_3_CSFE <- results$CSFE[["SEADE12_TDTGSP12_3"]]
UNRATE_6_CSFE <- results$CSFE[["SEADE12_TDTGSP12_6"]]
UNRATE_9_CSFE <- results$CSFE[["SEADE12_TDTGSP12_9"]]
UNRATE_12_CSFE <- results$CSFE[["SEADE12_TDTGSP12_12"]]

# Get best models according to RANK_M of MCS
UNRATE_1_CSFE_df <- UNRATE_1_CSFE %>%
  as.data.frame()%>%
  mutate(dates = seq(as.Date("2009-01-01"), by = "month", length.out = 125))%>%
  select(b3_rw_cv_lasso, b3_en_rw_cv_en, b3_rw_cv_ridge, b2_rw_cv_lasso, b2_rw_cv_ridge, b1_rw_cv_lasso, svr_rbf_rich_rw_cv, shrink_rich_en_2_rw_cv_en, shrink_rich_rw_cv_lasso, shrink_poor_rw_cv_lasso,dates) %>%
  rename(b3_rw_cv_en = b3_en_rw_cv_en) %>%
  rename(shrink_rich_rw_cv_en = shrink_rich_en_2_rw_cv_en) %>%
  gather(key = "Model", value = "value", -dates)
UNRATE_3_CSFE_df <- UNRATE_3_CSFE %>%
  as.data.frame()%>%
  mutate(dates = seq(as.Date("2009-01-01"), by = "month", length.out = 125))%>%
  select(b3_rw_cv_lasso, b3_en_rw_cv_en, b3_rw_cv_ridge, b2_rw_cv_lasso, b2_en_rw_cv_en, b2_rw_cv_ridge, b1_rw_cv_lasso, shrink_rich_en_2_rw_cv_en, shrink_rich_rw_cv_lasso, shrink_poor_rw_cv_lasso,dates) %>%
  rename(b3_rw_cv_en = b3_en_rw_cv_en) %>%
  rename(b2_rw_cv_en = b2_en_rw_cv_en) %>%
  rename(shrink_rich_rw_cv_en = shrink_rich_en_2_rw_cv_en) %>%
  gather(key = "Model", value = "value", -dates)
UNRATE_6_CSFE_df <- UNRATE_6_CSFE %>%
  as.data.frame()%>%
  mutate(dates = seq(as.Date("2009-01-01"), by = "month", length.out = 125))%>%
  select(b3_en_rw_cv_en, b3_rw_cv_ridge, b3_rw_cv_boost, b2_en_rw_cv_en , b2_rw_cv_ridge, b2_rw_cv_boost, b1_rw_cv_lasso, svr_rbf_rich_rw_cv, shrink_rich_rw_cv_lasso, shrink_poor_en_2_rw_cv_en,dates) %>%
  rename(b3_rw_cv_en = b3_en_rw_cv_en) %>%
  rename(b2_rw_cv_en = b2_en_rw_cv_en) %>%
  rename(shrink_poor_rw_cv_en = shrink_poor_en_2_rw_cv_en) %>%
  gather(key = "Model", value = "value", -dates)
UNRATE_9_CSFE_df <- UNRATE_9_CSFE %>%
  as.data.frame()%>%
  mutate(dates = seq(as.Date("2009-01-01"), by = "month", length.out = 125))%>%
  select(b3_en_rw_cv_en , b3_rw_cv_ridge , b2_en_rw_cv_en, b2_rw_cv_ridge, b2_rw_cv_boost, b1_rw_cv_lasso , svr_rbf_rich_rw_cv  , shrink_rich_rw_cv_lasso, shrink_poor_en_2_rw_cv_en , shrink_poor_rw_cv_ridge,dates) %>%
  rename(b3_rw_cv_en = b3_en_rw_cv_en) %>%
  rename(b2_rw_cv_en = b2_en_rw_cv_en) %>%
  rename(shrink_poor_rw_cv_en = shrink_poor_en_2_rw_cv_en) %>%
  gather(key = "Model", value = "value", -dates)
UNRATE_12_CSFE_df <- UNRATE_12_CSFE %>%
  as.data.frame()%>%
  mutate(dates = seq(as.Date("2009-01-01"), by = "month", length.out = 125))%>%
  select(b3_en_rw_cv_en , b3_rw_cv_ridge, b2_rw_cv_lasso, b2_en_rw_cv_en , b2_rw_cv_ridge , b1_rw_cv_lasso, svr_rbf_rich_rw_cv, shrink_rich_en_2_rw_cv_en, shrink_rich_rw_cv_lasso , shrink_poor_en_2_rw_cv_en,dates) %>%
  rename(b3_rw_cv_en = b3_en_rw_cv_en) %>%
  rename(b2_rw_cv_en = b2_en_rw_cv_en) %>%
  rename(shrink_rich_rw_cv_en = shrink_rich_en_2_rw_cv_en) %>%
  rename(shrink_poor_rw_cv_en = shrink_poor_en_2_rw_cv_en) %>%
  gather(key = "Model", value = "value", -dates)


UNRATE_1_plot <- ggplot(UNRATE_1_CSFE_df, aes(x=dates,y= value)) + 
  geom_line(aes(color = Model),linetype = "solid", size=0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_npg()+
  labs(x = "Time", title = "h=1")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 15))
UNRATE_3_plot <- ggplot(UNRATE_3_CSFE_df, aes(x=dates,y= value)) + 
  geom_line(aes(color = Model),linetype = "solid",size=0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_npg()+
  labs(x = "Time", title = "h=3")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 15))
UNRATE_6_plot <- ggplot(UNRATE_6_CSFE_df, aes(x=dates,y= value)) + 
  geom_line(aes(color = Model),linetype = "solid",size=0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_npg()+
  labs(x = "Time", title = "h=6")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 15))
UNRATE_9_plot <- ggplot(UNRATE_9_CSFE_df, aes(x=dates,y= value)) + 
  geom_line(aes(color = Model),linetype = "solid",size=0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_npg()+
  labs(x = "Time", title = "h=9")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 15))
UNRATE_12_plot <- ggplot(UNRATE_12_CSFE_df, aes(x=dates,y= value)) + 
  geom_line(aes(color = Model),linetype = "solid",size=0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_npg()+
  labs(x = "Time", title = "h=12")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 15))

plot_UNRATE=ggarrange(UNRATE_1_plot , UNRATE_3_plot , UNRATE_6_plot , UNRATE_9_plot, UNRATE_12_plot, 
                      ncol = 2, nrow = 3)
