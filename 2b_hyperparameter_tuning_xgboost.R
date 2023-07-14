#!/usr/bin/env Rscript 
args = commandArgs(trailingOnly=TRUE)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Veerendra P Gadekar 
# Email: ic36871@imail.iitm.ac.in; gpveerendra09@gmail.com

# script for hyperparameter tuning of xgboost using grid search

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

library(xgboost)
library(rsample)
library(openxlsx)
library(data.table)
library(tidyverse)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

data_dir <- args[1];
out_dir <- args[2];

log_ga_file <- paste(data_dir, "/extended_input_set/dbscan_train_log_ga.xlsx", sep = "")
log_ga <- subset(read.xlsx(log_ga_file), select = -c(X1));                   
train_file <- paste(data_dir, "/extended_input_set/dbscan_train_data_classic_log_sqrt.xlsx", sep = "")
train_data <- subset(read.xlsx(train_file),  select = -c(X1))
train <- setDT(cbind(log_ga, train_data)); setnames(train, old = "x", new = "logGA")
figures <- paste(out_dir, "/figures/", sep = "")

features <- colnames(train_data)
# Create the treatment plan from the training data
treatplan <- vtreat::designTreatmentsZ(train, features, verbose = FALSE)

train_data <- subset(read.xlsx(train_file),  select = -c(X1))
train = setDT(cbind(log_ga, train_data)); setnames(train, old = "x", new = "logGA")

# Get the "clean" variable names from the scoreFrame
new_vars <- treatplan %>%
  magrittr::use_series(scoreFrame) %>%        
  dplyr::filter(code %in% c("clean", "lev")) %>% 
  magrittr::use_series(varName)   

# Prepare the training data
features_train <- vtreat::prepare(treatplan, train, varRestriction = new_vars) %>% as.matrix()
response_train <- train$logGA

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# create hyperparameter grid
hyper_grid <- expand.grid(
  eta = c( .1, .3),
  max_depth = c(5, 7),
  min_child_weight = c( 5, 7),
  subsample = c( .8, 1), 
  colsample_bytree = c(.9, 1),
  optimal_trees = 0,               # a place to dump results
  min_RMSE = 0                     # a place to dump results
)

# grid search 
for (i in 1:nrow(hyper_grid)) {
  params <- list(
    eta = hyper_grid$eta[i],
    max_depth = hyper_grid$max_depth[i],
    min_child_weight = hyper_grid$min_child_weight[i],
    subsample = hyper_grid$subsample[i],
    colsample_bytree = hyper_grid$colsample_bytree[i]
  )
  
  # for reproducibility
  set.seed(123)
  
  # train model
  xgb.tune <- xgb.cv(
    params = params,
    data = features_train,
    label = response_train,
    nrounds = 5000,
    nfold = 5,
    objective = "reg:linear",  # for regression models
    verbose = 0,               # silent,
    early_stopping_rounds = 10 # stop if no improvement for 10 consecutive trees
  )
  
  # add min training error and trees to grid
  hyper_grid$optimal_trees[i] <- which.min(xgb.tune$evaluation_log$test_rmse_mean)
  hyper_grid$min_RMSE[i] <- min(xgb.tune$evaluation_log$test_rmse_mean)
}

hyper_grid <- hyper_grid %>% dplyr::arrange(min_RMSE)

write.table(hyper_grid, paste(out_dir,"hyper_grid_xgboost.tsv", sep = "/"), sep = "\t", row.names = F, quote = F)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

