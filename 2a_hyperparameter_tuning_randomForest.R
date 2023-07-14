#!/usr/bin/env Rscript 
args = commandArgs(trailingOnly=TRUE)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Veerendra P Gadekar 
# Email: ic36871@imail.iitm.ac.in; gpveerendra09@gmail.com

# script for hyperparameter tuning of RandomForest using grid search

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

library(openxlsx)
library(data.table)
require(ranger) #for Random Forest model
library(randomForest)
library(tidyverse)
library(gridExtra)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

seed_value <- 666
seed_range <- 100:599

data_dir <- args[1];
out_dir <- args[2];

log_ga_file <- paste(data_dir, "/extended_input_set/dbscan_train_log_ga.xlsx", sep = "")
log_ga <- subset(read.xlsx(log_ga_file), select = -c(X1));                   
train_file <- paste(data_dir, "/extended_input_set/dbscan_train_data_classic_log_sqrt.xlsx", sep = "")
train_data <- subset(read.xlsx(train_file),  select = -c(X1))
train <- setDT(cbind(log_ga, train_data)); setnames(train, old = "x", new = "logGA")
figures <- paste(out_dir, "/figures/", sep = "")

hyper_grid <- expand.grid(
  mtry       = seq(3, 15, by = 1),
  node_size  = seq(3, 10, by = 1),
  sample_size = c(.55, .632, .70, .80),
  OOB_RMSE   = 0
)

RFmodel <- randomForest(
  formula = logGA ~ .,
  data    = train
)

pdf(paste(figures,"Out-Of-Bag_MSE_across_ntrees.pdf", sep = ""))
plot(RFmodel, main = "Out-Of-Bag error plot Random Forest")
# we select num.trees 100 because 
# the performance gains with adding more trees is very less
# as we have only a few features to work with
abline(v = 100, lty = 2, col = "red")
dev.off()

model <- 
 apply(hyper_grid, 1, 
   function(y){      
      ranger(
        formula         = logGA ~ ., 
        data            = train, 
        num.trees       = 100,
        mtry            = y["mtry"],
        min.node.size   = y["node_size"],
        sample.fraction = y["sample_size"],
        seed            = 123)
    })

hyper_grid$OOB_RMSE <- unlist(lapply(model, function(x) x$prediction.error))
hyper_grid <- setDT(hyper_grid)[order(OOB_RMSE)]

# distribution of RMSE values in 100 iterations
# values of mtry, min.node.size, sample.fraction resulting the lowest OOB RMSE is
# taken from the results from the hyper grid from the above step

print("Selected Hyperparameters:")
print(paste("mtry: ", hyper_grid$mtry[1], sep = ""))
print(paste("min.node.size: ", hyper_grid$node_size[1], sep = ""))
print(paste("sample.fraction: ", hyper_grid$sample_size[1], sep = ""))

OOB_RMSE <- vector(mode = "numeric", length = 100)
for (i in seq_along(OOB_RMSE)) {
  optimal_ranger <- ranger(
    formula         = logGA ~ ., 
    data            = train, 
    num.trees       = 100,
    mtry            = hyper_grid$mtry[1],
    min.node.size   = hyper_grid$node_size[1],
    sample.fraction = hyper_grid$sample_size[1],
    importance      = 'impurity'
  )
  OOB_RMSE[i] <- sqrt(optimal_ranger$prediction.error)
}
pdf(paste(figures,"RMSE_distribution_100_iteration_RF.pdf", sep = ""))
hist(OOB_RMSE,breaks = 20)
dev.off()

write.table(hyper_grid, paste(out_dir,"hyper_grid_RF.tsv", sep = "/"), sep = "\t", row.names = F, quote = F)
write.table(train, paste(data_dir,"RF_XGB_train_data.tsv", sep = "/"), sep = "\t", row.names = F, quote = F)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 

