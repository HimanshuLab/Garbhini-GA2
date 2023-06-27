#!/usr/bin/env Rscript 
args = commandArgs(trailingOnly=TRUE)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Veerendra P Gadekar 
# Email: ic36871@imail.iitm.ac.in; gpveerendra09@gmail.com

# script to remove the outliers from the data using DBSCAN (density-based spatial clustering) method 
# only the training set with removed outlier is used in the downstream analysis

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# library

library(data.table)
library(haven)
library(tidyverse)
library(dbscan)
library(readxl)
library(StatMatch)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# user inputs

data_dir    <- args[1]; # path to the directory containing data
results_dir <- args[2]; # path to the directory where results are to be stored

# input files

train_file <- paste(data_dir, 'train_23.dta', sep = "");
test_file  <- paste(data_dir, 'test_23.dta', sep = "");

train   <- read_dta(train_file) %>% select(-sfh, -visit)
test    <- read_dta(test_file) %>% select(-sfh, -visit)
test_m  <- apply(as.matrix(test)[,2:6], 2, as.numeric)
train_m <- apply(as.matrix(train)[,2:6], 2, as.numeric)

# dbscan flag the outliers

pdf(paste(results_dir, "/figures/DBSCAN_train_knn.pdf", sep = ""), h=10, w=10)
  kNNdistplot(train_m, minPts = 5, k = 5)
  abline(h=0.7, col = "red", lty = 2)
  train_res <- dbscan(train_m, eps = 0.7, minPts = 5)
  train$dbscan_outlier = train_res$cluster
  pairs(train_m, col = train_res$cluster + 1L)
dev.off()

pdf(paste(results_dir, "/figures/DBSCAN_test_knn.pdf", sep=""), h=10, w=10)
  kNNdistplot(test_m, minPts = 5, k = 5)
  abline(h=0.7, col = "red", lty = 2)
  test_res <- dbscan(test_m, eps = 0.7, minPts = 5)
  test$dbscan_outlier = test_res$cluster
  pairs(test_m, col = test_res$cluster + 1L)
dev.off()

train1 = subset(subset(train, dbscan_outlier != 0), select = -c(dbscan_outlier))
test1 = subset(subset(test, dbscan_outlier != 0), select = -c(dbscan_outlier))

write.table(test1,  paste(data_dir, 'test_23_dbscan.tsv', sep = ""), sep = "\t", row.names = F, quote = F)
write.table(train1, paste(data_dir, 'train_23_dbscan.tsv', sep = ""), sep = "\t", row.names = F, quote = F)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# prepare the CMCV validation data set
# batch1

cmcv_data_part1  <- read_excel(paste(data_dir, "Dataset_for_validation_of_dating_models.xlsx", sep = ""))
cmcv_data_part1  <- as.data.table(cmcv_data_part1)
cmcv_data_part1a <- subset(cmcv_data_part1, 
                      select = c(recordid, v12, v13, v22, 
		                 biparietaldiameterbpdmeasurement, occipitofrontaldiameterofdmeasur, 
                                 headperimetercircumferencehpmeas, abdominalperimetercircumferencea, 
                                 femurlengthflmeasurementfrombest))

# found some wrong entries in the date of dating scan - correcting them here
cmcv_data_part1a$v13 <- gsub("2023", "2022", cmcv_data_part1a$v13)

# batch2

cmcv_data_part2  <- fread(paste(data_dir, "DevelopmentAndPerfor_DATA_2023-01-19_0634_VG.csv", sep = ""))
cmcv_data_part2  <- subset(cmcv_data_part2, !recordid %in% cmcv_data_part1$recordid)
cmcv_data_part2a <- subset(cmcv_data_part2, 
                      select = c(recordid, v12, v13, v22, 
				 biparietaldiameterbpdmeasurement, occipitofrontaldiameterofdmeasur, 
				 headperimetercircumferencehpmeas, abdominalperimetercircumferencea, 
				 femurlengthflmeasurementfrombest))

cmcv_data_part2a$v13 <- format(strptime(as.Date(cmcv_data_part2a$v13), "%d-%m-%Y"), "%Y-%m-%d")
cmcv_data_part2a$v22 <- format(strptime(as.Date(cmcv_data_part2a$v22), "%d-%m-%Y"), "%Y-%m-%d")

# concatenate batch1 and batch2

cmcv_data <- rbind(setDT(cmcv_data_part1a), cmcv_data_part2a)

setnames(cmcv_data, 
  c("enrid", "crl", "date_dating_scan" , "date_2_trimester_scan", "bpd",   "ofd",    "hp",    "ap",    "fl"))
cmcv_data_df1 <- subset(cmcv_data, crl < 100)

# compute gold standard GA using Garbhini-GA1 formula based on CRL measurements

cmcv_data_df1$t_weeks <- as.numeric(as.character((as.Date(cmcv_data_df1$date_2_trimester_scan)-as.Date(cmcv_data_df1$date_dating_scan))/7))
setDT(cmcv_data_df1)[, crl := crl/10]
cmcv_data_df1$Garbhini_GA1 <- (6.73526 + 1.15018*(cmcv_data_df1$crl) -0.02294*(cmcv_data_df1$crl^2)) # formula from GA1
cmcv_data_df1 <- subset(cmcv_data_df1, Garbhini_GA1 >= 7 & Garbhini_GA1 < 14)
cmcv_data_df1[, gold_standard_GA := Garbhini_GA1 + t_weeks]
cmcv_data_df1 <- subset(cmcv_data_df1, gold_standard_GA <= 40)
setnames(cmcv_data_df1, "gold_standard_GA", "ga")

features <- c('fl','ofd','bpd','hp','ap')
cmcv_test <- cmcv_data_df1[,c('enrid',features,'ga', 't_weeks'), with = FALSE]
cmcv_test_m <- as.matrix(cmcv_test[,2:6])

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# dbscan flag the outliers

pdf(paste(results_dir, "/figures/DBSCAN_CMCV_knn_NEW.pdf", sep=""), h=10, w=10)
kNNdistplot(cmcv_test_m, minPts = 5, k = 5)
abline(h=1.3, col = "red", lty = 2)
cmcv_test_res <- dbscan(cmcv_test_m, eps = 1.3, minPts = 5)
cmcv_test$dbscan_outlier = cmcv_test_res$cluster
pairs(cmcv_test_m, col = cmcv_test_res$cluster + 1L)
dev.off()

cmcv_test1 = subset(subset(cmcv_test, dbscan_outlier != 0), select = -c(dbscan_outlier))
write.table(cmcv_test1, paste(data_dir, 'cmcv_validation_dbscan.tsv', sep = ""), sep = "\t", row.names = F, quote = F)
write.table(cmcv_test, paste(data_dir, 'cmcv_validation.tsv', sep = ""), sep = "\t", row.names = F, quote = F)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

