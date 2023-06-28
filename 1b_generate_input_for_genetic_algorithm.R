#!/usr/bin/env Rscript 
args = commandArgs(trailingOnly=TRUE)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Veerendra P Gadekar 
# Email: ic36871@imail.iitm.ac.in; gpveerendra09@gmail.com

# script to generate input for Genetic Algorithm, Random Forest and XGBOOST

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# library

library(xlsx)
library(data.table)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# change the path according to the working directory

#data_dir    <- '/home/veer/scratch/IBSE/GarbhIni2/Garbhini-GA2-main/development/data/';
#out_dir <- '/home/veer/scratch/IBSE/GarbhIni2/Garbhini-GA2-main/development/data/input_for_gapolyfitn_matlab/';
data_dir    <- args[1];
out_dir <- args[2];


# input files
dbscan_train <- paste(data_dir, 'train_23_dbscan.tsv', sep = "");
dbscan_train = fread(dbscan_train)

# ga
ga <- dbscan_train$ga
write.xlsx(ga, paste(out_dir, "dbscan_train_ga.xlsx", sep = ""))
log_ga <- log(dbscan_train$ga)
write.xlsx(log_ga, paste(out_dir,"dbscan_train_log_ga.xlsx", sep = ""))


#classic
dbscan_train_data_classic = subset(dbscan_train, select = -c(enrid, ga_birth, ga))
dbscan_train_data_classic_bpd = subset(dbscan_train_data_classic, select = -c(bpd))
dbscan_train_data_classic_ofd = subset(dbscan_train_data_classic, select = -c(ofd))
dbscan_train_data_classic_hp  = subset(dbscan_train_data_classic, select = -c(hp))
dbscan_train_data_classic_ap  = subset(dbscan_train_data_classic, select = -c(ap))
dbscan_train_data_classic_fl  = subset(dbscan_train_data_classic, select = -c(fl))

write.xlsx(dbscan_train_data_classic, paste(out_dir, "dbscan_train_data_classic.xlsx", sep = ""))
write.xlsx(dbscan_train_data_classic_bpd, paste(out_dir, "dbscan_train_data_classic_bpd.xlsx", sep = ""))
write.xlsx(dbscan_train_data_classic_ofd, paste(out_dir, "dbscan_train_data_classic_ofd.xlsx", sep = ""))
write.xlsx(dbscan_train_data_classic_hp, paste(out_dir, "dbscan_train_data_classic_hp.xlsx", sep = ""))
write.xlsx(dbscan_train_data_classic_ap, paste(out_dir, "dbscan_train_data_classic_ap.xlsx", sep = ""))
write.xlsx(dbscan_train_data_classic_fl, paste(out_dir, "dbscan_train_data_classic_fl.xlsx", sep = ""))


#log
dbscan_train_data_log = as.data.table(apply(dbscan_train_data_classic, 2, log))
colnames(dbscan_train_data_log) = paste("log_", colnames(dbscan_train_data_log), sep = "")
dbscan_train_data_log_bpd = subset(dbscan_train_data_log, select = -c(log_bpd))
dbscan_train_data_log_ofd = subset(dbscan_train_data_log, select = -c(log_ofd))
dbscan_train_data_log_hp  = subset(dbscan_train_data_log, select = -c(log_hp))
dbscan_train_data_log_ap  = subset(dbscan_train_data_log, select = -c(log_ap))
dbscan_train_data_log_fl  = subset(dbscan_train_data_log, select = -c(log_fl))

write.xlsx(dbscan_train_data_log, paste(out_dir, "dbscan_train_data_log.xlsx", sep = ""))
write.xlsx(dbscan_train_data_log_bpd, paste(out_dir, "dbscan_train_data_log_bpd.xlsx", sep = ""))
write.xlsx(dbscan_train_data_log_ofd, paste(out_dir, "dbscan_train_data_log_ofd.xlsx", sep = ""))
write.xlsx(dbscan_train_data_log_hp, paste(out_dir, "dbscan_train_data_log_hp.xlsx", sep = ""))
write.xlsx(dbscan_train_data_log_ap, paste(out_dir, "dbscan_train_data_log_ap.xlsx", sep = ""))
write.xlsx(dbscan_train_data_log_fl, paste(out_dir, "dbscan_train_data_log_fl.xlsx", sep = ""))


#sqrt
dbscan_train_data_sqrt = as.data.table(apply(dbscan_train_data_classic, 2, sqrt))
colnames(dbscan_train_data_sqrt) = paste("sqrt_", colnames(dbscan_train_data_sqrt), sep = "")
dbscan_train_data_sqrt_bpd = subset(dbscan_train_data_sqrt, select = -c(sqrt_bpd))
dbscan_train_data_sqrt_ofd = subset(dbscan_train_data_sqrt, select = -c(sqrt_ofd))
dbscan_train_data_sqrt_hp  = subset(dbscan_train_data_sqrt, select = -c(sqrt_hp))
dbscan_train_data_sqrt_ap  = subset(dbscan_train_data_sqrt, select = -c(sqrt_ap))
dbscan_train_data_sqrt_fl  = subset(dbscan_train_data_sqrt, select = -c(sqrt_fl))

write.xlsx(dbscan_train_data_sqrt, paste(out_dir, "dbscan_train_data_sqrt.xlsx", sep = ""))
write.xlsx(dbscan_train_data_sqrt_bpd, paste(out_dir, "dbscan_train_data_sqrt_bpd.xlsx", sep = ""))
write.xlsx(dbscan_train_data_sqrt_ofd, paste(out_dir, "dbscan_train_data_sqrt_ofd.xlsx", sep = ""))
write.xlsx(dbscan_train_data_sqrt_hp, paste(out_dir, "dbscan_train_data_sqrt_hp.xlsx", sep = ""))
write.xlsx(dbscan_train_data_sqrt_ap, paste(out_dir, "dbscan_train_data_sqrt_ap.xlsx", sep = ""))
write.xlsx(dbscan_train_data_sqrt_fl, paste(out_dir, "dbscan_train_data_sqrt_fl.xlsx", sep = ""))


#classic + log
dbscan_train_data_classic_log = cbind(dbscan_train_data_classic, dbscan_train_data_log)
dbscan_train_data_classic_log_bpd = cbind(dbscan_train_data_classic_bpd, dbscan_train_data_log_bpd)
dbscan_train_data_classic_log_ofd = cbind(dbscan_train_data_classic_ofd, dbscan_train_data_log_ofd)
dbscan_train_data_classic_log_hp = cbind(dbscan_train_data_classic_hp, dbscan_train_data_log_hp)
dbscan_train_data_classic_log_ap = cbind(dbscan_train_data_classic_ap, dbscan_train_data_log_ap)
dbscan_train_data_classic_log_fl = cbind(dbscan_train_data_classic_fl, dbscan_train_data_log_fl)

write.xlsx(dbscan_train_data_classic_log, paste(out_dir, "dbscan_train_data_classic_log.xlsx", sep = ""))
write.xlsx(dbscan_train_data_classic_log_bpd, paste(out_dir, "dbscan_train_data_classic_log_bpd.xlsx", sep = ""))
write.xlsx(dbscan_train_data_classic_log_ofd, paste(out_dir, "dbscan_train_data_classic_log_ofd.xlsx", sep = ""))
write.xlsx(dbscan_train_data_classic_log_hp, paste(out_dir, "dbscan_train_data_classic_log_hp.xlsx", sep = ""))
write.xlsx(dbscan_train_data_classic_log_ap, paste(out_dir, "dbscan_train_data_classic_log_ap.xlsx", sep = ""))
write.xlsx(dbscan_train_data_classic_log_fl, paste(out_dir, "dbscan_train_data_classic_log_fl.xlsx", sep = ""))


#classic + log + sqrt
dbscan_train_data_classic_log_sqrt = cbind(dbscan_train_data_classic, dbscan_train_data_log, dbscan_train_data_sqrt)
dbscan_train_data_classic_log_sqrt_bpd = cbind(dbscan_train_data_classic_bpd, dbscan_train_data_log_bpd, dbscan_train_data_sqrt_bpd)
dbscan_train_data_classic_log_sqrt_ofd = cbind(dbscan_train_data_classic_ofd, dbscan_train_data_log_ofd, dbscan_train_data_sqrt_ofd)
dbscan_train_data_classic_log_sqrt_hp = cbind(dbscan_train_data_classic_hp, dbscan_train_data_log_hp, dbscan_train_data_sqrt_hp)
dbscan_train_data_classic_log_sqrt_ap = cbind(dbscan_train_data_classic_ap, dbscan_train_data_log_ap, dbscan_train_data_sqrt_ap)
dbscan_train_data_classic_log_sqrt_fl = cbind(dbscan_train_data_classic_fl, dbscan_train_data_log_fl, dbscan_train_data_sqrt_fl)

write.xlsx(dbscan_train_data_classic_log_sqrt, paste(out_dir, "dbscan_train_data_classic_log_sqrt.xlsx", sep = ""))
write.xlsx(dbscan_train_data_classic_log_sqrt_bpd, paste(out_dir, "dbscan_train_data_classic_log_sqrt_bpd.xlsx", sep = ""))
write.xlsx(dbscan_train_data_classic_log_sqrt_ofd, paste(out_dir, "dbscan_train_data_classic_log_sqrt_ofd.xlsx", sep = ""))
write.xlsx(dbscan_train_data_classic_log_sqrt_hp, paste(out_dir, "dbscan_train_data_classic_log_sqrt_hp.xlsx", sep = ""))
write.xlsx(dbscan_train_data_classic_log_sqrt_ap, paste(out_dir, "dbscan_train_data_classic_log_sqrt_ap.xlsx", sep = ""))
write.xlsx(dbscan_train_data_classic_log_sqrt_fl, paste(out_dir, "dbscan_train_data_classic_log_sqrt_fl.xlsx", sep = ""))
	   
