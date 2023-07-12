#!/usr/bin/env Rscript 
args = commandArgs(trailingOnly=TRUE)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Veerendra P Gadekar 
# Email: ic36871@imail.iitm.ac.in; gpveerendra09@gmail.com

# script to test the performance of all the GA generated formulas on the Garbhini Test-set and 
# CMCV validation set 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# library

library(data.table)
library(caret)
library(openxlsx)
library(data.table)
library(ranger) #for Random Forest model
library(xgboost)
library(tidyverse)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

infile1 <- args[1]
infile2 <- args[2]
infile3 <- args[3]
infile4 <- args[4]
infile5 <- args[5]
outdir  <- args[6]

genetic_algo <- fread(infile1, sep = "\t")
genetic_algo <- genetic_algo[order(RMSE)]

dbscan_train <- fread(infile2)
cmpt_test  <- fread(infile3)
cmpt_cmcv  <- fread(infile4)
cmpt_cmcv  <- subset(cmpt_cmcv, !enrid %in% c(955, 909, 60)) # negative ga in this case; not sure why; removed
rf_xgb_train     <- fread(infile5)
dat <- dbscan_train

data_test <- cmpt_cmcv
out <-
lapply(unique(genetic_algo$formula),
   function(f){
     model <- lm(f, data = dat);

     formula <- f;
     ga <- data_test$ga;
     t_weeks <- data_test$t_weeks;
     ga_predicted <- if(grepl("^log", f)){
	               exp(predict(model, data_test))}else{predict(model, data_test)};
     error <- ga_predicted - ga;
     eq  <- paste0(round(coefficients(model)[1],5), " + ", 
              paste(sprintf("%.5f * %s", coefficients(model)[-1],  
                names(coefficients(model)[-1])), collapse=" + "));
     RMSE <- sqrt(mean((ga - ga_predicted)^2));
     R2 <- R2(predict(model, data_test),ga);
     AIC <- AIC(model);
     BIC <- BIC(model);
     out <- data.table(
              enrid = data_test$enrid, 
              fl =  data_test$fl, 
              ofd  = data_test$ofd,  
              bpd = data_test$bpd,  
              hp = data_test$hp, 
              ap =data_test$ap,
	          formula = formula, 
	          ga = ga, 
	          ga_predicted = ga_predicted, 
	          error = error, 
	          t_weeks = t_weeks, 
	          RMSE = RMSE, 
	          R2 = R2, 
	          AIC = AIC, 
	          BIC = BIC, 
	          eq = eq)
 })

perf_cmcv <- rbindlist(out)
perf_cmcv <- perf_cmcv[order(RMSE)]
perf_cmcv[, `:=` (q75=quantile(error, c(0.75)), 
		  q50=quantile(error, c(0.50)),
		  q25=quantile(error, c(0.25)),
		  mean_ga= mean(ga),
		  ga= ga,		  
		  mean_ga_predicted = mean(ga_predicted),
		  median_error = median(error),
		  mean_error = mean(error)), by = formula]

test <- cmpt_cmcv
predict_ga <-
  function(test, method){
    ga <- numeric()
    hp <- test$hp
    fl <- test$fl
    ap <- test$ap
    bpd <- test$bpd
    ofd <- test$ofd
    if("Hadlock" %in% method) {
        ga <- 10.85 + (0.06*hp*fl) + (0.67*bpd) + (0.168*ap)
      }else if("Intergrowth" %in% method){
        ga <- exp((0.03243*(log(hp*10)^2)) + (0.001644*(fl*10)*log(hp*10)) + 3.813)/7
     }
   }  

seed_value <- 666
seed_range <- 100:599

cmpt_cmcv_RF_XGB = copy(cmpt_cmcv)
cmpt_cmcv_RF_XGB[, `:=`(log_bpd=log(bpd),  log_ofd=log(ofd),   log_hp=log(hp),   log_ap=log(ap), log_fl=log(fl), 
                   sqrt_bpd=sqrt(bpd), sqrt_ofd=sqrt(ofd),  sqrt_hp=sqrt(hp), sqrt_ap =sqrt(ap),  sqrt_fl=sqrt(fl))]

ranger_model <- ranger(
  formula         = logGA ~ ., 
  data            = rf_xgb_train,   
  num.trees       = 100,
  mtry            = 3,
  min.node.size   = 10,
  sample.fraction = .55,
  importance = "permutation", 
  local.importance = TRUE)


features <- setdiff(colnames(rf_xgb_train) , "logGA")
# Create the treatment plan from the training data
treatplan <- vtreat::designTreatmentsZ(rf_xgb_train, features, verbose = FALSE)

# Get the "clean" variable names from the scoreFrame
new_vars <- treatplan %>%
  magrittr::use_series(scoreFrame) %>%        
  dplyr::filter(code %in% c("clean", "lev")) %>% 
  magrittr::use_series(varName)   

# Prepare the training data
features_train <- vtreat::prepare(treatplan, rf_xgb_train, varRestriction = new_vars) %>% as.matrix()
response_train <- rf_xgb_train$logGA
# Prepare the test data
features_test <- vtreat::prepare(treatplan, subset(cmpt_cmcv_RF_XGB, select = -c(enrid, t_weeks, dbscan_outlier, ga)), varRestriction = new_vars) %>% as.matrix()
response_test <- subset(cmpt_cmcv_RF_XGB, select = -c(enrid, t_weeks, dbscan_outlier, ga))$logGA

# xgboost
# parameter list
params <- list(
  eta = 0.1,
  max_depth = 5,
  min_child_weight = 7,
  subsample = 1,
  colsample_bytree = 1
)
# train final model
xgb.fit.final <- xgboost(
  params = params,
  data = features_train,
  label = response_train,
  nrounds = 77,
  objective = "reg:linear",
  verbose = 0
)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

test$hadlock <- predict_ga(test,"Hadlock")
test$intergrowth <- predict_ga(test,"Intergrowth")
test$randomforest <- exp(predictions(predict(ranger_model,subset(cmpt_cmcv_RF_XGB, select = -c(enrid, t_weeks, ga, dbscan_outlier)))))
test$xgboost <- exp(predict(xgb.fit.final, features_test))

hadlock <- subset(test, select = -c(intergrowth, randomforest, xgboost))
setnames(hadlock, "hadl	ock", "ga_predicted")
hadlock$eq  = 'ga = 10.85 + 0.060(hp)(fl) + 0.67(bpd) + 0.168(ap)'
hadlock$formula = "Hadlock"
hadlock$error <- hadlock$ga_predicted - hadlock$ga
hadlock$median_error = median(hadlock$error)
hadlock$AIC <- AIC(lm("log(ga) ~ log(hp) * log(fl) + log(bpd) + log(ap)", data = dat))
hadlock$BIC <- BIC(lm("log(ga) ~ log(hp) * log(fl) + log(bpd) + log(ap)", data = dat))
hadlock$RMSE = sqrt(mean((hadlock$ga - hadlock$ga_predicted)^2))
hadlock$R2 <- R2(predict(lm("log(ga) ~ hp * fl + bpd + ap", data = dat), data_test),hadlock$ga);
hadlock[, `:=` (
		  median_error = median(error),
		  mean_error = mean(error)), by = formula]

intergrowth <- subset(test,select = -c(hadlock, randomforest, xgboost))
setnames(intergrowth, "intergrowth", "ga_predicted")
intergrowth$eq  = 'log(ga) = 0.03243 × (log(hp))^2 + 0.001644 × fl × log(hp) + 3.813'
intergrowth$formula = "Intergrowth"
intergrowth$error <- intergrowth$ga_predicted - intergrowth$ga
intergrowth$median_error = median(intergrowth$error)
intergrowth$AIC <- AIC(lm("log(ga) ~ (log(hp))^2 + fl + log(hp)", data = dat))
intergrowth$BIC <- BIC(lm("log(ga) ~ (log(hp))^2 + fl + log(hp)", data = dat))
intergrowth$RMSE = sqrt(mean((intergrowth$ga - intergrowth$ga_predicted)^2))
intergrowth$R2 <- R2(predict(lm("log(ga) ~ hp * fl + bpd + ap", data = dat), data_test),intergrowth$ga);
intergrowth[, `:=` (
		  median_error = median(error),
		  mean_error = mean(error)), by = formula]


randomforest <- subset(test,select = -c(hadlock, intergrowth, xgboost))
setnames(randomforest, "randomforest", "ga_predicted")
randomforest$eq  = "RandomForest"
randomforest$formula = "RandomForest"
randomforest$error <- randomforest$ga_predicted - randomforest$ga
randomforest$median_error = median(randomforest$error)
randomforest$RMSE = RMSE(exp(predictions(predict(ranger_model,subset(cmpt_cmcv_RF_XGB, select = -c(enrid, t_weeks,dbscan_outlier,  ga))))), cmpt_cmcv_RF_XGB$ga) 
randomforest$R2 = cor(exp(predictions(predict(ranger_model,subset(cmpt_cmcv_RF_XGB, select = -c(enrid, t_weeks, dbscan_outlier, ga))))),cmpt_cmcv_RF_XGB$ga)^2
randomforest[, `:=` (
		  median_error = median(error),
		  mean_error = mean(error)), by = formula]

xgboost <- subset(test,select = -c(hadlock, intergrowth, randomforest))
setnames(xgboost, "xgboost", "ga_predicted")
xgboost$eq  = "GradientBoosting"
xgboost$formula = "GradientBoosting"
xgboost$error <- xgboost$ga_predicted - xgboost$ga
xgboost$median_error = median(xgboost$error)
xgboost$RMSE = RMSE(exp(predict(xgb.fit.final, features_test)), cmpt_cmcv_RF_XGB$ga) 
xgboost$R2 = cor(exp(predict(xgb.fit.final, features_test)),cmpt_cmcv_RF_XGB$ga)^2
xgboost[, `:=` (
		  median_error = median(error),
		  mean_error = mean(error)), by = formula]

kn <- rbindlist(list(hadlock, intergrowth, randomforest, xgboost), fill = TRUE)

perf_cmcv = rbindlist(list(setDT(kn), perf_cmcv), fill=TRUE)
perf_cmcv2 <- unique(subset(perf_cmcv, select = -c(enrid, ga, ga_predicted, error, t_weeks)))
perf_cmcv2 <- perf_cmcv2[order(abs(median_error))]
perf_cmcv$formula1 <- factor(perf_cmcv$formula, levels = rev(unique(perf_cmcv2$formula)))

perf_cmcv_df <- unique(subset(perf_cmcv, select = c(eq, formula, RMSE, R2, AIC, BIC, mean_error, median_error)))
perf_cmcv_df <- perf_cmcv_df[order(abs(median_error))]

plot1 <- ggplot(perf_cmcv, aes(y=error, x= formula1))+
  geom_violin(trim= FALSE, alpha = 0.5, fill = "red") +
  geom_boxplot(width = 0.2,  outlier.size = .5) +
  ggtitle("CMCV Validation-set")+
  xlab("error(in weeks)")+
  scale_y_continuous(breaks = -8:8, limits = c(-8,8)) + theme_bw() +
  labs(caption = paste("n = ",dim(data_test)[1])) +   xlab("") +
  geom_hline(aes(yintercept=0), linetype="dashed", size=0.05, colour="black") +
  theme(legend.position = "none") +
  coord_flip()

pdf(paste(outdir, "CMCV_validation_error_plot.pdf", sep = ""),  h=9, w=12)
plot1
dev.off()

write.table(subset(perf_cmcv, select = -c(formula1)), 
  paste(outdir, "CMCV_validation_error.tsv", sep = ""), 
  sep = "\t", row.names = F, quote = F)

write.table(perf_cmcv_df, 
  paste(outdir, "CMCV_validation_error_table.tsv", sep = ""), 
  sep = "\t", row.names = F, quote = F)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

data_test <- cmpt_test
out <-
lapply(genetic_algo$formula,
   function(f){
     model <- lm(f, data = dat);

     formula <- f;
     ga <- data_test$ga;
     ga_predicted <- if(grepl("^log", f)){
                       exp(predict(model, data_test))}else{predict(model, data_test)};
     error <- ga_predicted - ga;
     eq  <- paste0(round(coefficients(model)[1],5), " + ",
              paste(sprintf("%.5f * %s", coefficients(model)[-1],
                names(coefficients(model)[-1])), collapse=" + "));
     RMSE <- sqrt(mean((ga - ga_predicted)^2));
     R2 <- R2(predict(model, data_test),ga);
     AIC <- AIC(model);
     BIC <- BIC(model);
     out <- data.table(
              enrid = data_test$enrid,
              fl =  data_test$fl, 
              ofd  = data_test$ofd,  
              bpd = data_test$bpd,  
              hp = data_test$hp, 
              ap = data_test$ap,
              formula = formula,
              ga = ga,
              ga_predicted = ga_predicted,
              error = error,
              RMSE = RMSE,
              R2 = R2,
              AIC = AIC,
              BIC = BIC,
              eq = eq)
    })

perf_test <- rbindlist(out)
perf_test <- perf_test[order(RMSE)]
perf_test[, `:=` (q75=quantile(error, c(0.75)), 
		  q50=quantile(error, c(0.50)),
		  q25=quantile(error, c(0.25)),
		  mean_ga= mean(ga),
		  mean_ga_predicted = mean(ga_predicted),
		  median_error = median(error),
		  mean_error = mean(error)), by = formula]

test <- cmpt_test

cmpt_test_RF_XGB = copy(cmpt_test)
cmpt_test_RF_XGB[, `:=`(log_bpd=log(bpd),  log_ofd=log(ofd),   log_hp=log(hp),   log_ap=log(ap), log_fl=log(fl), 
                   sqrt_bpd=sqrt(bpd), sqrt_ofd=sqrt(ofd),  sqrt_hp=sqrt(hp), sqrt_ap =sqrt(ap),  sqrt_fl=sqrt(fl))]

# Prepare the test data
features_test <- vtreat::prepare(treatplan, subset(cmpt_test_RF_XGB, select = -c(enrid, ga_birth, ga)), varRestriction = new_vars) %>% as.matrix()
response_test <- subset(cmpt_test_RF_XGB, select = -c(enrid, ga_birth, ga))$logGA

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

test$hadlock <- predict_ga(test,"Hadlock")
test$intergrowth <- predict_ga(test,"Intergrowth")
test$randomforest <- exp(predictions(predict(ranger_model,subset(cmpt_test_RF_XGB, select = -c(enrid, ga_birth, ga)))))
test$xgboost <- exp(predict(xgb.fit.final, features_test))

hadlock <- subset(test, select = -c(intergrowth, randomforest, xgboost))
setnames(hadlock, "hadlock", "ga_predicted")
hadlock$eq  = 'ga = 10.85 + 0.060(hp)(fl) + 0.67(bpd) + 0.168(ap)'
hadlock$formula = "Hadlock"
hadlock$error <- hadlock$ga_predicted - hadlock$ga
hadlock$median_error = median(hadlock$error)
hadlock$AIC <- AIC(lm("log(ga) ~ hp * fl + bpd + ap", data = dat))
hadlock$BIC <- BIC(lm("log(ga) ~ hp * fl + bpd + ap", data = dat))
hadlock$RMSE = sqrt(mean((hadlock$ga - hadlock$ga_predicted)^2))
hadlock$R2 <- R2(predict(lm("log(ga) ~ hp * fl + bpd + ap", data = dat), data_test),hadlock$ga);
hadlock[, `:=` (
		  median_error = median(error),
		  mean_error = mean(error)), by = formula]
		  
		  
intergrowth <- subset(test,select = -c(hadlock, randomforest, xgboost))
setnames(intergrowth, "intergrowth", "ga_predicted")
intergrowth$eq  = 'log(ga) = 0.03243 × (log(hp))^2 + 0.001644 × fl × log(hp) + 3.813'
intergrowth$formula = "Intergrowth"
intergrowth$error <- intergrowth$ga_predicted - intergrowth$ga
intergrowth$median_error = median(intergrowth$error)
intergrowth$AIC <- AIC(lm("log(ga) ~ (log(hp))^2 + fl + log(hp)", data = dat))
intergrowth$BIC <- BIC(lm("log(ga) ~ (log(hp))^2 + fl + log(hp)", data = dat))
intergrowth$RMSE = sqrt(mean((intergrowth$ga - intergrowth$ga_predicted)^2))
intergrowth$R2 <- R2(predict(lm("log(ga) ~ hp * fl + bpd + ap", data = dat), data_test),intergrowth$ga);
intergrowth[, `:=` (
		  median_error = median(error),
		  mean_error = mean(error)), by = formula]


randomforest <- subset(test,select = -c(hadlock, intergrowth, xgboost))
setnames(randomforest, "randomforest", "ga_predicted")
randomforest$eq  = "RandomForestModel"
randomforest$formula = "RandomForest"
randomforest$error <- randomforest$ga_predicted - randomforest$ga
randomforest$median_error = median(randomforest$error)
randomforest$RMSE = RMSE(exp(predictions(predict(ranger_model,subset(cmpt_test_RF_XGB, select = -c(enrid, ga_birth, ga))))), cmpt_test_RF_XGB$ga) 
randomforest$R2 = cor(exp(predictions(predict(ranger_model,subset(cmpt_test_RF_XGB, select = -c(enrid, ga_birth, ga))))),cmpt_test_RF_XGB$ga)^2
randomforest[, `:=` (
		  median_error = median(error),
		  mean_error = mean(error)), by = formula]


xgboost <- subset(test,select = -c(hadlock, intergrowth, randomforest))
setnames(xgboost, "xgboost", "ga_predicted")
xgboost$eq  = "GradientBoostingModel"
xgboost$formula = "GradientBoosting"
xgboost$error <- xgboost$ga_predicted - xgboost$ga
xgboost$median_error = median(xgboost$error)
xgboost$RMSE = RMSE(exp(predict(xgb.fit.final, features_test)), cmpt_test_RF_XGB$ga) 
xgboost$R2 = cor(exp(predict(xgb.fit.final, features_test)),cmpt_test_RF_XGB$ga)^2
xgboost[, `:=` (
		  median_error = median(error),
		  mean_error = mean(error)), by = formula]
		  
kn <- rbindlist(list(hadlock, intergrowth, randomforest, xgboost), fill = TRUE)

perf_test = rbindlist(list(setDT(kn), perf_test), fill=TRUE)
perf_test2 <- unique(subset(perf_test, select = -c(enrid, ga, ga_predicted, error, ga_birth)))
perf_test2 <- perf_test2[order(abs(median_error))]
perf_test$formula1 <- factor(perf_test$formula, levels = rev(unique(perf_test2$formula)))

perf_test_df <- unique(subset(perf_test, select = c(eq, formula, RMSE, R2, AIC, BIC, mean_error, median_error)))
perf_test_df <- perf_test_df[order(abs(median_error))]

plot2 <- ggplot(perf_test, aes(y=error, x= formula1))+
  geom_violin(trim= FALSE, alpha = 0.5, fill = "red") +
  geom_boxplot(width = 0.2,  outlier.size = .5) +
  ggtitle("Garbhini Test-set")+
  xlab("error(in weeks)")+
  scale_y_continuous(breaks = -8:8, limits = c(-8,8)) + theme_bw() +
  labs(caption = paste("n = ",dim(data_test)[1])) +   xlab("") +
  geom_hline(aes(yintercept=0), linetype="dashed", size=0.05, colour="black") +
  theme(legend.position = "none") +
  coord_flip()


pdf(paste(outdir, "Garbhini_test_error_plot.pdf", sep = ""),   h=9, w=12)
plot2
dev.off()

write.table(subset(perf_test, select = -c(formula1)),
  paste(outdir, "Garbhini_test_error.tsv", sep = ""), 
  sep = "\t", row.names = F, quote = F)

write.table(perf_test_df, 
  paste(outdir, "Garbhini_test_error_table.tsv", sep = ""), 
  sep = "\t", row.names = F, quote = F)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

