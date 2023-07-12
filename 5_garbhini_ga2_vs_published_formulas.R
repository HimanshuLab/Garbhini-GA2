#!/usr/bin/env Rscript 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Veerendra P Gadekar    
# Email: ic36871@imail.iitm.ac.in; gpveerendra09@gmail.com

# script to generate the error distribution plots for the GA formulas 
# along with the intergrowth and hadlock formulas

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# libraries
library(data.table)
library(tidyverse)
library(gridExtra)
library(DescTools)
library(ggpubr)
library(ggstatsplot)
library(dlookr)
library(car)
library(repr)
library(BlandAltmanLeh)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

cmcv_predictions <- "./results/figures/CMCV_validation_error.tsv"
garbhini_predictions <- "./results/figures/Garbhini_test_error.tsv"

outdir <- "./results/figures/";
outdir2 <- "./results/";

cmcv_predictions <- fread(cmcv_predictions)
garbhini_predictions <- unique(fread(garbhini_predictions))

cmcv_data <- subset(cmcv_predictions, formula %in% c("log(ga) ~ I((log(bpd)) * (log(ofd))) + I((log(hp))^2)", "Hadlock",  "Intergrowth"))
cmcv_data$formula2 = ifelse(cmcv_data$formula == "log(ga) ~ I((log(bpd)) * (log(ofd))) + I((log(hp))^2)", "Garbhini-GA2", cmcv_data$formula)
cmcv_data$formula2 = factor(cmcv_data$formula2, level = c("Garbhini-GA2", "Intergrowth", "Hadlock"))

garbhini_data <- unique(subset(garbhini_predictions, formula %in% c("log(ga) ~ I((log(bpd)) * (log(ofd))) + I((log(hp))^2)", "Hadlock",  "Intergrowth")))
garbhini_data$formula2 = ifelse(garbhini_data$formula == "log(ga) ~ I((log(bpd)) * (log(ofd))) + I((log(hp))^2)", "Garbhini-GA2", garbhini_data$formula)
garbhini_data$formula2 = factor(garbhini_data$formula2, level = c("Garbhini-GA2", "Intergrowth", "Hadlock"))
#garbhini_data = subset(garbhini_data, error < 10)

# Check for Normality (Shapiro-Wilk test of normality)
cmcv_data %>% group_by(formula) %>% normality(error)
garbhini_data %>% group_by(formula) %>% normality(error)

# Check for Homogeneity of variance (Levene's test for homogeneity of variance)
leveneTest(error ~ formula, cmcv_data)
leveneTest(error ~ formula, garbhini_data)

MedianCI(subset(cmcv_data, formula2 == "Hadlock")$error)
MedianCI(subset(cmcv_data, formula2 == "Intergrowth")$error)
MedianCI(subset(cmcv_data, formula2 == "Garbhini-GA2")$error)

MedianCI(subset(garbhini_data, formula2 == "Hadlock")$error)
MedianCI(subset(garbhini_data, formula2 == "Intergrowth")$error)
MedianCI(subset(garbhini_data, formula2 == "Garbhini-GA2")$error)


pdf(paste(outdir, "error_plot_CMCV_Validation.pdf", sep=""), h=8, w=7)
ggbetweenstats(
  data = cmcv_data, 
  x = formula2, 
  y = error,  
  type= "nonparametric",
  title = "VALIDATION set",
  xlab = "GA estimation model", 
  ylab = paste("Error in GA estimation", sep = ""),
  var.equal = FALSE,
  p.adjust.method = "bonferroni", 
  centrality.label.args = list(size  = 3),
  ggsignif.args = list(textsize = 4, tip_length = 0.01),  
  ggplot.component = list(theme(text = element_text(size = 15)))) + 
  scale_y_continuous(limits = c(-10, 10))
dev.off()


pdf(paste(outdir, "error_plot_Garbhini_Test.pdf", sep=""), h=8, w=7)
ggbetweenstats(
  data = garbhini_data, 
  x = formula2, 
  y = error,  
  type= "nonparametric",
  title = "Test set",
  xlab = "GA estimation model", 
  ylab = paste("Error in GA estimation", sep = ""),
  var.equal = FALSE,
  p.adjust.method = "bonferroni", 
  centrality.label.args = list(size  = 3),
  ggsignif.args = list(textsize = 4, tip_length = 0.01), 
  ggplot.component = list(theme(text = element_text(size = 15)))) #+ 
#  scale_y_continuous(limits = c(-10, 10))
dev.off()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
# absolute error plots
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

# cmcv
my_comparisons <- list( c( "log(ga) ~ I((log(bpd)) * (log(ofd))) + I((log(hp))^2)", "Hadlock"), c("log(ga) ~ I((log(bpd)) * (log(ofd))) + I((log(hp))^2)", "Intergrowth"))
pdf(paste(outdir, "error_plot_CMCV_Validation_updated.pdf", sep=""), h=8, w=12)
ggplot(subset(cmcv_predictions, formula %in% c("Hadlock","Intergrowth","log(ga) ~ I((log(bpd)) * (log(ofd))) + I((log(hp))^2)")), aes(x = formula, y =  error, fill = formula)) +
  geom_violin() + ggtitle("Distribution of errors on CMCV Validation-set ") +
  geom_boxplot(width=0.05, outlier.size = .5) + 
  scale_y_continuous(breaks = -8:8, limits = c(-8,8)) + theme_bw() +
  labs(caption = paste("n = ",uniqueN(cmcv_predictions$enrid))) +
  geom_hline(aes(yintercept=0), linetype="dashed", size=0.2, colour="black") +
  stat_summary(fun.y="mean", geom="point", size=2, color="black", shape = 8) +
  stat_compare_means(comparisons = my_comparisons, label.y = c(7, 6), method = "t.test") + coord_flip()

dev.off()


# garbhini
pdf(paste(outdir, "error_plot_Garbhini_Validation.pdf", sep=""), h=8, w=12)
ggplot(subset(garbhini_predictions, formula %in% c("Hadlock","Intergrowth","log(ga) ~ I((log(bpd)) * (log(ofd))) + I((log(hp))^2)")), aes(x = formula, y = error, fill = formula)) +
  geom_violin() + ggtitle("Distribution of errors on Garbhini Test-set ") +
  geom_boxplot(width=0.05, outlier.size = .5) + 
  scale_y_continuous(breaks = -8:8, limits = c(-8,8)) + theme_bw() +
  labs(caption = paste("n = ",length(subset(garbhini_predictions, formula == "Hadlock")$enrid))) +
  geom_hline(aes(yintercept=0), linetype="dashed", size=0.2, colour="black") + 
#  geom_vline(aes(xintercept=0), linetype="dashed", size=0.2, colour="black")# + 
  stat_summary(fun.y="mean", geom="point", size=2, color="black", shape = 8) +
  stat_compare_means(comparisons = my_comparisons, label.y = c(7, 6), method = "t.test") + coord_flip()
dev.off()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
# parity plots  
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

# cmcv
pdf(paste(outdir, "parity_plot_CMCV_Validation_hadlock_vs_gold.pdf", sep=""), h=8, w=8)
ggplot(subset(cmcv_predictions, formula %in% c("Hadlock")), aes(x = ga,y = ga_predicted)) + 
  geom_point(alpha = 0.8, shape = 1) +
  geom_abline(slope = 1,color = "red") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_bw() +
  ggtitle("CMCV Validation-set: Hadlock predicted GA vs gold standard") +
  ylab("Predicted GA (weeks)") +
  xlab("Gold Standard GA (weeks)")   
dev.off()

cmcv_hadlock_vs_ga =
subset(cmcv_predictions, formula %in% c("Hadlock"))[, .(enrid, diff_hadlock_ga = abs(ga_predicted - ga), ga_predicted, ga, t_weeks)]
cmcv_hadlock_vs_ga_df = subset(cmcv_hadlock_vs_ga, diff_hadlock_ga > 3)

write.table(cmcv_hadlock_vs_ga_df,
  paste(outdir2, "gold_standard_vs_Hadlock_cmcv_diff.tsv", sep = ""),
  sep = "\t", row.names = F, quote = F)

pdf(paste(outdir, "parity_plot_CMCV_Validation_intergrowth_vs_gold.pdf", sep=""), h=8, w=8)
ggplot(subset(cmcv_predictions, formula %in% c("Intergrowth")), aes(x = ga,y = ga_predicted)) + 
  geom_point(alpha = 0.8, shape = 1) +
  geom_abline(slope = 1,color = "red") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_bw() +
  ggtitle("CMCV Validation-set: Intergrowth predicted GA vs gold standard") +
  ylab("Predicted GA (weeks)") +
  xlab("Gold Standard GA (weeks)")   
dev.off()

cmcv_intergrowth_vs_ga =
subset(cmcv_predictions, formula %in% c("Intergrowth"))[, .(enrid, diff_intergrowth_ga = abs(ga_predicted - ga), ga_predicted, ga, t_weeks)]
cmcv_intergrowth_vs_ga_df = subset(cmcv_intergrowth_vs_ga, diff_intergrowth_ga > 3)

write.table(cmcv_intergrowth_vs_ga_df,
  paste(outdir2, "gold_standard_vs_Intergrowth_cmcv_diff.tsv", sep = ""),
  sep = "\t", row.names = F, quote = F)

pdf(paste(outdir, "parity_plot_CMCV_Validation_garbhiniga2_vs_gold.pdf", sep=""), h=8, w=8)
ggplot(subset(cmcv_predictions, formula %in% c("log(ga) ~ I((log(bpd)) * (log(ofd))) + I((log(hp))^2)")), aes(x = ga,y = ga_predicted)) + 
  geom_point(alpha = 0.8, shape = 1) +
  geom_abline(slope = 1,color = "red") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_bw() +
  ggtitle("CMCV Validation-set: GarbhiniGA2 predicted GA vs gold standard") +
  ylab("Predicted GA (weeks)") +
  xlab("Gold Standard GA (weeks)")   
dev.off()

cmcv_garbhiniga2_vs_ga =
subset(cmcv_predictions, formula %in% c("log(ga) ~ I((log(bpd)) * (log(ofd))) + I((log(hp))^2)"))[, .(enrid, diff_garbhiniga2_ga = abs(ga_predicted - ga), ga_predicted, ga, t_weeks)]
cmcv_garbhiniga2_vs_ga_df = subset(cmcv_garbhiniga2_vs_ga, diff_garbhiniga2_ga > 3)

write.table(cmcv_garbhiniga2_vs_ga_df,
  paste(outdir2, "gold_standard_vs_GarbhiniGA2_cmcv_diff.tsv", sep = ""),
  sep = "\t", row.names = F, quote = F)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# garbhini
pdf(paste(outdir, "parity_plot_Garbhini_Validation_hadlock_vs_gold.pdf", sep=""), h=8, w=8)
ggplot(subset(garbhini_predictions, formula %in% c("Hadlock")), aes(x = ga,y = ga_predicted)) + 
  geom_point(alpha = 0.8, shape = 1) +
  geom_abline(slope = 1,color = "red") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_bw() +
  ggtitle("Garnhini Test-set: Hadlock predicted GA vs gold standard") +
  ylab("Predicted GA (weeks)") +
  xlab("Gold Standard GA (weeks)")   
dev.off()


garbhini_hadlock_vs_ga =
subset(garbhini_predictions, formula %in% c("Hadlock"))[, .(enrid, diff_hadlock_ga = abs(ga_predicted - ga), ga_predicted, ga, ga_birth)]
garbhini_hadlock_vs_ga_df = subset(garbhini_hadlock_vs_ga, diff_hadlock_ga > 3)

write.table(garbhini_hadlock_vs_ga_df,
  paste(outdir2, "gold_standard_vs_Hadlock_Garbhini_diff.tsv", sep = ""),
  sep = "\t", row.names = F, quote = F)

pdf(paste(outdir, "parity_plot_Garbhini_Validation_intergrowth_vs_gold.pdf", sep=""), h=8, w=8)
ggplot(subset(garbhini_predictions, formula %in% c("Intergrowth")), aes(x = ga,y = ga_predicted)) + 
  geom_point(alpha = 0.8, shape = 1) +
  geom_abline(slope = 1,color = "red") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_bw() +
  ggtitle("Garnhini Test-set: Intergrowth predicted GA GA vs gold standard") +
  ylab("Predicted GA (weeks)") +
  xlab("Gold Standard GA (weeks)")   
dev.off()

garbhini_intergrowth_vs_ga =
subset(garbhini_predictions, formula %in% c("Intergrowth"))[, .(enrid, diff_intergrowth_ga = abs(ga_predicted - ga), ga_predicted, ga, ga_birth)]
garbhini_intergrowth_vs_ga_df = subset(garbhini_intergrowth_vs_ga, diff_intergrowth_ga > 3)

write.table(garbhini_intergrowth_vs_ga_df,
  paste(outdir2, "gold_standard_vs_Intergrowth_Garbhini_diff.tsv", sep = ""),
  sep = "\t", row.names = F, quote = F)

pdf(paste(outdir, "parity_plot_Garbhini_Validation_garbhiniga2_vs_gold.pdf", sep=""), h=8, w=8)
ggplot(subset(garbhini_predictions, formula %in% c("log(ga) ~ I((log(bpd)) * (log(ofd))) + I((log(hp))^2)")), aes(x = ga,y = ga_predicted)) + 
  geom_point(alpha = 0.8, shape = 1) +
  geom_abline(slope = 1,color = "red") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_bw() +
  ggtitle("Garnhini Test-set: GarbhiniGA2 predicted GA vs gold standard") +
  ylab("Predicted GA (weeks)") +
  xlab("Gold Standard GA (weeks)")   
dev.off()

garbhini_garbhiniga2_vs_ga =
unique(subset(garbhini_predictions, formula %in% c("log(ga) ~ I((log(bpd)) * (log(ofd))) + I((log(hp))^2)"))[, .(enrid, diff_garbhiniga2_ga = abs(ga_predicted - ga), ga_predicted, ga)])
garbhini_garbhiniga2_vs_ga_df = subset(garbhini_garbhiniga2_vs_ga, diff_garbhiniga2_ga > 3)

write.table(garbhini_garbhiniga2_vs_ga_df,
  paste(outdir2, "gold_standard_vs_Garbhiniga2_Garbhini_diff.tsv", sep = ""),
  sep = "\t", row.names = F, quote = F)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
# BLAND-ALTMAN ANALYSIS 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

tcmcv_hadlock_vs_ga = subset(cmcv_hadlock_vs_ga, select  = c(ga, ga_predicted));
tcmcv_intergrowth_vs_ga = subset(cmcv_intergrowth_vs_ga, select  = c(ga_predicted));
tcmcv_garbhiniga2_vs_ga = subset(cmcv_garbhiniga2_vs_ga, select  = c(ga_predicted)); 

setDF(setnames(tcmcv_hadlock_vs_ga, old = "ga_predicted", new = "ga_predicted_hadlock"))
setDF(setnames(tcmcv_intergrowth_vs_ga, old = "ga_predicted", new = "ga_predicted_intergrowth"))
setDF(setnames(tcmcv_garbhiniga2_vs_ga, old = "ga_predicted", new = "ga_predicted_garnhiniga2"))

cmcv_preterm_test <- cbind(tcmcv_hadlock_vs_ga, tcmcv_intergrowth_vs_ga, tcmcv_garbhiniga2_vs_ga)
cmcv_pairwise_BA_test <- 
sapply(cmcv_preterm_test, function(i){
  sapply(cmcv_preterm_test,function(j){
    BA_stats <- bland.altman.stats(i,j, conf.int = 0.95)
    paste0(BA_stats$mean.diffs %>% round(3),
           "( ",BA_stats$lower.limit %>% round(3),
           ", ",BA_stats$upper.limit %>% round(3),")")
    
  })
})  %>% (function(x) {
  x[upper.tri(x)] <- NA
  x
}) 


tgarbhini_hadlock_vs_ga = subset(garbhini_hadlock_vs_ga, select  = c(ga, ga_birth, ga_predicted));
tgarbhini_intergrowth_vs_ga = subset(garbhini_intergrowth_vs_ga, select  = c(ga_predicted));
tgarbhini_garbhiniga2_vs_ga = subset(garbhini_garbhiniga2_vs_ga, select  = c(ga_predicted)); 

setDF(setnames(tgarbhini_hadlock_vs_ga, old = "ga_predicted", new = "ga_predicted_hadlock"))
setDF(setnames(tgarbhini_intergrowth_vs_ga, old = "ga_predicted", new = "ga_predicted_intergrowth"))
setDF(setnames(tgarbhini_garbhiniga2_vs_ga, old = "ga_predicted", new = "ga_predicted_garbhiniga2"))

garbhini_preterm_test <- cbind(tgarbhini_hadlock_vs_ga, tgarbhini_intergrowth_vs_ga, tgarbhini_garbhiniga2_vs_ga)
garbhini_pairwise_BA_test <- 
sapply(garbhini_preterm_test, function(i){
  sapply(garbhini_preterm_test,function(j){
    BA_stats <- bland.altman.stats(i,j, conf.int = 0.95)
    paste0(BA_stats$mean.diffs %>% round(3),
           "( ",BA_stats$lower.limit %>% round(3),
           ", ",BA_stats$upper.limit %>% round(3),")")
    
  })
}) %>% (function(x) {
  x[lower.tri(x)] <- NA
  x
})

bland_altman_garbhini_cmcv <- matrix(NA, nrow = 4, ncol = 4)
bland_altman_garbhini_cmcv[upper.tri(bland_altman_garbhini_cmcv)] <- garbhini_pairwise_BA_test[upper.tri(garbhini_pairwise_BA_test)]
bland_altman_garbhini_cmcv[lower.tri(bland_altman_garbhini_cmcv)] <- cmcv_pairwise_BA_test[lower.tri(cmcv_pairwise_BA_test)]
dimnames(bland_altman_garbhini_cmcv) = list(c("Gold Standard" ,"Hadlock", "Intergrowth", "Garbhini-GA2"), c("Gold Standard" ,"Hadlock", "Intergrowth", "Garbhini-GA2"))
pdf(paste(outdir, "Bland_Altman.pdf", sep=""), h=10, w=10)
grid.table(bland_altman_garbhini_cmcv)
dev.off()

                  
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
# ptb rate plots
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

tgarbhini_hadlock_vs_ga = subset(garbhini_hadlock_vs_ga, select  = c(enrid, ga, ga_birth, ga_predicted));
tgarbhini_intergrowth_vs_ga = subset(garbhini_intergrowth_vs_ga, select  = c(ga_predicted));
tgarbhini_garbhiniga2_vs_ga = subset(garbhini_garbhiniga2_vs_ga, select  = c(ga_predicted)); 

setDF(setnames(tgarbhini_hadlock_vs_ga, old = "ga_predicted", new = "ga_predicted_hadlock"))
setDF(setnames(tgarbhini_intergrowth_vs_ga, old = "ga_predicted", new = "ga_predicted_intergrowth"))
setDF(setnames(tgarbhini_garbhiniga2_vs_ga, old = "ga_predicted", new = "ga_predicted_garbhiniga2"))

garbhini_preterm_test <- cbind(tgarbhini_hadlock_vs_ga, tgarbhini_intergrowth_vs_ga, tgarbhini_garbhiniga2_vs_ga)

#@@@@@@@@@@@@@@@@@@@@@@@@

# compute GA birth the original data
# ga1 gestional age assessed by GA1 model +  Date of delivery - Date of CRL USG scan 
# To each of the variables (GARBHINI-GA2, Hadlock, Intergrowth & Gold-standard GA), 
# we add (date of delivery - date of US scan) to derive GA at delivery using each of the dating model

full_train_data <- paste('/home/veer/scratch/IBSE/GarbhIni2/Garbhini-GA2-main/development/data/', 'full_train.dta', sep = "");
raw_data<- fread("/home/veer/scratch/IBSE/GarbhIni2/Garbhini-GA2-main/data/datasets/merged_final_dataset.csv")
raw_data_mod <- subset(raw_data, select = c(enrid, vdt_as1, vdt_fwb1, vdt_fwbv1, op_dt_event))
garbhini_preterm_test_df = merge(garbhini_preterm_test, raw_data_mod, by = "enrid", all.x = TRUE)
setDT(garbhini_preterm_test_df)
garbhini_preterm_test_df[order(ga),index:= 1:.N, by=enrid]
garbhini_preterm_test_df[, `:=`(vdt_as1_diff = (op_dt_event - vdt_as1)/7, vdt_fwb1_diff = (op_dt_event - vdt_fwb1)/7, vdt_fwbv1_diff = (op_dt_event - vdt_fwbv1)/7)]

garbhini_preterm_test_df[is.na(vdt_as1) & index == 1 & is.na(vdt_fwb1) & !is.na(vdt_fwbv1), index := 3]
garbhini_preterm_test_df[!is.na(vdt_as1) & index == 2 & is.na(vdt_fwb1) & !is.na(vdt_fwbv1), index := 3]
garbhini_preterm_test_df[is.na(vdt_as1) & index == 1 & !is.na(vdt_fwb1), index := 2]

garbhini_preterm_test_df[index == 1, `:=`(ga_mod = ga+vdt_as1_diff, ga_predicted_hadlock_mod = ga_predicted_hadlock+vdt_as1_diff, ga_predicted_intergrowth_mod = ga_predicted_intergrowth+vdt_as1_diff, ga_predicted_garbhiniga2_mod = ga_predicted_garbhiniga2+vdt_as1_diff)]
garbhini_preterm_test_df[index == 2, `:=`(ga_mod = ga+vdt_fwb1_diff, ga_predicted_hadlock_mod = ga_predicted_hadlock+vdt_fwb1_diff, ga_predicted_intergrowth_mod = ga_predicted_intergrowth+vdt_fwb1_diff, ga_predicted_garbhiniga2_mod = ga_predicted_garbhiniga2+vdt_fwb1_diff)]
garbhini_preterm_test_df[index == 3, `:=`(ga_mod = ga+vdt_fwbv1_diff, ga_predicted_hadlock_mod = ga_predicted_hadlock+vdt_fwbv1_diff, ga_predicted_intergrowth_mod = ga_predicted_intergrowth+vdt_fwbv1_diff, ga_predicted_garbhiniga2_mod = ga_predicted_garbhiniga2+vdt_fwbv1_diff)]

preterm_test <- subset( garbhini_preterm_test_df, select = c(ga_predicted_hadlock_mod, ga_predicted_intergrowth_mod, ga_predicted_garbhiniga2_mod, ga_mod))

# 95% Confidence intervals with percentage preterm births
conf_test <- preterm_test %>% sapply(function(i) BinomCI(sum(i < 37),length(i), conf.level = 0.95, method = "clopper-pearson")) %>% round(3)
apply(preterm_test,2,function(i) sum(i < 37)*100/length(i)) %>% range()

sapply(colnames(setDF(preterm_test)), function(i){
  sapply(colnames(preterm_test), function(j){
    (c((preterm_test[,i] < 37) %>% table , (preterm_test[,j] < 37) %>% table) %>%
       matrix(nrow = 2) %>% 
       fisher.test())$p.value
  })
}) %>% (function(x) {
  x[lower.tri(x)] <- NA
  x
}) %>%  melt() %>% na.omit() %>% mutate( value = p.adjust(value, method = "bonferroni")) %>% filter(value < 0.05)
colnames(preterm_test) <- c("Hadlock","Intergrowth", "log(ga) ~ I((log(bpd)) * (log(ofd))) + I((log(hp))^2)","Gold Standard")

table_preterm_test <- preterm_test %>% sapply(function(i) BinomCI(sum(i < 37),length(i), conf.level = 0.95, method = "clopper-pearson")) %>% 
  t() %>% `*`(100) %>% round(2) %>%
  {
    data.frame(formula = row.names(.),
               PTB_rate = .[,1],
               CI_lower = .[,2],
               CI_higher = .[,3],
               row.names = NULL)
  } 

ptbrate_plot <- table_preterm_test %>%
  mutate(formula = fct_relevel(formula, "Hadlock","Intergrowth","log(ga) ~ I((log(bpd)) * (log(ofd))) + I((log(hp))^2)")) %>%
  ggplot(aes(x = formula, y = PTB_rate)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_higher), width = .2, size = 1) +
  coord_flip() +
  labs(title = 'PTB rates calculated by models on GarbhIni Test-set', y = 'PTB rate', x = 'Model') +
  theme_bw(base_size = 13) 

pdf(paste(outdir, "Estimated_ptb_rate.pdf", sep=""), h=5, w=8)
ptbrate_plot
dev.off()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
# Jaccard Similarity Calculation
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

garbhini_preterm_class <- preterm_test
garbhini_preterm_class[garbhini_preterm_class < 37] <- 1
garbhini_preterm_class[!garbhini_preterm_class < 37] <- 0
garbhini_preterm_vals <- sapply(colnames(preterm_test), function(i){
  sapply(colnames(preterm_test), function(j){
    round((nrow(garbhini_preterm_class[ which(garbhini_preterm_class[,i] == 1 & garbhini_preterm_class[,j] == 1) , ])*100)/nrow(garbhini_preterm_class[ which(garbhini_preterm_class[,i] == 1 | garbhini_preterm_class[,j] == 1) , ]),3)
  })##Calculating Jaccard index for all columns
}) %>% data.frame()

garbhini_preterm_class_mat <- as.matrix(garbhini_preterm_vals)
dimnames(garbhini_preterm_class_mat) <- list(colnames(preterm_test),colnames(preterm_test))##Joining Jaccard and BA Analysis results

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
# performance against the Gold Standard
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

preterm_plot_h <- preterm_test[c("Hadlock", "Gold Standard")]
preterm_plot_h[preterm_plot_h$Hadlock>=37 & preterm_plot_h$`Gold Standard`>=37,'level'] <- 'Term'
preterm_plot_h[preterm_plot_h$Hadlock>=37 & preterm_plot_h$`Gold Standard`<37,'level'] <- 'Preterm, GA'
preterm_plot_h[preterm_plot_h$Hadlock<37 & preterm_plot_h$`Gold Standard`>=37,'level'] <- 'Preterm, Hadlock'
preterm_plot_h[preterm_plot_h$Hadlock<37 & preterm_plot_h$`Gold Standard`<37,'level'] <- 'Preterm'

plot_preterm_h <- ggplot(preterm_plot_h) +
  geom_point(aes(x=Hadlock, y= `Gold Standard`, color=level),size=1, alpha=0.5) +
  geom_hline(yintercept=37, linetype="dashed", size=0.5) +
  geom_vline(xintercept=37, linetype="dashed", size=0.5) +
  scale_colour_manual(values=c("Term" = "darkgreen", "Preterm" = "red", "Preterm, GA" = "blue", "Preterm, Hadlock" = "orange"))+
  labs(x="GA by Hadlock at birth (weeks)", y="True GA (weeks)", 
       color="PTB by formula", caption=paste("n = ",dim(garbhini_hadlock_vs_ga)[1],sep=""))+
  ggtitle('True GA vs Hadlock')+
  theme_bw()+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_preterm_h

pdf(paste(outdir, "gold_vs_hadlock.pdf", sep=""), h=6, w=7)
plot_preterm_h
dev.off()


preterm_plot_i <- preterm_test[c("Intergrowth", "Gold Standard")]
preterm_plot_i[preterm_plot_i$Intergrowth>=37 & preterm_plot_i$`Gold Standard`>=37,'level'] <- 'Term'
preterm_plot_i[preterm_plot_i$Intergrowth>=37 & preterm_plot_i$`Gold Standard`<37,'level'] <- 'Preterm, GA'
preterm_plot_i[preterm_plot_i$Intergrowth<37 & preterm_plot_i$`Gold Standard`>=37,'level'] <- 'Preterm, Intergrowth'
preterm_plot_i[preterm_plot_i$Intergrowth<37 & preterm_plot_i$`Gold Standard`<37,'level'] <- 'Preterm'

plot_preterm_i <- ggplot(preterm_plot_i) +
  geom_point(aes(x=Intergrowth, y= `Gold Standard`, color=level),size=1, alpha=0.5) +
  geom_hline(yintercept=37, linetype="dashed", size=0.5) +
  geom_vline(xintercept=37, linetype="dashed", size=0.5) +
  scale_colour_manual(values=c("Term" = "darkgreen", "Preterm" = "red", "Preterm, GA" = "blue", "Preterm, Intergrowth" = "orange"))+
  labs(x="GA by Intergrowth at birth (weeks)", y="True GA (weeks)", 
       color="PTB by formula", caption=paste("n = ",dim(preterm_test)[1],sep=""))+
  ggtitle('True GA vs Intergrowth')+
  theme_bw()+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_preterm_i

pdf(paste(outdir, "gold_vs_intergrowth.pdf", sep=""), h=6, w=7)
plot_preterm_i
dev.off()


preterm_plot_g <- preterm_test[c("log(ga) ~ I((log(bpd)) * (log(ofd))) + I((log(hp))^2)", "Gold Standard")]
preterm_plot_g[preterm_plot_g$`log(ga) ~ I((log(bpd)) * (log(ofd))) + I((log(hp))^2)`>=37 & preterm_plot_g$`Gold Standard`>=37,'level'] <- 'Term'
preterm_plot_g[preterm_plot_g$`log(ga) ~ I((log(bpd)) * (log(ofd))) + I((log(hp))^2)`>=37 & preterm_plot_g$`Gold Standard`<37,'level'] <- 'Preterm, GA'
preterm_plot_g[preterm_plot_g$`log(ga) ~ I((log(bpd)) * (log(ofd))) + I((log(hp))^2)`<37 & preterm_plot_g$`Gold Standard`>=37,'level'] <- 'Preterm, Garbhini-GA2'
preterm_plot_g[preterm_plot_g$`log(ga) ~ I((log(bpd)) * (log(ofd))) + I((log(hp))^2)`<37 & preterm_plot_g$`Gold Standard`<37,'level'] <- 'Preterm'


plot_preterm_g <- ggplot(preterm_plot_g) +
  geom_point(aes(x=`log(ga) ~ I((log(bpd)) * (log(ofd))) + I((log(hp))^2)`, y= `Gold Standard`, color=level),size=1, alpha=0.5) +
  geom_hline(yintercept=37, linetype="dashed", size=0.5) +
  geom_vline(xintercept=37, linetype="dashed", size=0.5) +
  scale_colour_manual(values=c("Term" = "darkgreen", "Preterm" = "red", "Preterm, GA" = "blue", "Preterm, Garbhini-GA2" = "orange"))+
  labs(x="GA by Garbhini-GA2 at birth (weeks)", y="True GA (weeks)", 
       color="PTB by formula", caption=paste("n = ",dim(preterm_test)[1],sep=""))+
  ggtitle('True GA vs Garbhini-GA2')+
  theme_bw()+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_preterm_g

pdf(paste(outdir, "gold_vs_garbhini_ga2.pdf", sep=""), h=6, w=7)
plot_preterm_g
dev.off()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

