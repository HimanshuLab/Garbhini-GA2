#!/usr/bin/env Rscript 
args = commandArgs(trailingOnly=TRUE)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Veerendra P Gadekar 
# Email: ic36871@imail.iitm.ac.in; gpveerendra09@gmail.com

# script to collect the results from gapolyfitn in matlab and generate terms
# with variable names 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# library

library(data.table)
library(stringr)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# input arguments

terms <- args[1];
rmse  <- args[2];
variables <- args[3];
out1 <- args[4];
out2 <- args[5];

terms <- fread(terms, header = F); terms$V1 = basename(terms$V1)
rmse <- fread(rmse, header = F); rmse$V1 = basename(rmse$V1)
variables <- fread(variables, header = F); variables$V1 = basename(variables$V1)
variables <- setDF(unique(variables))
rmse <- setDF(unique(rmse))
var_rmse <- merge(variables, rmse, by = c("V1","V2"))
setnames(var_rmse, old = c("V3.x", "V3.y"), new = c("variables","RMSE"))

terms_var_rmse <- merge(terms, var_rmse, by = c("V1","V2"))

terms_var_rmse$V4 =
  sapply(1:dim(terms_var_rmse)[1],
     function(x){
        t = str_trim(strsplit(terms_var_rmse$V4[x], "[*]")[[1]], side = "both");
        v1 = strsplit(terms_var_rmse$variables[x], " ")[[1]];
        v2 = v1[!v1 %in% ""];
        v3 = sapply(sub(")^", "))^",
                sapply(1:length(v2),
                    function(i) gsub("x",paste("(",v2[i],")",sep=""), t[i]))),
                       function(j) ifelse(grepl("_", j), gsub("\\^", ")^", gsub("_","(",j)), j));
        f = paste(v3, collapse = " * ")
    })

terms_var_rmse$terms =
  sapply(1:dim(terms_var_rmse)[1],
     function(x){
        t = str_trim(strsplit(terms_var_rmse$V4[x], "[*]")[[1]], side = "both");
        v1 = strsplit(terms_var_rmse$variables[x], " ")[[1]];
        v2 = v1[!v1 %in% ""];
        v3 = sapply(sub(")^", "))^",
                sapply(1:length(v2),
                    function(i) gsub("x",paste("(",v2[i],")",sep=""), t[i]))),
                       function(j) ifelse(grepl("_", j), gsub("\\^", ")^", gsub("_","(",j)), j));
        o = ifelse(grepl("0", v3), "", v3);
        o2 = o[!o %in% ""];
        o3 = ifelse(grepl("1", o2), gsub("\\^1", "", o2), o2);
        f = paste(o3, collapse = " * ")
    })

terms_var_rmse = terms_var_rmse[terms_var_rmse[order(terms), .I ,by= c("V1","V2")]$I]
terms_var_rmse = subset(terms_var_rmse, terms != "")
terms_var_rmse[, forms := paste("I(", terms,")", sep = "")]
terms_var_rmse[, formula := paste("log(ga) ~ ", paste(forms, collapse = " + "), sep = ""), by = c("V1","V2")]
terms_var_rmse[, formula := ifelse(grepl("I\\(\\) +", formula), gsub("I\\(\\) \\+", "", formula), formula), 1:nrow(terms_var_rmse)]
terms_var_rmse_main = unique(subset(terms_var_rmse, select = c(V1, V2, formula, RMSE)))
terms_var_rmse_main[, global_maxima := uniqueN(V2), by = c("V1","formula")]
setnames(terms_var_rmse_main, old = c("V1","V2"), new = c("input_data","iteration"))
terms_var_rmse_out = unique(subset(subset(terms_var_rmse_main, global_maxima >5), select = -c(iteration)))

write.table(terms_var_rmse_main, out1, sep = "\t", quote = F, row.names = F)
write.table(terms_var_rmse_out, out2, sep = "\t", quote = F, row.names = F)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
