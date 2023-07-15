# Garbhini-GA2
Script repository for the Garbhini-GA2 manuscript

### STEPS

#### 1. Remove outliers from the training data and generate the extended input used to train models 
```bash
Rscript 1a_data_preparation.R ./data/ ./results/ #outlier removal
Rscript 1b_generated_extended_input_data.R ./data/ ./data/extended_input_set/
```
#### 2. Randomforest and XGBoost hyperparameter tuning
```bash
Rscript 2a_hyperparameter_tuning_randomForest.R ./data/ ./results/
Rscript 2b_hyperparameter_tuning_xgboost.R ./data/ ./results/
```
#### 3. Run gapolyfitn in MATLAB
```Matlab
script1 = fullfile('./3a_gapolyfitn_CMDs_log_logGA.m');
run(script1)
script2 = fullfile('./3b_gapolyfitn_CMDs_classic_log_logGA.m');
run(script2)
script3 = fullfile('./3c_gapolyfitn_CMDs_classic_log_sqrt_logGA.m');
run(script3)
```
##### Process the results from gapolyfitn run
```bash
# pass the path to the directory containing results as an argument to the bash script
bash ./3d_gather_gapolyfitn_results.sh ./results/output_gapolyfitn_matlab_log
bash ./3d_gather_gapolyfitn_results.sh ./results/output_gapolyfitn_matlab_classic_log
bash ./3d_gather_gapolyfitn_results.sh ./results/output_gapolyfitn_matlab_classic_log_sqrt/

# gather results in one file
(cat ./results/gapolyfitn_formulas_consistent_1.txt; exec 0<./results/gapolyfitn_formulas_consistent_2.txt; read HEADER; cat; exec 0<./results/gapolyfitn_formulas_consistent_3.txt; read HEADER; cat) > ./results/gapolyfitn_formulas_consistent.txt
```

#### 4. Performance evaluation of all the models and selection of Garbhini-GA2
```bash
Rscript ./4_model_performance_evaluation.R ./results/gapolyfitn_formulas_consistent.txt ./data/train_23_dbscan.tsv ./data/test_23.tsv ./data/cmcv_validation.tsv ./data/RF_XGB_train_data.tsv ./results/figures/
```

#### 5. Compare Garbhini-GA2 against Hadlock and INTERGROWTH-21st formulae
```bash
Rscript ./5_garbhini_ga2_vs_published_formulas.R 
```
### DEPENDENCIES
[R version 4.3.1](https://cran.r-project.org/doc/manuals/r-patched/R-admin.html), R packages - [data.table](https://cran.r-project.org/web/packages/data.table/index.html), [optparse](https://cran.r-project.org/web/packages/optparse/index.html), [haven](https://cran.r-project.org/package=haven), 

