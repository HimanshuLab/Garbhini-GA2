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
[R version 4.3.1](https://cran.r-project.org/doc/manuals/r-patched/R-admin.html), R packages - [data.table](https://cran.r-project.org/web/packages/data.table/index.html), [tidyverse](https://cran.r-project.org/package=tidyverse), [dbscan](https://cran.r-project.org/package=dbscan), [readxl](https://cran.r-project.org/package=readxl), [StatMatch](https://cran.r-project.org/package=StatMatch), [xlsx](https://cran.r-project.org/package=xlsx), [openxlsx](https://cran.r-project.org/package=openxlsx), [ranger](https://cran.r-project.org/package=ranger), [randomForest](https://cran.r-project.org/package=randomForest), [gridExtra](https://cran.r-project.org/package=gridExtra), [xgboost](https://cran.r-project.org/package=xgboost), [rsample](https://cran.r-project.org/package=rsample), [caret](https://cran.r-project.org/package=caret), [DescTools](https://cran.r-project.org/web/packages/DescTools/index.html), [ggpubr](https://cran.r-project.org/package=ggpubr), [ggstatsplot](https://cran.r-project.org/package=ggstatsplot), [dlookr](https://cran.r-project.org/package=dlookr), [car](https://cran.r-project.org/package=car), [repr](https://cran.r-project.org/package=repr), [BlandAltmanLeh](https://cran.r-project.org/package=BlandAltmanLeh), [optparse](https://cran.r-project.org/web/packages/optparse/index.html), [haven](https://cran.r-project.org/package=haven)

[MATLAB R2022b](https://www.mathworks.com/products/new_products/r2022b-transition.html), [gapolyfitn](https://www.mathworks.com/matlabcentral/fileexchange/25499-gapolyfitn)

