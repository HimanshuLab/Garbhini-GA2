# Garbhini-GA2
Script repository for the Garbhini-GA2 manuscript

#### 1. Remove outliers from the training data and generate the extended input used to train models 
```{r, engine = 'bash', eval = FALSE}
Rscript 1a_data_preparation.R ./data/ ./results/ #outlier removal
Rscript 1b_generated_extended_input_data.R ./data/ ./data/extended_input_set/ #generate data used for training
```
#### 2. Randomforest and XGBoost hyperparameter tuning
```{r, engine = 'bash', eval = FALSE}
Rscript 2a_hyperparameter_tuning_randomForest.R ./data/ ./results/ 
```
