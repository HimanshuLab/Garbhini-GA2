# Garbhini-GA2
Script repository for the Garbhini-GA2 manuscript

#### Scripts used to remove outliers from the training data and generate the extended input used to train models 
```{r, engine = 'bash', eval = FALSE}
Rscript 1a_data_preparation.R ./data/ ./results/ #outlier removal
Rscript 1b_generated_extended_input_data.R ./data/ ./data/extended_input_set/ #generate extended data used to train
```
#### 2. Script used to generate the extended input data set used to train models 
```{r, engine = 'bash', eval = FALSE}
```
