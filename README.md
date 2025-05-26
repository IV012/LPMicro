## LP-Micro: Analyzing Time-Varying Microbial Effects using Machine Learning

## Introduction

This GitHub repository includes the R code for `LP-Micro`, a machine learning-based method to predict the complex traits (e.g., early childhood caries, obesity.) using longitudinal metatranscriptomics. Additionally, `LP-Micro` can be used to analyze the temporal association between disease onset and microbiome by identifying the most predictive microbial taxa and ideal prognostic time.

## Contact and Citation

Yifan Dai: yfd@unc.edu, Di Wu: did@email.unc.edu

Dai Y, Qian Y, Qu Y, et al. LP-Micro Offers Interpretable Disease Outcome Prediction by Leveraging Microbial Biomarkers and Their Time-Varying Effects. *bioRxiv.* **2024**

## Installation

- Available for Linux and Mac.
- For Windows users, please install Rtools prior to install LP-Micro.

```r
install.packages("devtools")
devtools::install_github("IV012/LPMicro")
```

## Usage

`LP-Micro` includes the following components:

- **Microbial Taxa Screening `screen_taxa`:** LP-Micro combines group lasso and polynomial splines to screen disease-correlated taxa from longitudinal microbiome data.

- **Visitwise Prediction `visit_predict`:** LP-Micro supports popular machine learning methods for predicting the future disease onset using microbial profiles at **a single time point**, namely visitwise prediction in our paper. User can idenfity the most predictive time point of microbiome for disease outcomes.

- **Cumulative Prediction `cumulative_predict`:** LP-Micro also supports the prediction of disease outcomes using microbial **up to a chosen visit,** namely cumulative prediction. Cumulative prediction is usually more accurate than visitwise prediction with proper feature engineering using `screen_taxa`.

- **Feature Interpretation `taxa_interpret`:** LP-Micro provides statistical p-values to test whether the following features are useful for prediction, including (i) a microbial taxon, (ii) all information from one time point, and (iii) a microbial taxon at a chosen time point.

Feel free to check more details using `help` and apply `LP-Micro` to your data!

```r
help(screen_taxa)
help(cumulative_predict)
help(visit_predict)
help(taxa_interpret)
```

## Examples

In the following script, we introduce the general input and output of the taxa screening function.

```r
# general input format
n <- 200 # 200 individuals
p <- 500 # 500  microbial taxa
q <- 5 # 5 repeated measurements
x <- matrix(rnorm(n*p*q), n, p*q) # simulated microbial matrix
y <- rnorm(n) # simulated response matrix
idx <- rep(1:p, q) # identifier of the microbial taxa across the columns of x

# microbial feature screening, lambda is the tuning parameter
screen_taxa(x, y, idx, lambda = c(1e-2, 1e-4)) # this returns a list of taxa of the same length as lambda
```
Next, we introduce the machine learning based feature importance. First, we need to specify the machine
learning model and its parameter. Below is a template of the model specification. We take random forest
as an example.

```r
# models and parameters supported by LP-Micro
n_ensemble = 100 #100
n_perm = 100 #100

n_epoch = 1000 #1000
n_tree = 1000
node_size = 3
esCtrl = list(
  n.hidden = c(50, 40, 30, 20), 
  activate = "relu", 
  l1.reg = 10**-4, 
  early.stop.det = 1000, 
  n.batch = 50, 
  n.epoch = n_epoch, 
  learning.rate.adaptive = "adam", 
  plot = FALSE)

args_list = list(
  svm = list(method="svm", n.ensemble=n_ensemble),
  rf = list(
    method = "random_forest",
    n.ensemble=n_ensemble, 
    ntree=n_tree, 
    nodesize=node_size
  ),
  xgb = list(
    method="xgboost", 
    params=list(
      booster="gbtree", 
      objective="reg:squarederror", 
      eta=0.3, 
      gamma=0, 
      max_depth=5, 
      min_child_weight=1, 
      subsample=1, 
      colsample_bytree=1
    )
  ),
  lasso = list(method="lasso"),
  dnn = list(
    method = "ensemble_dnnet",
    n.ensemble = n_ensemble,
    esCtrl = esCtrl,
    verbose = 0
  )
)

mod_args <- args_list[[2]]
```

For cumulative/visit-wise prediction, we also need to specify the cumulative dataset. This
records the predictive accuracy up-to or of each visit.

```r
# taxa_list represents the microbial abundance measured up to a time point
taxa_list <- list()
for(t in 1:q){
  taxa_list[[t]] <- list(x=x[, 1:(p*t)], y=y)
}

cumulative_predict(taxa_list, mod_args, type="regression")

# taxa_list represents the microbial abundance measured at a time point
taxa_list <- list()
for(t in 1:q){
  taxa_list[[t]] <- list(x=x[, (p*(t-1)+1):(p*t)], y=y)
}

visit_predict(taxa_list, mod_args, type="regression")
```

Finally, we present the feature and visit importance from permutation importance test.

```r
idx_tx <- rep(1:p, q) # identifier of the microbial taxa across the columns of x
idx_tp <- c() # identifier of the measured time point
for(t in 1:q){
  idx_tp <- c(idx_tp, rep(t, p))
}

taxa_interpret(x, y, idx_tx, idx_tp, mod_args) # This gives the taxa and visit importance
```
