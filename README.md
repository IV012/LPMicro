## LP-Micro: Analyzing Time-Varying Microbial Effects using Machine Learning

## Introduction

This GitHub repository includes the R code for `LP-Micro`, a machine learning-based method to predict the complex traits (e.g., early childhood caries, obesity.) using longitudinal metatranscriptomics. Additionally, `LP-Micro` can be used to analyze the temporal association between disease onset and microbiome by identifying the most predictive microbial taxa and ideal prognostic time.

## Contact and Citation

Yifan Dai: yfd@unc.edu, Di Wu: did@email.unc.edu

Dai Y, Qian Y, Qu Y, et al. LP-Micro Offers Interpretable Disease Outcome Prediction by Leveraging Microbial Biomarkers and Their Time-Varying Effects. *bioRxiv (to appear).* **2024**

## Installation

- Available for Linux and Mac.
- For Windows users, please install Rtools prior to install LP-Micro.

```r
install.packages("devtools")
devtools:install_github("IV012/LPMicro")
```

## Usage

`LP-Micro` includes the following components:

- **Microbial Taxa Screening `screen_taxa`:** LP-Micro combines group lasso and polynomial splines to screen disease-correlated taxa from longitudinal microbiome data.

- **Visitwise Prediction `visit_predict`:** LP-Micro supports popular machine learning methods for predicting the future disease onset using microbial profiles at **a single time point**, namely visitwise prediction in our paper. User can idenfity the most predictive time point of microbiome for disease outcomes.

- **Cumulative Prediction `cumulative_predict`:** LP-Micro also supports the prediction of disease outcomes using microbial **up to a chosen visit,** namely cumulative prediction. Cumulative prediction is usually more accurate than visitwise prediction with proper feature engineering using `screen_taxa`.

- **Feature Interpretation `taxa_interpret`:** LP-Micro provides statistical p-values to test whether the following features are useful for prediction, including (i) a microbial taxon, (ii) all information from one time point, and (iii) a microbial taxon at a chosen time point.

Feel free to check more details using `help` and apply `LP-Micro` to your data!

### Simulation

```r

```