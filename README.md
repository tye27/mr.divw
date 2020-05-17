# mr.divw: Debiased Inverse-Variance Weighted Estimator in Two-Sample Summary-Data Mendelian Randomization
Debiased Inverse-Variance Weighted Estimator in Two-Sample Summary-Data Mendelian Randomization

## Setup
*mr.divw* is an R package for two-sample summary-data Mendelian randomization using the Debiased Inverse-Variance Weighted (dIVW) Estimator. 
The dIVW estimator is a simple modification of the popular IVW estimator, that has an explicit form and is robust to many weak instruments. 

To install the most up-to-date version, run the following command in R

```
library(devtools)
install_github("tye27/mr.divw")
```

## Examples
The main function is *mr.divw*. You can find examples by running

```
library(mr.divw)
example(mr.divw) 
```





