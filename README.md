# mr.divw: Debiased Inverse-Variance Weighted Estimator in Summary-Data Mendelian Randomization and Multivariable Mendelian Randomization

*mr.divw* is an R package for summary-data Mendelian randomization (MR) and multivariable Mendelian randomization (MVMR) using the Debiased Inverse-Variance Weighted (dIVW) Estimator. The dIVW estimator is a simple modification of the popular IVW estimator, that has an explicit form and is robust to many weak instruments. The latest version allows for overlapping exposure and outcome datasets.

For details of the methods, please refer to our manuscripts [Ye. et al (2021)](https://arxiv.org/pdf/1911.09802.pdf) and [Wu. et al (2024)](https://arxiv.org/abs/2402.00307).

## Installation

To install the most up-to-date version, run the following command in R

```
library(devtools)
install_github("tye27/mr.divw")
```

## Examples
The main functions are *mr.divw* and *mvmr.divw* for respectively MR and MVMR analysis. You can find examples by running

```
library(mr.divw)
example(mr.divw)
example(mvmr.divw)
```

## Data Preprocessing

For two-sample summary-data univariable MR, summary statistics for SNP-exposure and SNP-outcome associations are enough for *mr.divw* to perform causal effect estimation. However, if the exposure and outcome dataset are overlapping, and/or if MVMR analysis is performed using *mvmr.divw*, we additionally require a shared correlation matrix as input (see the argument ```gen_cor```). This matrix can be assumed to be identity (i.e., assuming all exposure and outcome datasets are mutually independent), or fixed at identity, or estimated based on GWAS summary statistics. In the following, we use an example to illustrate how to use the [GRAPPLE](https://github.com/jingshuw/GRAPPLE) package to estimate the shared correlation matrix. For more details about installation of the [GRAPPLE](https://github.com/jingshuw/GRAPPLE) package and its usage, please visit https://github.com/jingshuw/GRAPPLE.

As an example, we are going to analyze the effects of certain HDL subfractions on the risk of coronary artery disease (CAD) with adjustments of the levels of traditional lipids, namely HDL, LDL and TG. One can download the summary statistics files for these traditional lipids, HDL subfractions, and CAD from datasets [link](https://www.dropbox.com/scl/fo/y0n8x81py5kxeiw97djsg/h?rlkey=dxeqkfvvrja02f2d0s52nyqs8&dl=0) and the reference panel [link](https://www.dropbox.com/scl/fo/cucd65mredj3kl5yukmo3/h?rlkey=zvzn3pc33zb0gt8e4tco9fdfd&dl=0). The data were from public GWASs (see our manuscript [link] for more details). We thank [Dr. Qingyuan Zhao](https://www.statslab.cam.ac.uk/~qz280/) for sharing the data. 

Then, we can call the function *getInput* from the [GRAPPLE](https://github.com/jingshuw/GRAPPLE) package and run the following code. 

```
library(GRAPPLE)
traits <- c("s_hdl_p", "s_hdl_tg", "m_hdl_p", "m_hdl_c", "l_hdl_p", "l_hdl_c")
sel.files <- c(paste0("data_csv/gera_", c("hdl","ldl", "tg"), ".csv"),
                   paste0("data_csv/davis_", traits, ".csv"))
exp.files <- c(paste0("data_csv/willer_", c("hdl","ldl", "tg"), ".csv"),
                   paste0("data_csv/kettunen_", traits, ".csv"))
out.file <- paste0("data_csv/cardiogramplusc4d_ukbb_cad.csv")
hdl_subfractions <- getInput(sel.files, exp.files, out.file,
                plink_exe = './plink', plink_refdat = "ld_files/data_maf0.01_rs", 
                max.p.thres = 1 * (10^(-4)), cal.cor = TRUE, get.marker.candidates = FALSE)
save(hdl_subfractions, file = "hdl_subfractions.rda")
```

According to [GRAPPLE](https://github.com/jingshuw/GRAPPLE), this function extracts all independent SNPs whose selection p-values do not exceed 1e-04. It returns three elements. One is 'data' which is a data frame of the summary statistics of the selected SNPs, the other is 'marker.data' which is a data frame for all candidate marker SNPs (this will be empty, because we set ```markder.data = FALSE```). The third element is the estimated correlation matrix for the GWAS cohorts, which can then be used as the input for *mr.divw* and *mvmr.divw* in the case with overlapping datasets.

Note that [GRAPPLE](https://github.com/jingshuw/GRAPPLE) uses PLINK to perform LD clumping. The default name/path of the plink command is "plink" (passed through the argument ```plink_exe```). For users with Linux / Windows systems, this command would not work and one need to specify the path of the exe file, like "./plink" depending on where they install plink. If R fails to run the plink command, there will be an error stating that the clumped file is not found.

The *hdl_subfractions.rda* is used as an example dataset for *mr.divw* and *mvmr.divw*.




