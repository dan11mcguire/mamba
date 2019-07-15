# MAMBA 
A **M**eta-**A**analysis **M**odel **B**ased **A**ssessment of Replicability for genome-wide association studies.

## Installation

    library(devtools)
    install_github("dan11mcguire/mamba")

## Quick Tutorial 

The MAMBA model is fit with the `mamba()` function, and takes as input a matrix of  effect size estimates and a matrix of their variances from GWA meta-analysis. 

For a quick start, one could use 

    d<-generateData_mamba()
 
To generate summary statistics according to the MAMBA model (SNPs are simulated independently).  The returned data `d$betajk` is a *M*x*K* matrix of effect size estimates, where *M* is the total number of SNPs and *K* is the total number of studies, with the *j*th row and *k*th column correspond to the effect size estimate of the *j*th SNP in the *k*th study.  Similarly, `d$sjk2` is a *M*x*k* matrix of the effect size estimate variances.  

The MAMBA model could then be fit by:

    mod<-mamba(betajk=d$beta, sjk2=d$sjk2)


The `mamba()` function then fits a two-level mixture model to the data, and calculates a posterior-probability-of-replicability (PPR) for each SNP.  For SNPs with low PPR, the model also estimates a posterior probability that a particular summary statistic is an ``outlier''. We use the term ``outlier'' here to indicate a summary statistic which is significant,  


## Analysis Pipeline 

**Step 0**: First perform a fixed-effect GWA meta-analysis to identify loci of interest with suggestive evidence of association (For example, could use a threshold such as *p < 1\times 10^-5* 

**Step 1a**: Prune variants with suggestive evidence of association using the clumping procedure implemented in Plink v1.9.

**Step 1b**: Given that the SNPs from Step 1a will all initially appear to have non-zero effects from the fixed-effects meta-analysis, we incorporate summary statistics from an independent set of variants randomly pruned based upon a reference panel.  These SNPs allow the non-replicable, zero-mean component of the MAMBA mixture model to be reliably estimated.

**Step 2**: Fit the MAMBA model










