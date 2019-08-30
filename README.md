Dan McGuire
2019-08-30

MAMBA
=====

A **M**eta-**A**nalysis **M**odel **B**ased **A**ssessment of Replicability for genome-wide association studies.

The main goal of the MAMBA model is to identify loci which have non-zero replicable associations with a phenotype based on summary statistics from genome-wide association meta-analysis (GWAMA).

Installation
------------

    library(devtools)
    install_github("dan11mcguire/mamba")

Quick Tutorial
--------------

The MAMBA model is fit with the `mamba()` function, and takes as input a matrix of effect size estimates and a matrix of their variances from GWA meta-analysis.

For a quick start, one could use

    d<-generate_data_mamba()

to generate summary statistics according to the MAMBA model (SNPs are simulated independently). The returned data `d$betajk` is a *M*x*K* matrix of effect size estimates, where *M* is the total number of SNPs and *K* is the total number of studies, with the *j*th row and *k*th column correspond to the effect size estimate of the *j*th SNP in the *k*th study. Similarly, `d$sjk2` is a *M*x*k* matrix of the effect size estimate variances.

The MAMBA model could then be fit by:

    mod<-mamba(betajk=d$beta, sjk2=d$sjk2)

The `mamba()` function then fits a two-level mixture model to the data, and calculates a posterior-probability-of-replicability (PPR) for each SNP. This is given in `mod$ppr`.

For SNPs with low PPR, the model also estimates a posterior probability that a particular summary statistic is an "outlier". We use the term "outlier" here to indicate a summary statistic which may be extreme or significant due only to artifacts. This is given in the `mod$outliermat` matrix.

Detailed example
----------------

### Step 0: Initial GWAMA

To start, we conduct an initial fixed-effects meta-analysis to identify loci of interest which we want to validate with the MAMBA model. This could be done using METAL, rareGWAMA, or some other standard-type meta-analysis software for GWAS. As an example, we have the helper function `make_mamba_summary_files()` which takes as input a set of summary statistic files, and conducts inverse-variance-weighted fixed effects meta-analysis. We can use some example summary statistic files available as external package datasets for the example analysis.

``` r

library(mamba)
#> Loading required package: TMB
#> Loading required package: RcppEigen
library(data.table)
library(parallel)

datadir<-system.file("extdata", package="mamba")
datafiles<-paste0(datadir, "/", grep("study_.*assoc", list.files(datadir),value=T))

fread(datafiles[1])[1:5] ### format of example summary statistic files 
#>    chr       bp ref alt    n        af        betaj         sj2          p
#> 1:  14 19012028   T   G 3873 0.2109570 -0.042544349 0.003150100 0.44844000
#> 2:  14 19077386   G   A 3873 0.4448340 -0.009459681 0.002123235 0.83734200
#> 3:  14 19108980   T   A 3873 0.0109496  0.590080092 0.048417489 0.00732496
#> 4:  14 19111033   G   A 3873 0.0215132 -0.057194393 0.024909152 0.71706200
#> 5:  14 19118923   A   C 3873 0.1010160  0.060027754 0.005774000 0.42954200

## holdout the first dataset to assess replicability
discovery_datasets<-setdiff(datafiles, grep("study_1.assoc", datafiles, value=T))
exmpl<-make_mamba_summary_files(summary_files=discovery_datasets, 
                                standardize_beta_and_s2_datasets=TRUE,
                                beta_name="betaj", s2_name="sj2")

assoc<-exmpl$assoc
beta<-exmpl$beta
s2<-exmpl$s2
rm(exmpl)
```

`assoc` contains the meta-analysis results (in `plink` .assoc file format, with some additional columns) which can be passed directly to `plink` in step2.

<https://www.cog-genomics.org/plink/1.9/formats>

``` r

assoc[1:3] 
#>    CHR       BP             SNP A1 A2        FRQ INF0         BETA
#> 1:  14 19012028 14:19012028_T_G  G  T 0.19705563    1 -0.003134722
#> 2:  14 19077386 14:19077386_G_A  A  G 0.43236777    1 -0.001605019
#> 3:  14 19108980 14:19108980_T_A  A  T 0.01007469    1 -0.022115559
#>             SE         P  k          z  tot_n     tot_ac
#> 1: 0.006685018 0.6391285 11 -0.4689176 304200 119888.647
#> 2: 0.005399108 0.7662567 11 -0.2972748 304200 263052.553
#> 3: 0.026866038 0.4104062  9 -0.8231791 295401   5952.149
```

The `beta` and `s2` objects will later be used as input to the `mamba()` function.
Each row represents the effect size estimates (or their variances) of a SNP, and each column represents a particular study. Not every SNP is analyzed in every study, missing values are represented by `NA`.

``` r

beta[1:2]
#>    chr       bp             snp ref alt          b_1       b_2        b_3
#> 1:  14 19012028 14:19012028_T_G   T   G -0.024547326 -4.477517 0.03258657
#> 2:  14 19077386 14:19077386_G_A   G   A -0.006648167 -3.235424 0.11786526
#>            b_4 b_5 b_6 b_7 b_8 b_9 b_10 b_11 b_12 b_13      b_14 b_15 b_16
#> 1: -0.01951356  NA  NA  NA  NA  NA   NA   NA   NA   NA -7.512437   NA   NA
#> 2: -0.03475506  NA  NA  NA  NA  NA   NA   NA   NA   NA -2.333415   NA   NA
#>    b_17       b_18 b_19 b_20 b_21 b_22 b_23       b_24         b_25 b_26
#> 1:   NA 0.04751072   NA   NA   NA   NA   NA -0.1918304 -0.001492210   NA
#> 2:   NA 0.02447689   NA   NA   NA   NA   NA  0.7427556 -0.000931085   NA
#>    b_27 b_28      b_29 b_30 b_31 b_32      b_33 b_34       b_35
#> 1:   NA   NA 0.2870923   NA   NA   NA 30.391214   NA -0.7272565
#> 2:   NA   NA 0.2439959   NA   NA   NA  1.476701   NA -0.3871208
s2[1:2]
#>    chr       bp             snp ref alt        s2_1     s2_2        s2_3
#> 1:  14 19012028 14:19012028_T_G   T   G 0.001048694 23.97333 0.004037369
#> 2:  14 19077386 14:19077386_G_A   G   A 0.001048694 22.75360 0.004037369
#>            s2_4 s2_5 s2_6 s2_7 s2_8 s2_9 s2_10 s2_11 s2_12 s2_13    s2_14
#> 1: 0.0006575349   NA   NA   NA   NA   NA    NA    NA    NA    NA 5.943251
#> 2: 0.0006575349   NA   NA   NA   NA   NA    NA    NA    NA    NA 4.809094
#>    s2_15 s2_16 s2_17       s2_18 s2_19 s2_20 s2_21 s2_22 s2_23     s2_24
#> 1:    NA    NA    NA 0.001979334    NA    NA    NA    NA    NA 0.9839675
#> 2:    NA    NA    NA 0.001979334    NA    NA    NA    NA    NA 0.4770850
#>           s2_25 s2_26 s2_27 s2_28     s2_29 s2_30 s2_31 s2_32      s2_33
#> 1: 1.499739e-05    NA    NA    NA 0.7574003    NA    NA    NA 8776.17643
#> 2: 1.499739e-05    NA    NA    NA 0.5295231    NA    NA    NA   20.84973
#>    s2_34    s2_35
#> 1:    NA 3.600483
#> 2:    NA 1.029959
```

We can subset the analysis to SNPs which were analzed in at least 4 studies. We write out the association file to filepath: `ivw_met_chr14.tsv`.

``` r
assoc<-assoc[k >=4]
beta<-beta[snp %in% assoc[,SNP]]
s2<-s2[snp %in% assoc[,SNP]]


fwrite(assoc, file="ivw_met_chr14.tsv", sep="\t")
```

### Step 1a: Clumping GWAMA results

<https://www.cog-genomics.org/plink/1.9/postproc#clump>

<https://www.cog-genomics.org/plink/1.9/formats>

An example command for clumping meta-analysis results in plink is below:

``` bash

plink \
  --bed "$refpanel".bed \
  --bim "$refpanel".bim \
  --fam "$refpanel".fam \
  --clump ivw_met_chr14.tsv \
  --clump-p1 1e-5 \
  --clump-kb 500 \
  --clump-p2 1 \
  --clump-r2 0.1 \
  --out example_chr14

## outputted file with clumped association statistics: example_chr14.clumped
```

An output file from clumping using the above command is available as an external package dataset.

``` r

datadir<-system.file("extdata", package="mamba")
clump_file<-paste0(datadir, "/", "example_chr14.clumped")
tail(fread(clump_file))
#>    CHR F              SNP        BP        P TOTAL NSIG S05 S01 S001 S0001
#> 1:  14 1  14:69803439_A_T  69803439 7.58e-06     9    0   0   0    4     5
#> 2:  14 1 14:105904300_G_A 105904300 8.13e-06     0    0   0   0    0     0
#> 3:  14 1 14:105906937_T_G 105906937 9.23e-06     0    0   0   0    0     0
#> 4:  14 1  14:27126790_A_C  27126790 9.49e-06    40    0   0   0   28    12
#> 5:  14 1 14:105893876_G_A 105893876 9.61e-06     0    0   0   0    0     0
#> 6:  14 1 14:101025848_C_T 101025848 9.79e-06     0    0   0   0    0     0
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        SP2
#> 1:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              14:69801864_C_A(1),14:69803200_G_A(1),14:69804924_C_T(1),14:69806055_A_T(1),14:69807319_G_A(1),14:69808627_G_C(1),14:69809624_G_A(1),14:69813370_A_G(1),14:69816344_G_A(1)
#> 2:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    NONE
#> 3:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    NONE
#> 4: 14:26935305_G_A(1),14:26937116_G_A(1),14:26954078_T_C(1),14:26958210_T_G(1),14:26961879_T_C(1),14:26963559_T_G(1),14:26981079_G_A(1),14:26983770_C_G(1),14:26985347_A_C(1),14:26985560_G_T(1),14:26986041_G_A(1),14:26994171_C_T(1),14:26996510_C_T(1),14:27009199_A_C(1),14:27010629_C_A(1),14:27016744_T_C(1),14:27017018_T_C(1),14:27017537_A_G(1),14:27021647_T_C(1),14:27030207_T_C(1),14:27038915_T_C(1),14:27051340_A_T(1),14:27051498_A_T(1),14:27056033_T_C(1),14:27056252_C_T(1),14:27057492_T_C(1),14:27105035_C_G(1),14:27105474_A_T(1),14:27107166_G_A(1),14:27107260_A_G(1),14:27107921_T_A(1),14:27111797_G_T(1),14:27112102_A_G(1),14:27113632_G_C(1),14:27114262_C_G(1),14:27114773_T_C(1),14:27114785_C_T(1),14:27115116_C_T(1),14:27120136_G_A(1),14:27128033_A_T(1)
#> 5:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    NONE
#> 6:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    NONE
```

### Step 1b: Prune variants from a reference panel, combine with clumped SNPs.

<https://www.cog-genomics.org/plink/1.9/ld>

An example command for pruning from a reference panel using plink is shown below.

``` bash

plink \
  --bed "$refpanel".bed \
  --bim "$refpanel".bim \
  --fam "$refpanel".fam \
  --indep-pairwise 500kb 1 0.1 \
  --maf 0.01 \
  --out example_chr14 

## outputted file with pruned snps:  example_chr14.prune.in
```

A set of SNPs pruned from HRC reference panel using the below commands (in *CHR:BP\_REF\_ALT* format) are available as a default in the MAMBA package.

``` r
data(prunedsnps_HRC, package="mamba")
prunevars[1:5]
#>              SNP
#> 1:  10:65030_C_A
#> 2:  10:94263_C_A
#> 3:  10:95248_C_G
#> 4:  10:99932_T_G
#> 5: 10:112184_A_G
```

We can use the helper function `select_mamba_loci()` to combine the pruned and clumped variants, removing any pruned SNPs within 500kb of a clumped SNP.

``` r

mamba_snps<-select_mamba_loci(meta_file="ivw_met_chr14.tsv", clump_file=clump_file)
#> [1] "input read.  creating rema dataset.."
#> Time difference of 0.01489162 secs

betajk<-beta[snp %in% mamba_snps[,SNP]]
sjk2<-s2[snp %in% mamba_snps[,SNP]]
```

### Step 2: Fit the MAMBA model.

``` r

mod<-mamba(betajk=betajk[,-c("chr", "bp", "snp", "ref", "alt"),with=F], 
           sjk2=sjk2[,-c("chr", "bp", "snp", "ref", "alt"),with=F],
           snpids=betajk[,snp])
#> [1] "e step finished."
#> [1] "m step finished."
#> [1] "p=0.005 lambda=0.959 tau2=0.00014673 alpha=9.657"
#> [1] "iteration 1 complete."
#> [1] "e step finished."
#> [1] "m step finished."
#> [1] "p=0.007 lambda=0.969 tau2=0.00010031 alpha=12.592"
#> [1] "relative ll increase: 0.00289249858605164"
#> [1] "ll increase: 691.81104935435"
#> [1] "iteration 2 complete."
#> [1] "e step finished."
#> [1] "m step finished."
....

....

#> [1] "e step finished."
#> [1] "m step finished."
#> [1] "p=0.032 lambda=0.995 tau2=3.835e-05 alpha=68.676"
#> [1] "relative ll increase: 1.0518675706234e-08"
#> [1] "ll increase: 0.00253462721593678"
#> [1] "iteration 45 complete."
#> [1] "e step finished."
#> [1] "m step finished."
#> [1] "p=0.032 lambda=0.995 tau2=3.826e-05 alpha=68.676"
#> [1] "relative ll increase: 9.5648508132523e-09"
#> [1] "ll increase: 0.00230478931916878"
```

### Assessing the output

The estimated posterior-probability-of-replicability (PPR) and effect sizes are contained in the `mod$ppr` and `mod$mu.hat` objects, which we reattach to our initial fixed effects meta-analysis table for comparison.

``` r

mamba_snps[,ppr:=mod$ppr]
mamba_snps[,mu_hat:=mod$mu.hat]
```

We can use the `mod$outliermat` object to identify the "most probable outlier study" for non-replicable SNPs. The columns in `mod$outliermat` assess posterior probability of an outlier *given* that true SNP association is zero and non-replicable. The number is only interpretable for SNPs with low PPR, where higher values indicate a summary statistic which dominated the meta-analysis result.

``` r
mod$outliermat[ppr < 0.5][1:3]
#>                snp         Oj_1         Oj_2         Oj_3         Oj_4
#> 1: 14:19012028_T_G 0.0007593869 0.0008637992 0.0006513838 0.0007610746
#> 2: 14:19077386_G_A 0.0005842565 0.0007177489 0.0031102118 0.0014135876
#> 3: 14:19108980_T_A 0.0194177552 0.0011617249 0.0008604968 0.0040570477
#>    Oj_5 Oj_6 Oj_7 Oj_8 Oj_9 Oj_10 Oj_11 Oj_12 Oj_13        Oj_14 Oj_15
#> 1:   NA   NA   NA   NA   NA    NA    NA    NA    NA 0.0580568694    NA
#> 2:   NA   NA   NA   NA   NA    NA    NA    NA    NA 0.0009992573    NA
#> 3:   NA   NA   NA   NA   NA    NA    NA    NA    NA 0.0006687855    NA
#>    Oj_16 Oj_17        Oj_18 Oj_19 Oj_20 Oj_21 Oj_22 Oj_23        Oj_24
#> 1:    NA    NA 0.0010033110    NA    NA    NA    NA    NA 0.0005828923
#> 2:    NA    NA 0.0006642334    NA    NA    NA    NA    NA 0.0010112191
#> 3:    NA    NA 0.0006028209    NA    NA    NA    NA    NA 0.0006064193
#>           Oj_25 Oj_26 Oj_27 Oj_28        Oj_29 Oj_30 Oj_31 Oj_32
#> 1: 0.0006156617    NA    NA    NA 0.0006037579    NA    NA    NA
#> 2: 0.0005887791    NA    NA    NA 0.0006048313    NA    NA    NA
#> 3: 0.0008110691    NA    NA    NA 0.0006554179    NA    NA    NA
#>           Oj_33 Oj_34        Oj_35        ppr
#> 1: 0.0006026946    NA 0.0006151846 0.01787158
#> 2: 0.0006025006    NA 0.0006147624 0.01805945
#> 3:           NA    NA           NA 0.02179809
top_outlier<-melt(mod$outliermat, 
                  id.vars=c("ppr", "snp"), 
                  variable.factor=F,
                  value.name="outlier_prob",
                  variable.name="outlier_study")[order(-outlier_prob)][,head(.SD,1),by=.(snp)]

mamba_snps<-merge(mamba_snps, top_outlier[,-c("ppr"),with=F], by.x="SNP", by.y="snp",sort=FALSE)
```

#### Estimate FDR from PPR

Dr. Brad Efron has written extensively on the connection between empirical bayes posterior probability and false discovery rates, for example here: [Empirical Bayes Analysis of a Microarray Experiment](https://pdfs.semanticscholar.org/21b4/c0579e8c7398cfa2f5c476b659090fea9209.pdf)

Below we calculate an estimated FDR based on posterior probability (PPR) from the MAMBA model.

``` r

mamba_snps[,ppr_rank:=rank(-ppr,ties.method="max")]
mamba_snps[order(-ppr),Fdr:=cumsum(1-ppr)/ppr_rank]
mamba_snps[,Fdr:=max(Fdr),by=.(ppr_rank)]
mamba_snps[,ppr_rank:=NULL]

mamba_snps[order(-ppr),.(SNP, p_fe=P, ppr, Fdr)][1:5]
#>                SNP         p_fe       ppr          Fdr
#> 1: 14:29273634_G_A 5.653125e-09 0.9998325 0.0001674948
#> 2: 14:28356733_G_A 4.674038e-08 0.9996620 0.0002527490
#> 3: 14:29795213_G_A 1.577426e-07 0.9992481 0.0004191414
#> 4: 14:32438484_G_A 1.349447e-07 0.9986866 0.0006427105
#> 5: 14:29791660_G_C 8.369951e-07 0.9978895 0.0009362671
```

#### Calculate p-values

We use parametric bootstrap procedure to calculate p-values based on estimated PPR.
Basically, this means

1.  Simulating data based on the parameters we estimated in `mod`.
    -   The function `get_null_scores()` using the call below outputs PPR for non-replicable SNPs into the file `mamba_pval/chr14.nullscores.txt`.

    -   *Note*: The argument `total_null_scores` will affects the minimum p-value you can possibly obtain. For assessing a SNP at genome-wide significance, you need at least 20 × 10<sup>6</sup> null scores, $\\frac{1}{(20\\times 10^6+1)} \\approx 5\\times 10^{-8}$

    -   *Note*: `clean_any_existing_scores=FALSE`. When set to FALSE, simultaneous calls to `get_null_scores()` can dump output to the same file `mamba_pval/chr14.nullscores.txt`. This can speed up computation when additional computing resources are available.

2.  Comparing the PPR from the simulated non-replicable SNPs (data generated under the null hypothesis), with the estimated PPR from our original data. For a SNP *i* with PPR score=*y* and a set of simulated null SNP PPR scores *x*, the bootstrap p-value can be calculated as
    $$p\_i = \\frac{1 + \\sum (y &lt; x)}{1 + |x|}$$
     , where |*x*| is the number of simulated null SNP PPR scores.

For this short example, we settle for 1 million (10^6) scores. With `clean_any_existing_scores=FALSE`, we see some have already been generated by a prior call to the function.

``` r
system("if [ ! -d mamba_pvals ]; then mkdir mamba_pvals/ ; fi ")
null_scores<-mamba:::get_null_scores(model=mod, 
                                     s2=sjk2[,-c("chr", "bp", "ref", "alt", "snp"),with=F], 
                                     total_null_scores=10^6, 
                                     out="mamba_pvals/chr14",
                                     clean_any_existing_scores=FALSE)
#> [1] "990000 scores already generated in mamba_pvals/chr14.nullscores.txt."
#> [1] "Fitting 2 models to generate approximately 1e+06 non-replicable SNPs."
#> [1] "Model 1 of 2 complete."
#> [1] "Model 2 of 2 complete."
```

``` r
mod$pval<-sapply(mod$ppr, 
                 function(score) (1 + sum(score < null_scores[[1]]))/(1 + length(null_scores[[1]])))

mamba_snps[,mamba_p:=mod$pval]
mamba_snps[order(mamba_p),.(SNP, mamba_p, ppr, Fdr)][1:10]
#>                 SNP      mamba_p       ppr          Fdr
#>  1: 14:28356733_G_A 9.985312e-07 0.9996620 0.0002527490
#>  2: 14:29273634_G_A 9.985312e-07 0.9998325 0.0001674948
#>  3: 14:29795213_G_A 9.985312e-07 0.9992481 0.0004191414
#>  4: 14:32438484_G_A 9.985312e-07 0.9986866 0.0006427105
#>  5: 14:29791660_G_C 1.997062e-06 0.9978895 0.0009362671
#>  6: 14:79580218_C_T 1.997062e-06 0.9960663 0.0017957679
#>  7: 14:90057852_C_T 1.997062e-06 0.9977493 0.0011553457
#>  8: 14:98594361_A_C 1.997062e-06 0.9964996 0.0014903512
#>  9: 14:77592879_G_A 2.995593e-06 0.9953230 0.0021159029
#> 10: 14:98811932_A_G 2.995593e-06 0.9942386 0.0024804534
```

#### Comparison with holdout dataset

We can make some further comparisons with the initial fixed effects result with the results from the MAMBA model by looking at the holdout dataset excluded from our analysis.

``` r
### attach held-out dataset to results
holdout<-fread(paste0(datadir, "/study_1.assoc"))

mamba_snps<- merge(mamba_snps, 
             holdout[,.(snp=paste0(chr, ":", bp, "_", ref, "_", alt), beta_rep=betaj, p_rep=p)], 
             by.x="SNP", by.y="snp", all.x=T, sort=F)


### Proportion of SNPs with PPR > 0.5 which replicate direction of effect in the holdout cohort 
mamba_snps[ppr > 0.5][!is.na(p_rep),mean(sign(beta_rep)==sign(mu_hat))]
#> [1] 0.8965517

### Proportion of SNPs with PPR > 0.5 which replicate direction of effect and have 
### replication p-value < 0.01 in the holdout cohort 
mamba_snps[ppr > 0.5][!is.na(p_rep),mean(p_rep < .01 & (sign(beta_rep)==sign(mu_hat)))]
#> [1] 0.3793103


### Comparing replication using a fixed effects p-value threshold corresponding to the 
### PPR cutoff=0.5 used above.
max_fe_p<-mamba_snps[ppr > 0.5][,max(P)]
mamba_snps[P < max_fe_p][!is.na(p_rep),mean(sign(beta_rep)==sign(mu_hat))]
#> [1] 0.5849515
mamba_snps[P < max_fe_p][!is.na(p_rep),mean(p_rep < .01 & (sign(beta_rep)==sign(mu_hat)))]
#> [1] 0.03883495


### Replication for SNPs with likely outlier summary statistics 
mamba_snps[outlier_prob > 0.5][ppr < 0.5][index_snp==1][!is.na(p_rep),mean(sign(beta_rep)==sign(mu_hat))]
#> [1] 0.5639687
mamba_snps[outlier_prob > 0.5][ppr < 0.5][index_snp==1][!is.na(p_rep),mean(p_rep < .05 & (sign(beta_rep)==sign(mu_hat)))]
#> [1] 0.04699739


### Correlation between mamba effects size estimates 
### and replication effect size estimates, 
mamba_snps[!is.na(beta_rep),cor(mu_hat, beta_rep)] 
#> [1] 0.02926117

### Correlation between fixed effects size estimates 
### and replication effect size estimates, 
mamba_snps[!is.na(beta_rep), cor(BETA, beta_rep)]
#> [1] 0.004871806
```
