# Zero-Inflated Beta Binomial PGLMM for Analyzing the Cross-Species Cancer Data 

This repository contains code for fitting zero-inflated beta-binomial phylogenetic generalized linear mixed models for cancer data across species. The models are implemented in R using the `brms` package with the `zero_inflated_beta_binomial()` family. The workflow is built to reanalyze the [Compton et al. 2024](https://doi.org/10.1158/2159-8290.CD-24-0573) dataset and model neoplasia and malignancy prevalence while accounting for overdispersion, excess zeros, variation in sampling effort across species, and phylogenetic structure.

## Model summary

### Purpose

The standard binomial model shows severe overdispersion when applied to cross-species cancer data. This extra variance is not only due to differences in necropsy numbers but also to phylogenetic structure, excess zeros, ecological factors, and unmeasured biological differences across species.

This workflow uses a zero-inflated beta-binomial model with an explicit overdispersion parameter that can capture the excess variance in the data and account for extra zeros, sampling differences, and phylogenetic relatedness across species. The model includes

* a beta-binomial distribution to handle overdispersion

* a zero-inflation component to model extra zeros

* a phylogenetic random effect for shared ancestry

* life history predictors to test associations with cancer risk 

### Response and trials

The model uses a binomial-style structure

* `cases` is the number of cancer cases in each species  
* `Trials` is the number of necropsies in that species  

Two responses should be modelled separately (see [Why Multivariate Model Would Not Work Here](#why-multivariate-model-would-not-work-here) below).

* `NeoplasiaCases`  
* `MalignancyCases`

## Model Structure
### Beta binomial component

The beta-binomial component models the underlying cancer probability for each species. The dispersion parameter widens the variance beyond the binomial limit. A variance inflation factor (VIF) is calculated from this parameter to summarize how much additional variance is present relative to a binomial model.

### Zero-Inflation Component

Many zeros in the dataset cannot be explained by the beta-binomial variance alone. These extra zeros often reflect low sampling effort. Two zero-inflation modes are available:

- `zi ~ log_trials_s` to link zero inflation to standardized log necropsy counts  
- `zi ~ 1` for a global zero-inflation level  

The choice is controlled by the `zi_mode_global` setting in the script.

### Life History Predictors

Life history traits are used as predictors in the mean model. These include:

- adult body mass  
- maximum longevity  
- gestation length  

They are log-transformed and standardized before modeling.

The workflow fits four predictor sets for each response:

1. mass only  
2. longevity only  
3. gestation only  
4. all three traits together  

### Phylogenetic Random Effect

Species are evolutionarily related and cannot be treated as independent. The workflow reads a phylogenetic tree, prunes it to species present in the data, constructs a correlation matrix under a Brownian motion model, and passes this matrix to the species-level random effect. This produces a proper phylogenetic generalized linear mixed model.

## Priors

Weakly informative priors are used for:

- intercepts in the mean and zero-inflation parts  
- slopes for life history predictors  
- the standard deviation of the species random effect  
- the beta-binomial dispersion parameter  

These priors stabilize estimation while remaining data driven.

## Model Fitting

The models are fitted in R with `brms` using the `zero_inflated_beta_binomial()` family. The workflow uses four chains and conservative sampler settings:

- increased warmup  
- `adapt_delta = 0.99`  
- `max_treedepth = 14`  

These settings reduce divergent transitions and improve convergence.

## Diagnostics and Output

For each fitted model, the workflow calculates:

- posterior odds ratios with 95 percent credible intervals  
- maximum Rhat values across all parameter groups  
- the variance inflation factor summary  
- the overdispersion statistic based on Pearson residuals  
- a summary of the zero-inflation component  

Output files include:

- `summary_<model_label>_<response>.txt`  
- `overdispersion_<model_label>_<response>.txt`  
- `BetaBinom_ZI_PGLMM_LHT_summary.csv` summarizing all models  

## Data Preparation and Tree Alignment

Before modeling, the script:

- reads `Zach_data.csv`  
- cleans species names  
- aligns species with the phylogenetic tree  
- prunes the tree to the intersecting species  
- reorders rows to match the tree tip order  
- constructs log-transformed and standardized predictors  
- constructs `log_trials_s` for zero-inflation  

Only species with complete data for a given model are included.

## Main Script Structure

The main script defines the function `fit_lht_model()` which performs:

- data subsetting and alignment  
- formula construction for both the mean and zero-inflation parts  
- model fitting via `brm`  
- extraction of posterior summaries  
- writing of output files  

The script loops over all predictor sets and both responses. To run everything:

```r
source("your_script_name.R")
