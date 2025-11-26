# Zero-Inflated-Beta-Binomial-PGLMM-for-Analyzing-the-Cross-Species-Cancer-Data
A workflow for fitting zero-inflated beta binomial phylogenetic GLMMs that model cross species cancer data while handling overdispersion, excess zeros, sampling effort, and phylogenetic nonindependence.

# Zero Inflated Beta Binomial PGLMM  
Cancer prevalence across species

This repository contains code for fitting zero inflated beta binomial phylogenetic generalized linear mixed models for cancer data across species. The goal is to model neoplasia and malignancy prevalence while accounting for overdispersion, extra zeros, sampling effort, and phylogenetic structure.

## Model summary

### Purpose

Cancer counts across species show strong overdispersion and many zero values. Simple binomial models and non phylogenetic models do not fit these data.

This workflow uses

* a beta binomial distribution to handle extra variation  
* a zero inflation component to model extra zeros  
* a phylogenetic random effect for shared ancestry  
* life history predictors to test associations with cancer risk  

### Response and trials

The model uses a binomial style structure

* `cases` is the number of cancer cases in each species  
* `Trials` is the number of necropsies in that species  

Two responses are modelled separately

* `NeoplasiaCases`  
* `MalignancyCases`  

### Beta binomial component

The beta binomial part models the underlying cancer probability for each species. It includes a dispersion parameter `phi` that allows the variance to be larger than in a simple binomial model.

From `phi` the script computes a variance inflation factor

* VIF shows how much wider the variance is relative to a binomial model  
* VIF is based on `phi` and the median number of necropsies per species  

### Zero inflation component

The data contain more zero values than expected under the beta binomial variance. This is often due to low sampling effort for some species.

The zero inflation part accounts for this and supports two modes

* `zi ~ log_trials_s`  
  * zero inflation depends on standardized log necropsy counts  
  * species with fewer necropsies can have a higher probability of extra zeros  

* `zi ~ 1`  
  * a single global zero inflation level without covariates  

The choice is controlled by the `zi_mode_global` setting in the script.

### Predictors and life history traits

Life history traits are used as predictors in the mean part of the model

* body mass  
* maximum longevity  
* gestation length  

In the script these are

* log transformed  
* standardized to mean zero and unit variance  

The code fits four predictor sets for each response

* mass only  
* longevity only  
* gestation only  
* all three traits together  

This allows direct comparison of single trait and multi trait models.

### Phylogenetic random effect

Species are not independent. Closely related species tend to have similar cancer probabilities.

The script builds a phylogenetic covariance matrix `A` from the pruned tree using Brownian motion and uses it in a random effect term

```r
(1 | gr(Species, cov = A_model))
