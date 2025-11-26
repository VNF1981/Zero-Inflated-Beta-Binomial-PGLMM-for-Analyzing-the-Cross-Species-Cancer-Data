# Zero-Inflated Beta Binomial PGLMM for Analyzing the Cross-Species Cancer Data 

This repository contains code for fitting zero-inflated beta-binomial phylogenetic generalized linear mixed models for cancer data across species. The models are implemented in R using the `brms` package with the `zero_inflated_beta_binomial()` family. The workflow reanalyzes the [Compton et al. 2024](https://doi.org/10.1158/2159-8290.CD-24-0573) dataset and model neoplasia and malignancy prevalence while accounting for overdispersion, excess zeros, variation in sampling effort across species, and phylogenetic structure.

## Model summary

### Purpose

The standard binomial model shows severe overdispersion when applied to current cross-species cancer data. This extra variance is not only limited to differences in necropsy numbers but also to innate inaccuracies in phylogenetic covariances, excess zeros, ecological and environmental factors, and other unmeasured biological differences across species.

This workflow uses a zero-inflated beta-binomial model with an explicit overdispersion parameter that can capture the excess variance in the data and account for extra zeros, sampling differences, and phylogenetic relatedness across species. The model includes

* a beta-binomial distribution to handle overdispersion

* a zero-inflation component to model extra zeros

* a phylogenetic random effect for shared ancestry

* life history predictors to test associations with cancer risk 

### Response and trials

The model uses a binomial-style structure

* `cases` is the number of naoplasia or malignancy cases in each species  
* `Trials` is the number of necropsies in that species  

The two responses should be modelled separately (see [Why Multivariate Model Would Not Work Here](#why-multivariate-model-would-not-work-here) below).

* `NeoplasiaCases`  
* `MalignancyCases`

## Model Structure
### Beta binomial component

The beta-binomial component models the underlying cancer probability for each species. The dispersion parameter widens the variance beyond the binomial limit. A variance inflation factor (VIF) is calculated from this parameter to summarize how much additional variance is present relative to a binomial model.

### Zero-Inflation Component

Many zeros in the dataset cannot be explained by the beta-binomial variance alone. These extra zeros often reflect low sampling effort. Two zero-inflation modes are available:

- `zi ~ log_trials_s` to link zero inflation to standardized log necropsy counts (default mode) 
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
3. gestation only (1-3 are simple regression models) 
4. all three traits together (multiple regression model)  

### Phylogenetic Random Effect

Because species are evolutionarily related through shared ancestry, their trait values cannot be treated as independent. The workflow reads a phylogenetic tree, prunes it to the species present in the data, constructs a variance-covariance matrix under a Brownian motion model, and passes this matrix to the species-level random effect. This produces a proper phylogenetic generalized linear mixed model.


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
- a custom variance inflation factor that summarizes **how much wider the actual variance is relative to what a standard binomial model could handle**. Please note that this is not the multicollinearity VIF commonly used in linear models and is calculated as `1 + (n_median − 1) / (phi + 1)`, which corresponds to `Var(beta-binomial) / Var(binomial)`.
- the formal overdispersion statistic `c`, calculated from Pearson residuals under a standard binomial assumption; **this quantifies how severely overdispersed the data would be if a binomial model were used instead of a beta-binomial model**
- a summary of the zero-inflation component  


Output files include:

- `summary_<model_label>_<response>.txt`  
- `standard_binomial_overdispersion_<model_label>_<response>.txt`  
- `BetaBinom_ZI_PGLMM_LHT_summary.csv` summarizing all models  

## Data Preparation and Tree Alignment

Before modeling, the script:

- reads `Compton_data.csv`  
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

#### *This script fits both simple models (single predictors) and the full multiple predictor model, and it runs each of them independently for the neoplasia and malignancy responses.*

## Why Multivariate Model Would Not Work Here

**1. The two response variables are not statistically independent**  
Malignancy is a strict subset of neoplasia, which creates a deterministic dependency between the two responses. Multivariate models assume stochastic dependence, not deterministic ones, so modeling both outcomes jointly violates the core assumptions of multivariate modeling and mis-specifies the joint likelihood and covariance structure.

**2. The scale and variance of the counts differ**  
Neoplasia counts can be an order of magnitude larger than malignancy counts and thus have a different dispersion level. A multivariate model assumes comparable variance scales and tries to estimate a shared covariance, random effects structure, and shared residual correlation. This assumes both responses operate on similar variance scales, which is false in the case of neoplasia and malignancy. 

**3. Overdispersion differs between neoplasia and malignancy**  
Neoplasia shows milder overdispersion, while malignancy is often extremely overdispersed (e.g., when using the full dataset or a subset of it). A multivariate model must link dispersion parameters and therefore forces an incorrect joint dispersion structure, which can cause non-convergence, parameter inflation, and unstable covariance estimates.

**4. The phylogenetic signal may differ between neoplasia and malignancy**  
Neoplasia may track general life-history traits, while malignancy may reflect deeper adaptations such as cancer resistance mechanisms (e.g., in elephants). A shared phylogenetic random effect misrepresents these differences.

**5. Multivariate zero-inflation is biologically incorrect**  
Zero inflation is not shared. A species can have many benign tumors and zero malignancies. Forcing a joint zero-inflation process assumes identical structural-zero mechanisms, which is false.

**6. Interpretation becomes incoherent**  
Shared predictor effects, shared covariance, and shared zero inflation have no biological meaning when one variable is always less than the other. The model becomes impossible to interpret cleanly. For example, we cannot interpret correlations when one variable is always less than the other, or we cannot explain species with neoplasia > 0 but malignancy = 0.
  

**7. The two outcomes have different biological processes**  
Neoplasia is shaped by general tumor initiation factors (e.g., mutation accumulation), while malignancy depends on invasion, metastasis, immune escape, and progression. Combining them forces the model to assume a shared covariance structure that does not fully match the actual biology.

***If we want to model neoplasia and malignancy together, the correct statistical structure is sequential, not multivariate. In that case, the proper framework would be conditional: first model tumor occurrence, then model malignancy given neoplasia. A multivariate model does not capture this hierarchy.***

## Why We Do Not Transform or Remove Outliers in a Beta-Binomial PGLMM

**1. Outlier removal is unnecessary and harmful for discrete count data**  
Outlier removal matters in Gaussian models where extreme values distort the mean and variance. In binomial and beta-binomial models, the response is bounded, discrete, and naturally includes many zeros or values near the boundaries. These are not statistical outliers. Removing species with high or low prevalence removes the real biological signal rather than noise.

**2. Transformations do not apply to count-based responses**  
Transforms like log, sqrt, Box-Cox, or arcsine are used for continuous responses to stabilize variance. They should not be applied to binomial counts or proportions because the distribution already accounts for boundedness and variance. Transforming counts breaks the likelihood and destroys interpretability. Predictors (such as body mass) may be transformed, but the response must remain untransformed.

**3. Misconceptions come from using linear models on proportions**  
Transformations and outlier filtering are common in linear models that treat prevalence as a continuous variable. In a beta-binomial PGLMM, the likelihood already handles overdispersion, boundedness, and skew. Zero inflation and phylogenetic structure absorb the remaining heterogeneity, making transformations unnecessary.

**4. The real challenges in cancer data are structural, not outliers**  
Issues such as overdispersion, phylogenetic correlation, excess zeros, and uneven sampling are inherent to cross-species cancer data. These are correctly handled by the dispersion parameter, the phylogenetic covariance matrix, and the zero-inflation component, and not by transforming the response or deleting observations.


```r
