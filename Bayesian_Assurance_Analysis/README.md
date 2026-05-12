# Bayesian Assurance Analysis for ZiBBPGLMM

## Overview

This repository contains R scripts for fitting a zero inflated beta binomial phylogenetic generalized linear mixed model, abbreviated here as ZiBBPGLMM, to cross species cancer data. The main goal is to test associations between cancer occurrence and life history traits while accounting for uneven sampling effort, excess zeros, overdispersion, and phylogenetic non independence among species.

The repository also includes a Bayesian assurance analysis. The assurance analysis uses simulation to ask whether the current dataset has enough information to detect a given effect size under the same model structure.

The current implementation focuses on two cancer related responses:

1. Neoplasia cases
2. Malignancy cases

The main predictors are:

1. Adult body mass
2. Maximum longevity
3. Gestation length

The current assurance script is set up to test body mass by default, while controlling for longevity and gestation in the full three predictor model.

## Repository structure

```text
ZiBBPGLMM/
├── 03_zero_inflated_beta_binomial_PGLMM_Monsoon.R
├── 08_Bayesian_Assurance_Analysis_Monsoon_body_mass.R
├── Zach_data.csv
├── min20Fixed516.nwk
├── README.md
└── Zero_inflated_BetaBinom_PGLMM_LHT_Allspecies/
    ├── fit_all_three_NeoplasiaCases_zi_log_trials.rds
    ├── fit_all_three_MalignancyCases_zi_log_trials.rds
    ├── dat_lht.rds
    ├── A.rds
    ├── BetaBinom_ZI_PGLMM_LHT_summary.csv
    └── Bayesian_assurance_body_mass_all_three_model.csv
```

The output directory is created automatically by the main model script.

## Required input files

The main model script requires two input files.

### 1. Zach_data.csv

This file contains the species level cancer and life history data.

The script expects the following columns:

```text
Species
NeoplasiaWithDenominators
Malignant
Necropsies
max_longevity.months.
adult_weight.g.
Gestation.months.
```

These columns are renamed internally as:

```text
NeoplasiaCases
MalignancyCases
Trials
max_longevity_months
adult_weight_g
gestation_months
```

### 2. min20Fixed516.nwk

This is the phylogenetic tree used to construct the Brownian motion covariance matrix.

Species names in the tree and the data are cleaned by trimming spaces and replacing spaces with underscores. Only species present in both the data and the tree are retained.

## Required R packages

The scripts require the following R packages:

```r
tidyverse
ape
phytools
brms
rstan
MASS
```

On Monsoon, the following library path may need to be exported before running the scripts:

```bash
export R_LIBS_USER=/scratch/vn229/Rlibs
```

A quick package check can be run with:

```bash
Rscript -e 'packages <- c("brms","rstan","tidyverse","ape","phytools","MASS"); print(sapply(packages, requireNamespace, quietly=TRUE))'
```

All packages should return `TRUE`.

## Main model

The main script is:

```text
03_zero_inflated_beta_binomial_PGLMM_Monsoon.R
```

This script fits zero inflated beta binomial PGLMMs for neoplasia and malignancy.

### Model structure

The response is modeled as cancer cases out of necropsy trials:

```r
cases | trials(Trials)
```

The model uses a zero inflated beta binomial family:

```r
zero_inflated_beta_binomial()
```

The general model structure is:

```r
cases | trials(Trials) ~ predictors + (1 | gr(Species, cov = A_model))
zi ~ log_trials_s
```

where:

```text
cases
```

is the number of cancer cases.

```text
Trials
```

is the number of necropsies or sampled individuals.

```text
predictors
```

are standardized life history traits.

```text
Species
```

is a phylogenetic random effect.

```text
A_model
```

is the Brownian motion phylogenetic covariance matrix.

```text
zi
```

is the zero inflation component.

```text
log_trials_s
```

is standardized log sampling effort.

## Terms explained

### Beta binomial model

A beta binomial model is used when the response is a number of successes out of a number of trials, but the variation is larger than expected under a simple binomial model.

In this project:

```text
successes = species with cancer diagnosis
trials = necropsied or sampled individuals
```

The beta binomial part allows cancer prevalence to vary more across species than expected under a simple binomial model.

### Zero inflation

Zero inflation accounts for extra zeros beyond what the beta binomial model expects.

In this dataset, many species may have zero observed cancer cases. Some zeros may reflect true low cancer occurrence, but others may reflect limited sampling. The zero inflation component helps model this process.

The default script uses:

```r
zi ~ log_trials_s
```

This means the probability of being an excess zero is allowed to depend on sampling effort.

### Phylogenetic generalized linear mixed model

Species are not statistically independent because related species share evolutionary history. The model includes a phylogenetic random effect using a covariance matrix derived from the tree.

This part of the model is:

```r
(1 | gr(Species, cov = A_model))
```

### Brownian motion covariance matrix

The covariance matrix `A_model` is built from the phylogenetic tree under a Brownian motion model of trait evolution. Closely related species have higher expected covariance than distantly related species.

### Standardized predictors

The life history predictors are log transformed and standardized:

```r
log10(adult body mass)
log10(maximum longevity)
log10(gestation length)
```

Then each transformed predictor is scaled to have mean 0 and standard deviation 1.

This means model coefficients represent the expected change in log odds of cancer occurrence for a one standard deviation increase in the predictor.

### Odds ratio

The model coefficients are on the log odds scale. The script converts each coefficient to an odds ratio using:

```r
exp(beta)
```

Interpretation:

```text
OR > 1
```

indicates a positive association with cancer occurrence.

```text
OR < 1
```

indicates a negative association with cancer occurrence.

```text
OR = 1
```

indicates no association.

### Phi

`phi` is the beta binomial precision parameter. Lower values of `phi` indicate stronger overdispersion. Higher values indicate behavior closer to a binomial model.

### Variance inflation factor

The script calculates an approximate variance inflation factor relative to a binomial model:

```r
F = 1 + (n_med - 1) / (phi + 1)
```

where `n_med` is the median number of trials.

A larger value indicates that the data are more variable than expected under a simple binomial model.

### Overdispersion statistic c

The script also calculates an approximate overdispersion statistic based on Pearson residuals:

```r
c = sum(residuals^2) / residual degrees of freedom
```

Values much larger than 1 suggest residual overdispersion.

### Rhat

`Rhat` is a convergence diagnostic for Bayesian MCMC chains. Values close to 1 indicate good mixing across chains. Values above 1.01 may suggest convergence problems.

## Models fitted by the main script

The main script fits four predictor sets for each response.

```text
mass_only
longevity_only
gestation_only
all_three
```

The two responses are:

```text
NeoplasiaCases
MalignancyCases
```

So the full script fits eight models:

```text
4 predictor sets × 2 responses = 8 models
```

## Main output files

The main model script saves fitted model objects as `.rds` files. These can be reloaded later for summaries, plotting, or simulation.

Examples:

```text
fit_all_three_NeoplasiaCases_zi_log_trials.rds
fit_all_three_MalignancyCases_zi_log_trials.rds
```

The script also saves:

```text
dat_lht.rds
```

This contains the cleaned and transformed life history data.

```text
A.rds
```

This contains the phylogenetic covariance matrix.

```text
BetaBinom_ZI_PGLMM_LHT_summary.csv
```

This contains the summarized model results, including posterior coefficient summaries, odds ratios, Rhat, variance inflation factor, and number of species.

## Bayesian assurance analysis

The assurance script is:

```text
08_Bayesian_Assurance_Analysis_Monsoon_body_mass.R
```

Bayesian assurance is a simulation based approach for estimating the probability that a study design will provide strong evidence for an effect, assuming that effect is real.

In this project, the assurance analysis asks:

```text
If the true body mass effect had a given size, how often would this model detect it with high posterior probability?
```

The current script tests body mass effects while keeping the full model structure:

```r
pred_vars <- c("mass_s", "longevity_s", "gestation_s")
target_predictor <- "mass_s"
```

This means the analysis tests whether the current data structure can detect a body mass effect while controlling for longevity and gestation.

## Assurance effect sizes

The current body mass assurance script uses:

```r
effect_grid_positive <- c(0.1, 0.3, 0.5, 0.75)
```

These are assumed true beta values on the log odds scale.

The corresponding odds ratios are:

```r
exp(beta)
```

Approximate values:

```text
beta = 0.10  OR = 1.11
beta = 0.30  OR = 1.35
beta = 0.50  OR = 1.65
beta = 0.75  OR = 2.12
```

## What the assurance script does

For each assumed body mass effect size, the script performs the following steps.

### Step 1. Load the real fitted model

The script loads the fitted all three predictor neoplasia model:

```r
fit_all_three_NeoplasiaCases_zi_log_trials.rds
```

This model provides realistic nuisance parameter values, including the intercept, beta binomial precision, zero inflation parameters, and phylogenetic random effect variation.

### Step 2. Load the cleaned data and phylogenetic covariance matrix

The script loads:

```text
dat_lht.rds
A.rds
```

Then it recreates the model data frame and aligns species with the covariance matrix.

### Step 3. Define the assumed effect size

For each effect size in the grid, the target body mass coefficient is set to that value.

Non target predictors are set to zero during simulation to isolate the power to detect the target predictor.

### Step 4. Simulate phylogenetic random effects

The script simulates species level phylogenetic random effects from a multivariate normal distribution using the phylogenetic covariance matrix.

### Step 5. Simulate cancer counts

The script simulates zero inflated beta binomial cancer counts using:

```text
Trials
mu
phi
zi probability
```

### Step 6. Refit the same model

For every simulated dataset, the script refits the same type of brms model:

```r
zero_inflated_beta_binomial()
```

with the phylogenetic random effect and zero inflation structure.

### Step 7. Check whether the effect was detected

The script extracts posterior draws for the target body mass coefficient.

For a positive expected effect, it calculates:

```r
mean(beta_draws > 0)
```

A simulation is counted as successful if:

```r
posterior probability > 0.95
```

### Step 8. Calculate assurance

Assurance is calculated as:

```r
number of successful simulations / number of valid simulations
```

For example, if 80 out of 100 simulated datasets successfully detect the effect, the assurance is:

```text
0.80
```

This means that under the assumed effect size and current data structure, the model detects the effect in 80 percent of simulations.

## Assurance output file

The assurance script saves:

```text
Bayesian_assurance_body_mass_all_three_model.csv
```

This file includes:

```text
target_predictor
direction
assumed_beta
assumed_OR
n_sim_requested
n_sim_valid
n_success
assurance
posterior_prob_cutoff
median_estimated_beta
median_posterior_prob
```

## How to interpret assurance results

### target_predictor

The predictor being tested.

For the current script:

```text
mass_s
```

### assumed_beta

The true simulated effect size on the log odds scale.

Larger values represent stronger assumed effects.

### assumed_OR

The odds ratio corresponding to `assumed_beta`.

```r
assumed_OR = exp(assumed_beta)
```

### n_sim_requested

The number of simulations requested for that effect size.

The current script uses:

```r
n_sim_per_effect = 100
```

### n_sim_valid

The number of simulations that successfully finished and produced a fitted model.

If this is much lower than `n_sim_requested`, some models failed and the result should be interpreted with caution.

### n_success

The number of valid simulations where the posterior probability for the effect direction exceeded the cutoff.

For body mass, success means:

```r
Pr(beta_mass > 0) > 0.95
```

### assurance

The estimated probability of successful detection.

General interpretation:

```text
assurance near 0.20
```

means weak detection probability.

```text
assurance near 0.50
```

means moderate detection probability.

```text
assurance near 0.80 or higher
```

means strong detection probability.

These thresholds are practical guidelines, not strict rules.

### median_estimated_beta

The median estimated body mass coefficient across simulations.

This can be compared with `assumed_beta` to assess whether the model tends to recover the simulated effect size.

### median_posterior_prob

The median posterior probability in the expected direction across simulations.

For body mass, this is the median value of:

```r
Pr(beta_mass > 0)
```

## Running the main model on Monsoon

Example Slurm script:

```bash
#!/bin/bash

#SBATCH --job-name=ZiBBPGLMM
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --output=./ZiBBPGLMM.out
#SBATCH --error=./ZiBBPGLMM.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=vn229@nau.edu

export R_LIBS_USER=/scratch/vn229/Rlibs

/packages/spack/1.1.0/var/spack/environments/r-2601/.spack-env/view/bin/Rscript 03_zero_inflated_beta_binomial_PGLMM_Monsoon.R
```

Submit with:

```bash
sbatch run_ZiBBPGLMM.sh
```

## Running the assurance analysis on Monsoon

Example Slurm script:

```bash
#!/bin/bash

#SBATCH --job-name=ZiBB_Assurance
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --output=./ZiBB_Assurance.out
#SBATCH --error=./ZiBB_Assurance.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=vn229@nau.edu

export R_LIBS_USER=/scratch/vn229/Rlibs

/packages/spack/1.1.0/var/spack/environments/r-2601/.spack-env/view/bin/Rscript 08_Bayesian_Assurance_Analysis_Monsoon_body_mass.R
```

Submit with:

```bash
sbatch run_assurance.sh
```

## Changing the assurance script to test all traits

The current script runs body mass only by default.

To also run longevity and gestation, change:

```r
run_longevity_gestation <- FALSE
```

to:

```r
run_longevity_gestation <- TRUE
```

This will run assurance for:

```text
mass_s
longevity_s
gestation_s
```

Body mass and longevity use positive effect sizes. Gestation uses negative effect sizes.

## Computational notes

The assurance analysis is computationally expensive because it repeatedly refits full Bayesian models.

For the current body mass setup:

```text
4 effect sizes × 100 simulations = 400 brms model refits
```

Each refit uses MCMC sampling. This is why the script can take a long time.

For larger runs, it is better to split the analysis using Slurm job arrays, for example one job per effect size.

## Recommended workflow

1. Place the data file, tree file, and R scripts in the project directory.
2. Run the main ZiBBPGLMM script.
3. Confirm that the output folder contains `.rds` model objects, `dat_lht.rds`, and `A.rds`.
4. Run the assurance script for body mass.
5. Check the assurance CSV.
6. If needed, modify the script to run longevity and gestation.
7. For large runs, split the assurance analysis by effect size using Slurm job arrays.

## Troubleshooting

### R packages not found in Slurm

Make sure the Slurm script includes:

```bash
export R_LIBS_USER=/scratch/vn229/Rlibs
```

### Module `r` not found

Do not use:

```bash
module load r
```

Instead, use the full Rscript path:

```bash
/packages/spack/1.1.0/var/spack/environments/r-2601/.spack-env/view/bin/Rscript
```

### Assurance script cannot find `A.rds`

Make sure the main model script includes:

```r
saveRDS(A, file.path(out_dir, "A.rds"))
```

and that the main script finished successfully.

### Assurance script cannot find the fitted model

Check that the output folder contains:

```text
fit_all_three_NeoplasiaCases_zi_log_trials.rds
```

If you want to run assurance for malignancy instead, change the loaded model and response name in the assurance script.

### Some simulations fail

Some Bayesian model refits may fail. The script records these as missing and reports `n_sim_valid`.

If many simulations fail, check the error file and consider testing with fewer simulations or a smaller effect grid.

## Citation and reuse

When using this code, describe the method as a zero inflated beta binomial phylogenetic generalized linear mixed model implemented in `brms`, followed by a simulation based Bayesian assurance analysis.

A concise description is:

```text
We fit zero inflated beta binomial phylogenetic generalized linear mixed models to cancer counts with necropsy denominators, using a Brownian motion phylogenetic covariance matrix to account for shared ancestry among species. We then performed a Bayesian assurance analysis by simulating datasets under assumed effect sizes and refitting the same model to estimate the probability of detecting each effect with high posterior probability.
```

## Notes

This code was developed for cross species cancer data where sampling effort varies among species and many species have zero observed cancer cases. The modeling approach is intended to account for several important features of these data, including trial based responses, overdispersion, zero inflation, and phylogenetic non independence.
