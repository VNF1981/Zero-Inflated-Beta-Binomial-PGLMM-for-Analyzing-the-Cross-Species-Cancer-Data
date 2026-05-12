# Bayesian Assurance Analysis for a ZiBBPGLMM

## Overview

This repository contains an R script for running a Bayesian assurance analysis for a fitted ZiBBPGLMM. The general idea for this simulation-based analysis is adapted from [bayesassurance](https://journal.r-project.org/articles/RJ-2023-066/) R package. 

### *** I will soon expand this framework to Gompertz models as well *** 

The goal of this analysis is to evaluate whether the current dataset and model structure are capable of detecting a body mass effect of different assumed sizes.

In simple terms, the script asks:
#### If the true body mass effect were small, moderate, or large, how often would this model detect it using the data structure we currently have?

The current version is set up to test the body mass effect while keeping longevity and gestation in the model as covariates.

## What is Bayesian assurance?
Bayesian assurance is similar in spirit to power analysis in the frequentist statistics, but it is framed in a Bayesian way.

A traditional power analysis usually asks:
- If an effect exists, how often would a statistical test reject the null hypothesis?

Bayesian assurance asks:
- If an effect of a given size exists, how often would the Bayesian model produce strong posterior support for that effect?

In this script, strong posterior support is defined as:
```r
Pr(beta_mass > 0) > 0.95
```

This means that, for a given model fit, more than 95% of the posterior samples for the body mass coefficient must be greater than zero. In other words, the model must show strong posterior support that the body mass effect is positive. We say positive here because body mass is expected to have a positive association with cancer occurrence. Obviously, this can be negative for traits like gestation length, such as:
```r
Pr(beta_gestation < 0) > 0.95
```

The assurance value is the proportion of simulations in which the model successfully detects the assumed effect.
For example:
```text
100 simulated datasets
80 have Pr(beta_mass > 0) > 0.95
Assurance = 80 / 100 = 0.80
```

This means the model successfully detected the assumed body mass effect in 80% of the simulated datasets.

## Main question addressed by this script

The script evaluates the following question:
- Given the current species, necropsy counts, predictor values, phylogenetic structure, and model specification, how reliably can we detect different assumed true body mass effects?

The script does not manually alter the observed cancer counts. Instead, it changes the assumed true beta value for body mass. That beta changes the expected cancer probability for each species. The script then simulates new cancer counts from those expected probabilities under the zero inflated beta binomial model.

### Because this analysis repeatedly refits full Bayesian models to simulated datasets, it is computationally time consuming. The full workflow is designed to run on a high performance computing cluster, such as Agave, Monsoon, or another Slurm based system. Running the full analysis on a local machine is not a good idea because it can take many hours and may be limited by CPU availability, memory, or thermal throttling.

## Required input files
#### This script is designed to run after the main model has already been fitted.
#### It requires output files from the fitted model analysis. The required files are:
```text
Zero_inflated_BetaBinom_PGLMM_LHT_Allspecies/
├── fit_all_three_NeoplasiaCases_zi_log_trials.rds
├── dat_lht.rds
└── A.rds
```

### fit_all_three_NeoplasiaCases_zi_log_trials.rds
This is the fitted Bayesian model object from the real data.
The assurance script uses this object to extract realistic nuisance parameters, including:
```text
intercept
phi
zero inflation parameters
phylogenetic random effect standard deviation
```
These values are used to simulate realistic datasets.

### dat_lht.rds
This file contains the cleaned and transformed dataset used by the model. It includes:
```text
Species
NeoplasiaCases
MalignancyCases
Trials
mass_s
longevity_s
gestation_s
log_trials_s
```

### A.rds
This file contains the phylogenetic covariance matrix aligned to the species in the dataset.
It is used to simulate phylogenetically structured species effects and to refit the model with the same phylogenetic structure.

## Required R packages
The assurance script requires:
```r
brms
rstan
tidyverse
MASS
```
It also uses objects produced by the main model workflow, which requires packages such as:
```r
ape
phytools
```

A quick check to confirm that the required packages are available in your active R or conda environment can be run with:
```bash
Rscript -e 'packages <- c("brms","rstan","tidyverse","ape","phytools","MASS"); print(sapply(packages, requireNamespace, quietly=TRUE))'
```
*** All packages should return `TRUE`.

## Script

The main script is:
```text
Bayesian_Assurance_Analysis.R
```

This script is currently configured to run assurance for body mass only.

## Effect sizes tested
The script tests a grid of assumed true body mass beta values:
```r
effect_grid_positive <- c(0.1, 0.3, 0.5, 0.7, 0.9)
```

These values are on the log odds scale.
Their approximate odds ratios (calculated as exp(x)) are:

```text
beta = 0.10    OR = 1.11
beta = 0.30    OR = 1.35
beta = 0.50    OR = 1.65
beta = 0.75    OR = 2.12
```
These effect sizes represent increasingly stronger assumed positive associations between body mass and cancer occurrence.

## Important terms

### Assumed true beta
The assumed true beta is the effect size used to simulate new datasets.
For example:
```r
beta_mass = 0.5
```
means that the simulation assumes a true positive body mass effect of 0.5 on the log odds scale.

### Target predictor
The target predictor is the predictor being evaluated in the assurance analysis.
In the current script:
```r
target_predictor = "mass_s"
```
This means the assurance analysis is focused on body mass.

### Non target predictors
The model still includes all three predictors:
```r
mass_s
longevity_s
gestation_s
```
However, during simulation, only the target predictor is assigned the assumed effect size. The non target predictors are set to zero in the data generating process. This isolates the model's ability to detect the body mass effect.

### Trials
`Trials` is the number of sampled individuals or necropsies for each species. The script keeps the observed trial counts fixed. This preserves the real sampling structure of the dataset.

### Mu
`mu` is the expected cancer probability for each species. The script calculates `mu` from the linear predictor:
```r
eta_mu = intercept + beta_mass * mass_s + phylogenetic effect
mu = plogis(eta_mu)
```

### Phi
`phi` is the beta binomial overdispersion parameter. It controls the amount of overdispersion in the simulated cancer counts.

- Lower `phi` means more overdispersion.
- Higher `phi` means the model behaves more like a binomial model.

### Zero inflation probability
The zero inflation probability is the probability that a species is assigned to the extra zero process.
The default zero inflation model is:
```r
zi ~ log_trials_s
```
This means the probability of an excess zero depends on sampling effort.

### Phylogenetic random effect
The script simulates a phylogenetic random effect using the covariance matrix `A_model`.
This preserves the idea that related species are not statistically independent. 

### NOTE: In this script, the phylogenetic covariance matrix is based on a Brownian motion model. If a different evolutionary framework is desired, such as Pagel's lambda, lambda should first be estimated and then applied to the off-diagonal elements of the phylogenetic variance-covariance matrix. The resulting VCV matrix should then be used in the simulations and model refitting steps.

### Posterior probability cutoff
The current cutoff is:
```r
posterior_prob_cutoff = 0.95
```

For body mass, a simulation is counted as successful if:
```r
Pr(beta_mass > 0) > 0.95
```

### Assurance
Assurance is the proportion of valid simulations that successfully detect the assumed effect.
```r
assurance = n_success / n_sim_valid
```

## Step by step explanation of what the script does
### Step 1. Load packages
The script loads the required packages:

```r
library(MASS)
library(brms)
library(rstan)
library(tidyverse)
```

`MASS` is used to simulate phylogenetic random effects from a multivariate normal distribution.
`brms` and `rstan` are used to refit the Bayesian model.
`tidyverse` is used for data manipulation and output formatting.

### Step 2. 
The output directory is:
```r
Zero_inflated_BetaBinom_PGLMM_LHT_Allspecies
```

### Step 3. Load the fitted real model
The script loads the fitted model from the real data:
```r
fit_real <- readRDS(
  file.path(out_dir, "fit_all_three_NeoplasiaCases_zi_log_trials.rds")
)
```
This model is not being refitted at this step. It is used to extract realistic parameter values for simulation.

### Step 4. Load the cleaned data and phylogenetic covariance matrix
The script loads:
```r
dat_lht <- readRDS(file.path(out_dir, "dat_lht.rds"))
A <- readRDS(file.path(out_dir, "A.rds"))
```
These objects provide the species level data and phylogenetic covariance structure.

### Step 5. Recreate the model dataframe
The script selects the columns needed for the assurance analysis:
```r
Species
Trials
log_trials_s
mass_s
longevity_s
gestation_s
NeoplasiaCases
```
Then it renames the response column to:
```r
cases
```
The observed case counts are needed to build the dataframe, but the script later replaces them with simulated counts during each simulation.

### Step 6. Align the dataframe and covariance matrix
The script makes sure that the species order in the dataframe matches the species order in the covariance matrix. This is essential because the phylogenetic covariance matrix must correspond exactly to the order of species in the model data.

### Step 7. Define the effect size grid
The script defines the assumed true body mass effects:
```r
effect_grid_positive <- c(0.1, 0.3, 0.5, 0.75)
```
The script will then run simulations separately for each value.

### Step 8. Extract nuisance parameters from the fitted real model
The script extracts posterior medians from `fit_real`.
These include:
```text
intercept_real
phi_real
zi_intercept_real
zi_log_trials_real
sd_species_real
```
These are called nuisance parameters because they are not the main target of the assurance analysis, but they are needed to simulate realistic data.

### Step 9. Start the loop over assumed effect sizes
For each assumed beta value, the script runs a set number of simulations. The current default is:
```r
n_sim_per_effect = 100
```
With five effect sizes, this gives:
```text
5 × 100 = 500 simulated model refits
```

### Step 10. Build the expected cancer probability
For each simulation, the script starts with the real model intercept:
```r
eta_mu <- intercept_real
```
Then it adds the assumed body mass effect:
```r
eta_mu <- eta_mu + beta_mass * mass_s
```
Then it adds a simulated phylogenetic random effect. Finally, the script converts the linear predictor (i.e., logit scale) to probability scale:
```r
mu <- plogis(eta_mu)
```

### Step 11. Build the zero inflation probability
If the zero inflation mode is:
```r
zi_mode = "log_trials"
```
then the script calculates:
```r
eta_zi <- zi_intercept_real + zi_log_trials_real * log_trials_s
zi_prob <- plogis(eta_zi)
```
This gives the probability that each species belongs to the excess zero process.

### Step 12. Simulate new cancer counts
The script simulates new cancer counts using:
```r
simulate_zibb_counts(
  trials = Trials,
  mu = mu,
  phi = phi_real,
  zi_prob = zi_prob
)
```
This step creates simulated case counts that are consistent with:
```text
the observed number of trials
the assumed true body mass effect (e.g., 0.1, 0.3, 0.5, etc)
the fitted overdispersion level
the fitted zero inflation structure
the fitted phylogenetic variance
```
The observed cancer counts are not manually changed. They are replaced by newly simulated counts in each simulation.

### Step 13. Refit the same Bayesian model
For each simulated dataset, the script refits the same model structure:
```r
cases | trials(Trials) ~ mass_s + longevity_s + gestation_s + (1 | gr(Species, cov = A_model))
zi ~ log_trials_s
```
The model is refitted using `brms`.

### Step 14. Extract the posterior body mass effect
After each refit, the script extracts posterior draws for:
```r
b_mass_s
```

### Step 15. Calculate posterior support in the expected direction
Because body mass is expected to have a positive effect, the script calculates:
```r
mean(beta_draws > 0)
```
This is the posterior probability that the body mass effect is positive.

### Step 16. Classify the simulation as successful or not
A simulation is counted as successful if:
```r
Pr(beta_mass > 0) > 0.95
```
Otherwise, it is counted as unsuccessful. If a model fails to fit, that simulation is recorded as missing and excluded from the assurance calculation.

### Step 17. Summarize assurance for each effect size
For each assumed beta value, the script reports:
```text
number of requested simulations
number of valid simulations
number of successful detections
assurance
median estimated beta
median posterior probability
```

### Step 18. Save the results
The script saves the final assurance table as:
```text
Bayesian_assurance_body_mass_all_three_model.csv
```

## Output columns and interpretation

### target_predictor

The predictor being tested. For the current script:
```text
mass_s
```

### direction
The expected direction of the effect. For body mass:
```text
positive
```

### assumed_beta
The assumed true beta value used to simulate datasets. This is the effect size being tested.

### assumed_OR
The odds ratio corresponding to the assumed beta:
```r
assumed_OR = exp(assumed_beta)
```

### n_sim_requested
The number of simulations requested for that effect size.

### n_sim_valid
The number of simulations that successfully finished. If this number is much lower than `n_sim_requested`, the assurance estimate may be less reliable.

### n_success
The number of valid simulations where the model detected the effect with high posterior probability.

### assurance
The proportion of valid simulations that successfully detected the effect.
```r
assurance = n_success / n_sim_valid
```

### posterior_prob_cutoff
The posterior probability threshold used to define success. The default is:
```text
0.95
```

### median_estimated_beta
The median estimated beta across valid simulations.
This can be compared with `assumed_beta`.
If the median estimated beta is much smaller or larger than the assumed beta, the model may be biased under that simulation setting.

### median_posterior_prob
The median posterior probability in the expected direction across valid simulations. For body mass, this is:
```r
Pr(beta_mass > 0)
```

## How to interpret the final results
A possible output might look like this:
```text
assumed_beta    assumed_OR    assurance
0.10            1.11          0.15
0.30            1.35          0.42
0.50            1.65          0.76
0.75            2.12          0.93
```

This would mean:
```text
A very small effect is difficult to detect.
A moderate effect is detected sometimes.
A larger effect is detected more reliably.
A strong effect is detected with high reliability.
```

The key interpretation is not whether one simulation is significant. The key interpretation is how often the model detects each assumed effect size across many simulated datasets.

## Computational notes
The script is computationally expensive because it refits a full Bayesian model for each simulated dataset.
The current body mass setup uses:
```text
5 effect sizes × 100 simulations = 500 model refits
```

Each model refit uses MCMC sampling. 
For larger analyses, I need to consider splitting the run across Slurm job arrays, for example one job per effect size.

## Changing the number of simulations
To reduce runtime for testing, change:
```r
n_sim_per_effect = 100
```
to whatever number of simulations you have envisioned, depending on available computational resources.

## Changing the effect size grid
To test different assumed body mass effects, edit:
```r
effect_grid_positive <- c(0.1, 0.3, 0.5, 0.7, 0.9)
```
For example:
```r
effect_grid_positive <- c(0.1, 0.4, 0.9)
```
More effect sizes give a smoother assurance curve but increase runtime.

## Changing the target predictor
The current target predictor is body mass:
```r
target_predictor = "mass_s"
```
To test longevity, use:
```r
target_predictor = "longevity_s"
direction = "positive"
```

To test gestation, use:
```r
target_predictor = "gestation_s"
direction = "negative"
effect_grid = effect_grid_negative
```
NOTE: For traits expected to have negative associations with neoplasia or malignancy, such as gestation length, we need to use negative input beta values. For example: -c(0.1, 0.3, 0.6, 0.9). If testing all traits, interpret each assurance result separately.

### What if some simulations fail?
Some model refits may fail because Bayesian models can have convergence or sampling problems for certain simulated datasets. 
The script records failed simulations as missing and reports:
```text
n_sim_valid
```

## Summary
This script performs a Bayesian assurance analysis for body mass effects in cross species cancer data. It varies the assumed true body mass beta, simulates datasets under the fitted model structure, refits the model to each simulated dataset, and estimates how often the model detects the effect with high posterior probability. The main output is an assurance value for each assumed effect size.
