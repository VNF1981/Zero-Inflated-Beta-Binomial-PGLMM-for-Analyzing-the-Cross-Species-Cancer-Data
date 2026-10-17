# ZiBBPGLMM V2: Joint Zero-Inflated Beta-Binomial Phylogenetic Generalized Linear Mixed Models

**An R + Stan workflow for modeling cross-species cancer prevalence while accounting for varying necropsy sample sizes, excess zero counts, extra-binomial variation, and shared evolutionary ancestry.**

The workflow fits either a **single-predictor (simple)** or a **multiple-predictor** regression, using one of three phylogenetic covariance models: **Brownian motion (BM; default), Pagel's λ (lambda), or Ornstein–Uhlenbeck (OU)**. Users select the outcome, predictors, transformations, paths, and sampling settings at the beginning of a single R script. The neoplasia and malignancy outcomes are fitted separately.

## What's new in the October 2026 implementation?

This version replaces the earlier `brms` implementation with **custom Stan models fitted through `rstan`**.

| Feature | Updated implementation |
| --- | --- |
| Phylogenetic models | BM (default), Pagel's λ, or OU, chosen by the user |
| Phylogenetic inference | Species-level phylogenetic effect jointly estimated with regression coefficients, zero inflation, dispersion, and (where applicable) λ or OU α |
| Overdispersion | Estimate **ρ (rho)** with **ρ ~ Beta(1, 9)**; derive **φ = (1 − ρ)/ρ** |
| Zero-inflation prior | Intercept **Normal(−2, 1.5)** on the logit scale; standardized log-necropsy slope **Normal(0, 1)** |
| Sampling | Revised variance-coordinate parameterization and partially non-centered phylogenetic latent effects for more stable Hamiltonian Monte Carlo sampling; intended to preserve the model's likelihood and priors |
| Regression | Automatically handles one or more user-selected predictors |
| Results | CSV summaries, an Excel model summary, fitted model object, and run metadata |

**Overdispersion in brief.** In a conventional binomial model, the count variance is determined by the success probability and number of trials. The beta-binomial permits additional variation. Here, **ρ ≈ 0** approaches the binomial limit; larger ρ means greater extra-binomial variation. The previous implementation put a half-Normal prior on **φ (precision)**; the new implementation places a **Beta(1, 9)** prior directly on **ρ**, favoring modest overdispersion while allowing larger values when supported by data. These **are different priors**, not merely a change of notation. The later sampling reparameterization is a separate computational improvement.

## Repository contents

```text
ZiBBPGLMM.R                          # Main user-configurable R script
stan/
  ZiBBPGLMM_BM_main.stan              # Brownian motion
  ZiBBPGLMM_joint_lambda_main.stan    # Pagel's lambda
  ZiBBPGLMM_joint_OU_alpha_main.stan  # Ornstein–Uhlenbeck
Compton_data.csv                      # Curated cancer and life-history data
min20Fixed516.nwk                     # Phylogenetic tree
species_aliases.csv                   # Taxonomic name reconciliation
README.md
```

The CSV, tree, and aliases are example/reference inputs and must be included when publishing this repository version. For other datasets, supply equivalent files and update the paths/column names in the R script.

## Requirements

- R with packages **`ape`**, **`rstan` (≥ 2.26)**, and **`writexl`**.
- A working C++ toolchain compatible with `rstan`/Stan (**Rtools** on supported Windows installations; suitable C++ compiler on Linux).
- Enough RAM and compute time for Bayesian sampling. Runtime can be substantial, especially with large phylogenetic trees.

Install the R packages (after setting up the compiler/toolchain):

```r
install.packages(c("ape", "rstan", "writexl"))
```

Follow the platform-specific [Stan / RStan installation instructions](https://mc-stan.org/users/interfaces/rstan). On Windows ARM-based devices, confirm the availability and compatibility of the RStan compilation toolchain before attempting a full run.

## Getting started

1. Download or clone this repository. Place the R script, data, tree, species alias file, and `stan/` folder together as shown above.
2. Open **`ZiBBPGLMM.R`** in R or RStudio and edit the **USER SETTINGS** block near the top.
3. Set the **working directory** to the repository folder, or use full paths for each input and output setting.
4. Run the script, for example:

```r
setwd("C:/Users/yourname/Desktop/Github")  # Example for Windows; edit this path
source("ZiBBPGLMM.R")
```

Alternatively, from a terminal in the repository directory:

```bash
Rscript ZiBBPGLMM.R
```

The same R code is intended for desktop R and Linux/Slurm; on Slurm, request sufficient memory/CPUs and invoke the R script from a batch job. **Paths alone do not install R, RStan, or a C++ compiler.** Relative paths are resolved from R's current working directory.

## User settings and model selection

The script begins with a configuration block similar to:

```r
data_csv    <- "Compton_data.csv"
tree_nwk    <- "min20Fixed516.nwk"
aliases_csv <- "species_aliases.csv"   # NULL if no reconciliation is needed
stan_dir    <- "stan"
output_dir  <- "results"

phylo_model <- "BM"                  # "BM", "lambda", "OU"
species_col <- "Species"
trials_col  <- "Necropsies"
cases_col   <- "NeoplasiaWithDenominators"  # Or "Malignant"
predictors  <- c("adult_weight.g.")
transforms  <- c("log10")
min_trials  <- 20L

chains  <- 4L
warmup  <- 2000L
samples <- 8000L                    # Post-warmup draws per chain
seed    <- 1234L
```

### Simple regression: one predictor

```r
predictors <- c("adult_weight.g.")
transforms <- c("log10")
```

### Multiple regression: two or more predictors

```r
predictors <- c("adult_weight.g.", "residLongevity", "Gestation.months.")
transforms <- c("log10", "identity", "log10")
```

**Specify exactly one transformation per predictor, in the same order.** Available transformations in this version are:

- `"log10"`: base-10 logarithm; requires positive values.
- `"identity"`: no transformation.

After transformation, each predictor is **centered and standardized** (one standard deviation = one unit). The standardization is determined across tree-matched species *before* filtering for complete cases for the selected regression. The script also standardizes the natural logarithm of necropsy counts for the zero-inflation model.

The number of predictors determines the regression type automatically; changing the vector does **not** require a different R script or Stan model. **Each execution fits one selected response and one selected predictor set.** To analyze multiple combinations or both cancer outcomes, run it again with different settings.

### Evolutionary model choices

| Setting | Phylogenetic structure |
| --- | --- |
| `phylo_model <- "BM"` | **Default.** Brownian-motion covariance derived from the supplied tree. |
| `phylo_model <- "lambda"` | Estimates **Pagel's λ** jointly with the other model parameters, adjusting phylogenetic correlation. |
| `phylo_model <- "OU"` | Estimates an **OU α** parameter jointly with the other model parameters, controlling how evolutionary covariance decays with phylogenetic distance. |

These are **alternative structures**, not three models that are all fitted in a single execution. The main R script chooses the matching Stan file automatically. All choices require a rooted phylogenetic tree with valid branch lengths.

### Cancer outcomes and minimum sample size

- Set `cases_col <- "NeoplasiaWithDenominators"` for **neoplasia cases**, or `cases_col <- "Malignant"` for **malignant cancer cases**.
- `trials_col <- "Necropsies"` supplies the denominator (number of necropsied animals per species).
- `min_trials <- 20L` requires at least 20 necropsies per species; **the current script stops with an error** if an input row is below this threshold rather than silently dropping it. Use a suitable filtered dataset if needed.
- Counts remain counts (cases and trials); the cancer outcomes are **not** log-transformed.
- **Residual longevity** (`residLongevity`) is the curated predictor used by the updated comparative analysis. It must already be present in the input table; the R script **does not calculate residual longevity from raw longevity**. The value `-1` is treated as missing for this column.

## Model overview

For species *i*, the observed number of cases is modeled relative to the number of necropsies with a **zero-inflated beta-binomial** likelihood. Its components are:

1. **Conditional cancer prevalence:** a logit-linked regression on the selected, standardized life-history predictors, plus a phylogenetic random effect.
2. **Overdispersion:** a beta-binomial parameter **ρ**, with **ρ ~ Beta(1, 9)** and derived precision **φ = (1 − ρ)/ρ**.
3. **Zero inflation:** a separate probability of an *excess/structural zero*, with logit linked to standardized **log(necropsies)**. The intercept prior is **Normal(−2, 1.5)**. This is **not** the same as the total probability of observing zero cancer cases, since the beta-binomial component can also generate zeros.
4. **Shared ancestry:** BM, Pagel's λ, or OU species covariance, selected as described above.

The custom Stan models retain the mean-intercept prior **Normal(0, 1.5)**, life-history regression coefficient priors **Normal(0, 1)**, phylogenetic standard-deviation prior **half-Normal(0, 1)**, and zero-inflation slope prior **Normal(0, 1)**. The computational sampling changes are designed to leave these priors and the generative model unchanged.

## Updated species reconciliation

The original cross-species cancer and life-history dataset required name reconciliation against the phylogenetic tree. The updated matching procedure standardizes scientific names and applies **nine documented taxonomic aliases** (e.g., historical genus names or spelling variants), thereby recovering species previously excluded solely because their names did not match.

| Reconciliation quantity | Count |
| --- | ---: |
| Species records in original input dataset | 327 |
| Tips in the supplied phylogenetic tree | 304 |
| Matches after name standardization | 293 |
| Present in the Cancer Discovery paper | 292 |
| Additional matches recovered from nine aliases | 9 |
| **Total matched species** | **302** |
| Data species without a tree match | 25 |

**302 represents the matched data/tree species pool, not the sample size of every model.** Each individual regression excludes species lacking the required outcome, denominator, or predictor measurements. The final count depends on the selected variables. The current curated CSV should be checked against these reconciliation counts before the repository release.

The mapping is recorded in `species_aliases.csv` with columns **`data_species`** and **`tree_species`**. Keep the file alongside the data and tree to reproduce the corrected alignments. The workflow also checks for duplicate or invalid scientific names and prunes the supplied tree to the relevant species.

These taxonomic corrections are documented here for researchers using the original dataset from Zach et al.; to incorporate the additional matched species without a separate alias file, the corresponding species names must be reconciled directly in either the cancer/LHT dataset or the phylogenetic tree so that they match exactly.

## Outputs

Each execution creates a **separate results subfolder** inside `output_dir`, named for the evolutionary model, response, predictors, and transformations. Its primary outputs are:

| File | Contents |
| --- | --- |
| **`model_summary.xlsx`** | Five worksheets: **Model summary**, **Regression effects**, **Zero inflation**, **All parameters**, **Diagnostics**. Includes the regression type, selected model, estimates, 95% intervals, zero-inflation details, and sampling checks. |
| `posterior_summary.csv` | Posterior mean, median, intervals, standard deviation, effective sample size, and R-hat for main parameters. |
| `sampler_diagnostics.csv` | Compact sampling checks (divergences, tree-depth hits, BFMI, R-hat, effective sample size). |
| `analysis_species.csv` | Species retained in this specific regression. |
| `scaling_parameters.csv` | Transformations and the centering/scaling parameters used for predictors and log-necropsy effort. |
| `fit.rds` | Complete saved Stan fit for further analysis. |
| `run_settings.rds` | Selected user settings and Stan-file checksum. |

In the Excel report, **Regression effects** includes log-odds coefficients and their exponentiated odds ratios. **Zero inflation** provides the excess-zero intercept, the sampling-effort coefficient, and the estimated excess-zero probability at average standardized log-necropsy effort. The report does not estimate a separate zero-inflation coefficient for each life-history trait: **necropsy effort is the zero-inflation predictor**.

**The public script fits and summarizes the selected model. It does not automatically produce corrected integrated leave-one-species-out (LOSO) comparisons, exact refits, or the complete production manuscript tables.** Low effective sample size or other sampler warnings require additional sampling and assessment before scientific interpretation.

## Short functionality test

Before a full fit, users may temporarily select smaller values, for example:

```r
chains  <- 2L
warmup  <- 200L
samples <- 400L
```

This is a **software execution check only**, not sufficient for reliable inference. For substantive work, use adequate sampling and inspect the resulting diagnostics.

## Data source and reproducibility

This workflow extends the cross-species cancer data published by **Compton et al. (2025)** and incorporates later taxonomic reconciliation and curated life-history information. Credit the original data creators and cite the relevant dataset/publication when reusing the example data. Source: [Compton et al. (2025), linked from the original repository](https://github.com/VNF1981/Zero-Inflated-Beta-Binomial-PGLMM-for-Analyzing-the-Cross-Species-Cancer-Data).

Older implementations and historical analyses should remain in a clearly labeled archive; this README documents only the **October 2026 custom-Stan implementation**. The repository's licensing and citation information should be maintained in their own files if applicable.
