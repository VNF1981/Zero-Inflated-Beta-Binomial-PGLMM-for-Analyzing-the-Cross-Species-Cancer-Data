######################################################################################################################################
# 
# This script fits Zero-Inflated Beta-Binomial Phylogenetic Generalized Linear Mixed Models (PGLMM)
# to study cancer prevalence (neoplasia and malignancy) across species.
# The goal is to model cancer counts while accounting for:
# - phylogenetic relatedness
# - overdispersion
# - excess zeros
# - variation in sampling effort
# - life history traits
# 
######################################################################################################################################

rm(list = ls(all = TRUE))   
ls()                        

# allow Stan to use multiple CPU cores
options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE)

# load libraries quietly
suppressPackageStartupMessages({
  library(tidyverse)   # data manipulation
  library(ape)         # phylogenetic tools
  library(phytools)    # extra phylogenetic methods
  library(brms)        # Bayesian modeling interface to Stan
})

###################################################################################################
# 0. Choose the zero-inflation model
# Zero inflation helps the model handle extra zeros in the data.
# "log_trials" uses sampling effort (log of necropsies) to predict zero inflation (default mode).
# "intercept" uses a single global zero-inflation level (alternative simplified mode).
###################################################################################################

zi_mode_global <- "log_trials"
# zi_mode_global <- "intercept"

###################################################################################################
# 1. Load the Compton dtaset and the phylogenetic tree.
# Clean species names, match species between the data and the tree, and prune both accordingly.
###################################################################################################

setwd("C: ... Path to the folder that includes the data and tree...your results will also be saved here ...")

dat_raw <- read.csv("Compton_data.csv", header = TRUE)   # load raw data

# clean species names (trim spaces and replace spaces with underscores)
dat_raw$Species <- trimws(dat_raw$Species)
dat_raw$Species <- gsub(" ", "_", dat_raw$Species)

# for now, keep all species
dat <- dat_raw

cat("Species rows:", nrow(dat), "\n")
cat("Unique species:", length(unique(dat$Species)), "\n")

# load the phylogenetic tree and clean tip labels
tree <- ape::read.tree("min20Fixed516.nwk")
tree$tip.label <- trimws(tree$tip.label)
tree$tip.label <- gsub(" ", "_", tree$tip.label)

# find species present in both the dataset and the tree
common_species <- intersect(tree$tip.label, dat$Species)
cat("Overlap species:", length(common_species), "\n")

if (length(common_species) == 0) {
  stop("No species match between tree and data even after cleaning")
}

# prune tree to species with data
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, common_species))

# keep only matched species in the data
dat <- dat %>% dplyr::filter(Species %in% tree_pruned$tip.label)

# reorder data rows to match the tree order
dat <- dat[match(tree_pruned$tip.label, dat$Species), ]

# sanity check: data and tree must be perfectly aligned
stopifnot(all(dat$Species == tree_pruned$tip.label))

#view(dat)

###################################################################################################
# 2. Rename important columns and check basic properties.
# These columns are used directly by the beta-binomial model.
###################################################################################################

dat <- dat %>%
  dplyr::rename(
    NeoplasiaCases        = NeoplasiaWithDenominators,   # number of neoplasia cases
    MalignancyCases       = Malignant,                   # number of malignancy cases
    Trials                = Necropsies,                  # number of necropsies (sample size)
    max_longevity_months  = max_longevity.months.,
    adult_weight_g        = adult_weight.g.,
    gestation_months      = Gestation.months.
  )

###################################################################################################
# 3. Build the phylogenetic correlation matrix using Brownian motion.
# This matrix tells the model how related species are.
###################################################################################################

A <- ape::vcv(tree_pruned, corr = TRUE)   # Brownian-motion correlation matrix
A <- A[dat$Species, dat$Species]          # reorder to match the data

# check again that matrix and data match perfectly
stopifnot(identical(rownames(A), dat$Species))
stopifnot(identical(colnames(A), dat$Species))

###################################################################################################
# 4. Prepare predictors: log-transform and standardize life history traits.
# Also compute log_trials_s as a zero-inflation predictor.
###################################################################################################

dat_lht <- dat %>%
  dplyr::mutate(
    # clean non-positive values
    adult_weight_g_clean       = ifelse(adult_weight_g > 0, adult_weight_g, NA_real_),
    max_longevity_months_clean = ifelse(max_longevity_months > 0, max_longevity_months, NA_real_),
    gestation_months_clean     = ifelse(gestation_months > 0, gestation_months, NA_real_),
    
    # log-transform continuous predictors
    log_mass       = log10(adult_weight_g_clean),
    log_longevity  = log10(max_longevity_months_clean),
    log_gestation  = log10(gestation_months_clean),
    
    # standardize predictors (mean 0, SD 1)
    mass_s         = as.numeric(scale(log_mass)),
    longevity_s    = as.numeric(scale(log_longevity)),
    gestation_s    = as.numeric(scale(log_gestation)),
    
    # log of necropsies, standardized (used only for zero inflation)
    log_trials     = log(Trials),
    log_trials_s   = as.numeric(scale(log_trials))
  ) %>%
  dplyr::select(
    Species, NeoplasiaCases, MalignancyCases, Trials,
    mass_s, longevity_s, gestation_s, log_trials_s
  )

###################################################################################################
# 5. Set prior distributions for the Bayesian model.
# Priors help the model remain stable and regularized.
###################################################################################################

base_priors <- c(
  set_prior("normal(0, 1.5)", class = "Intercept"),
  set_prior("normal(0, 1)",   class = "b"),
  set_prior("normal(0, 1)",   class = "sd",  group = "Species"),
  set_prior("normal(0, 1)",   class = "phi"),       # beta-binomial dispersion
  
  set_prior("normal(0, 1.5)", class = "Intercept", dpar = "zi"),     # baseline zero-inflation level
  set_prior("normal(0, 1)",   class = "b",         dpar = "zi")      # effect of predictors on zero-inflation
)

###################################################################################################
# 6. Create an output directory for saving results.
###################################################################################################

# Change the output file name to mammals, birds, etc if you are running the script for each subsets to avoid overwriting the results

out_dir <- file.path(getwd(), "Zero_inflated_BetaBinom_PGLMM_LHT_AllSpecies") 
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

results_list <- list()

###################################################################################################
# 7. Function to fit one model (one response and one predictor set).
# This function:
# - prepares the data
# - builds the model formulas
# - fits the Bayesian model using brms
# - extracts posterior summaries
# - computes diagnostics (VIF and binomial overdispersion c)
# - writes results to file
###################################################################################################

fit_lht_model <- function(response_name,
                          pred_vars,
                          model_label,
                          dat_lht,
                          A,
                          priors,
                          out_dir,
                          zi_mode = c("log_trials", "intercept")) {
  
  zi_mode <- match.arg(zi_mode)
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  message("Fitting model ", model_label,
          " for response ", response_name,
          " with zero inflation mode: ", zi_mode)
  
  # select columns requred for this model
  cols_needed <- c("Species", "Trials", "log_trials_s", response_name, pred_vars)
  
  df <- dat_lht %>%
    dplyr::select(all_of(cols_needed)) %>%
    tidyr::drop_na()   # drop rows with missing value
  
  if (nrow(df) < 20) {
    warning("Skipping model (not enough species)")
    return(NULL)
  }
  
  # reorder matrix and data to match the species used
  species_use <- intersect(rownames(A), df$Species)
  df <- df[df$Species %in% species_use, ]
  df <- df[match(species_use, df$Species), ]
  df$Species <- factor(df$Species, levels = species_use)
  A_model <- A[species_use, species_use]
  
  # build formulas for brms
  fixed_term <- paste(pred_vars, collapse = " + ")
  
  form_str <- paste(
    "cases | trials(Trials) ~",
    fixed_term,
    "+ (1 | gr(Species, cov = A_model))"
  )
  
  zi_form_str <- ifelse(zi_mode == "log_trials",
                        "zi ~ log_trials_s",
                        "zi ~ 1")
  
  bf_mod <- brms::bf(
    as.formula(form_str),
    as.formula(zi_form_str),
    family = zero_inflated_beta_binomial()
  )
  
  # rename response variable to "cases" for brms
  df_model <- df[, c("Species", "Trials", "log_trials_s", pred_vars, response_name)]
  names(df_model)[ncol(df_model)] <- "cases"
  
  # fit the Bayesian model
  fit <- brms::brm(
    formula = bf_mod,
    data    = df_model,
    data2   = list(A_model = A_model),
    prior   = priors,
    chains  = 4,
    iter    = 4000,
    warmup  = 2000,
    control = list(adapt_delta = 0.99, max_treedepth = 14),
    seed    = 1234
  )
  
  summ <- summary(fit)
  post <- as_draws_df(fit)
  
  ############################################
  # 1. Extract odds ratios for each predictor
  ############################################
  rows <- list()
  
  for (pv in pred_vars) {
    par_name <- paste0("b_", pv)
    if (!par_name %in% names(post)) next
    
    beta <- post[[par_name]]
    OR   <- exp(beta)
    
    beta_q <- quantile(beta, probs = c(0.025, 0.5, 0.975))
    OR_q   <- quantile(OR,   probs = c(0.025, 0.5, 0.975))
    
    rows[[pv]] <- tibble::tibble(
      model_label  = model_label,
      response     = response_name,
      predictor    = pv,
      n_species    = nrow(df_model),
      beta_l95     = beta_q[1],
      beta_med     = beta_q[2],
      beta_u95     = beta_q[3],
      OR_l95       = OR_q[1],
      OR_med       = OR_q[2],
      OR_u95       = OR_q[3]
    )
  }
  
  ############################################
  # 2. Rhat convergence diagnostic
  ############################################
  rhat_vals <- summ$fixed[, "Rhat"]
  if (!is.null(summ$random$Species)) {
    rhat_vals <- c(rhat_vals, summ$random$Species[, "Rhat"])
  }
  if (!is.null(summ$spec_pars)) {
    rhat_vals <- c(rhat_vals, summ$spec_pars[, "Rhat"])
  }
  max_rhat <- max(rhat_vals, na.rm = TRUE)
  
  ############################################
  # 3. Custom variance inflation factor (VIF)
  # VIF tells us how much larger the variance is 
  # under the beta-binomial model compared with a
  # standard binomial model.
  ############################################
  if (!"phi" %in% names(post)) {
    vif_q <- c(NA, NA, NA)
    n_med <- median(df_model$Trials)
  } else {
    phi_draws <- post$phi
    n_med     <- median(df_model$Trials)
    vif_draws <- 1 + (n_med - 1) / (phi_draws + 1)
    vif_q     <- quantile(vif_draws, probs = c(0.025, 0.5, 0.975))
  }
  
  res_df <- dplyr::bind_rows(rows) %>%
    dplyr::mutate(
      max_rhat = max_rhat,
      n_median = n_med,
      VIF_l95  = vif_q[1],
      VIF_med  = vif_q[2],
      VIF_u95  = vif_q[3]
    )
  
  ############################################
  # 4. Formal binomial overdisperssion statistic (c)
  # This shows how badly a *standard binomial* model
  # would fit the data.
  ############################################
  r_vec <- as.numeric(residuals(fit, type = "pearson"))
  
  p_fixed <- nrow(brms::fixef(fit))   # number of fixed parameters
  p_extra <- 3                        # random variance + phi + zi
                                      # Note: some scholars also include the number of random-effect parameters,
                                      # both are correct but that approach is more conservative.
                                      # The impact of using either approach on the results is negligible
  p_use   <- p_fixed + p_extra
  
  n_obs <- nrow(df_model)
  df_res <- n_obs - p_use
  
  c_hat <- sum(r_vec^2) / df_res
  
  ############################################
  # 5. Zero-inflation summary (zi)
  ############################################
  if (zi_mode == "intercept") {
    zi_q <- if ("zi" %in% names(post))
      quantile(post$zi, probs = c(0.025, 0.5, 0.975))
    else
      c(NA, NA, NA)
  } else {
    zi_q <- c(NA, NA, NA)   # no single global zi under log_trials mode
  }
  
  ############################################
  # Write overdispersion diagnostic to file
  ############################################
  c_file <- file.path(
    out_dir,
    paste0("standard_binomial_overdispersion_", model_label, "_", response_name, ".txt")
  )
  
  write.table(
    data.frame(
      model_label = model_label,
      response    = response_name,
      n_obs       = n_obs,
      df_resid    = df_res,
      c_hat       = c_hat
    ),
    file = c_file,
    row.names = FALSE,
    quote = FALSE
  )
  
  ############################################
  # Write text summary (priors, VIF, zi, odds ratios)
  ############################################
  txt_file <- file.path(
    out_dir,
    paste0("summary_", model_label, "_", response_name, ".txt")
  )
  
  sink(txt_file)
  cat("Model", model_label, "response", response_name, "\n\n")
  print(summ)
  
  cat("\nCustom variance inflation factor relative to binomial\n")
  cat("Median Trials (n_med):", n_med, "\n")
  cat("VIF 95 percent CrI:",
      signif(vif_q[1], 3), ";",
      signif(vif_q[2], 3), ";",
      signif(vif_q[3], 3), "\n\n")
  
  cat("Standard binomial overdispersion statistic (c)\n")
  cat("Number of observations:", n_obs, "\n")
  cat("Residual df:", df_res, "\n")
  cat("c =", signif(c_hat, 3), "\n\n")
  
  if (zi_mode == "intercept") {
    cat("Zero-inflation parameter zi (95 percent CrI)\n")
    cat("zi CrI:",
        signif(zi_q[1], 3), ";",
        signif(zi_q[2], 3), ";",
        signif(zi_q[3], 3), "\n\n")
  } else {
    cat("Zero-inflation model: zi ~ log_trials_s\n")
    cat("No single global zi parameter under this mode.\n\n")
  }
  
  cat("Posterior odds ratios and VIF summary\n")
  print(res_df %>% dplyr::select(predictor, OR_l95, OR_med, OR_u95, VIF_med))
  sink()
  
  return(res_df)
}

###################################################################################################
# 8. Define model settings and fit all models.
# Here we run models for:
# - mass only
# - longevity only
# - gestation only
# - all three predictors
# And for two responses separately:
# - NeoplasiaCases
# - MalignancyCases
###################################################################################################

model_specs <- list(
  list(label = "mass_only",      preds = c("mass_s")),
  list(label = "longevity_only", preds = c("longevity_s")),
  list(label = "gestation_only", preds = c("gestation_s")),
  list(label = "all_three",      preds = c("mass_s", "longevity_s", "gestation_s"))
)

responses <- c("NeoplasiaCases", "MalignancyCases")

for (ms in model_specs) {
  for (resp in responses) {
    res_df <- fit_lht_model(
      response_name = resp,
      pred_vars     = ms$preds,
      model_label   = ms$label,
      dat_lht       = dat_lht,
      A             = A,
      priors        = base_priors,
      out_dir       = out_dir,
      zi_mode       = zi_mode_global
    )
    if (!is.null(res_df)) {
      key <- paste(ms$label, resp, sep = "_")
      results_list[[key]] <- res_df
    }
  }
}

###################################################################################################
# 9. Combine all model summaries and save them to CSV.
###################################################################################################

if (length(results_list) == 0) {
  warning("No models were fitted. Check data coverage and thresholds.")
} else {
  results_all <- dplyr::bind_rows(results_list)
  print(results_all)
  
  out_csv <- file.path(out_dir, "BetaBinom_ZI_PGLMM_LHT_summary.csv")
  write.csv(results_all, out_csv, row.names = FALSE)
}
