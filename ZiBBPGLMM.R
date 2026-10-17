################################################################################
# ZiBBPGLMM v2: joint zero-inflated beta-binomial phylogenetic GLMM (VNF-October 2026)
# Fit ONE simple or multiple regression using BM (default), Pagel's lambda, or OU.
# Requires R packages ape, rstan (>= 2.26), writexl and the three Stan files 
# for three evolutionary models.

rm(list=ls(all=TRUE))
ls()

# -------------------- USER SETTINGS --------------------
data_csv    <- "Compton_data.csv"       # Full or relative paths; use / on Windows
tree_nwk   <- "min20Fixed516.nwk"
aliases_csv <- "species_aliases.csv"    # Set to NULL if names already match
stan_dir    <- "stan"
output_dir  <- "results"

phylo_model <- "BM"                  # "BM" (default), "lambda", or "OU"
species_col <- "Species"
trials_col  <- "Necropsies"
cases_col   <- "NeoplasiaWithDenominators"  # Or "Malignant"

# One predictor = simple regression; two or more = multiple regression.
predictors  <- c("adult_weight.g.", "residLongevity")   # e.g. c("adult_weight.g.", "residLongevity", "Gestation.months.")
transforms  <- c("log10","identity")                   # e.g. c("log10", "identity", "log10"); one per predictor
# Note: You must specify one transformation for each predictor.
# The number and order of transformations must match the predictors.
# Example: predictors = c("adult_weight.g.", "residLongevity")
#          transforms = c("log10", "identity")
min_trials  <- 20L #minimum number of necropsies for species to be included in the analysis

chains <- 4L
warmup <- 2000L
samples <- 4000L                        # Post-warmup iterations per chain
seed <- 1234L
# -------------------------------------------------------

suppressPackageStartupMessages({
  library(ape)
  library(rstan)
})
if (packageVersion("rstan") < "2.26.0") stop("rstan >= 2.26 is required.")
if (!requireNamespace("writexl", quietly = TRUE))
  stop("Install the Excel writer first: install.packages('writexl')")
phylo_model <- match.arg(phylo_model, c("BM", "lambda", "OU"))
if (!length(predictors) || length(predictors) != length(transforms) ||
    any(!transforms %in% c("log10", "identity"))) {
  stop("Supply at least one predictor and one valid transform per predictor.")
}
regression_type <- if (length(predictors) == 1L) "Simple regression" else "Multiple regression"
message(regression_type, " (", length(predictors), " predictor(s)): ", paste(predictors, collapse = ", "))
stan_names <- c(BM = "ZiBBPGLMM_BM_main.stan",
                lambda = "ZiBBPGLMM_joint_lambda_main.stan",
                OU = "ZiBBPGLMM_joint_OU_alpha_main.stan")
stan_file <- file.path(stan_dir, unname(stan_names[phylo_model]))
needed <- c(data_csv, tree_nwk, stan_file, if (!is.null(aliases_csv)) aliases_csv)
if (any(!file.exists(needed))) stop("Missing file(s): ", paste(needed[!file.exists(needed)], collapse = ", "))

# Match species names, optionally using documented taxonomic aliases.
clean_species <- function(x) {
  x <- gsub("_+", "_", gsub("[[:space:]]+", "_", trimws(as.character(x))))
  vapply(strsplit(x, "_", fixed = TRUE), function(z) {
    if (length(z) < 2L) return(z[1])
    paste(c(paste0(toupper(substr(z[1], 1, 1)),
                   tolower(substring(z[1], 2))), tolower(z[-1])), collapse = "_")
  }, character(1))
}

dat <- read.csv(data_csv, check.names = FALSE, stringsAsFactors = FALSE)
tree <- ape::read.tree(tree_nwk)
required_cols <- c(species_col, trials_col, cases_col, predictors)
if (!all(required_cols %in% names(dat)))
  stop("Missing data column(s): ", paste(setdiff(required_cols, names(dat)), collapse = ", "))
if (anyNA(dat[[species_col]]) || anyNA(tree$tip.label) ||
    any(!nzchar(trimws(as.character(dat[[species_col]])))) ||
    any(!nzchar(trimws(tree$tip.label)))) stop("Missing or blank species names.")
dat$.species <- clean_species(dat[[species_col]])
tree$tip.label <- clean_species(tree$tip.label)
if (!is.null(aliases_csv)) {
  aliases <- read.csv(aliases_csv, stringsAsFactors = FALSE)
  if (!all(c("data_species", "tree_species") %in% names(aliases)))
    stop("Alias file needs data_species and tree_species columns.")
  keys <- clean_species(aliases$data_species)
  if (anyDuplicated(keys)) stop("Duplicate species aliases.")
  mapping <- setNames(clean_species(aliases$tree_species), keys)
  hit <- dat$.species %in% names(mapping)
  dat$.species[hit] <- unname(mapping[dat$.species[hit]])
}
if (anyDuplicated(dat$.species) || anyDuplicated(tree$tip.label))
  stop("Duplicate species after name reconciliation.")
if (!ape::is.rooted(tree) || is.null(tree$edge.length) ||
    any(!is.finite(tree$edge.length)) || any(tree$edge.length < 0))
  stop("Tree must be rooted and have valid, nonnegative branch lengths.")

# Preserve the validated global transformation/scaling policy:
# standardize on ALL tree-matched species BEFORE predictor complete-case filtering.
common <- intersect(tree$tip.label, dat$.species)
if (length(common) < 20L) stop("Too few data/tree species matches.")
dat <- dat[dat$.species %in% common, , drop = FALSE]
to_drop <- setdiff(tree$tip.label, common)
if (length(to_drop)) tree <- ape::drop.tip(tree, to_drop)

n <- as.numeric(dat[[trials_col]])
y <- as.numeric(dat[[cases_col]])
valid <- function(z) is.finite(z) & abs(z - round(z)) < 1e-8
if (any(!is.na(n) & (!valid(n) | n < 0)) ||
    any(!is.na(y) & (!valid(y) | y < 0)) ||
    any(!is.na(y) & !is.na(n) & y > n)) stop("Invalid cases or necropsy counts.")
if (any(n < min_trials, na.rm = TRUE))
  stop("Input includes species below min_trials; filter input or change min_trials.")

scaling <- data.frame(predictor = c(predictors, "log_trials"),
                      transform = c(transforms, "ln"), center = NA_real_, sd = NA_real_)
for (j in seq_along(predictors)) {
  x <- suppressWarnings(as.numeric(dat[[predictors[j]]]))
  if (predictors[j] == "residLongevity") x[is.finite(x) & x == -1] <- NA_real_
  z <- rep(NA_real_, length(x))
  good <- is.finite(x) & (transforms[j] == "identity" | x > 0)
  z[good] <- if (transforms[j] == "log10") log10(x[good]) else x[good]
  m <- mean(z, na.rm = TRUE); s <- sd(z, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) stop("Cannot standardize: ", predictors[j])
  dat[[paste0(".x", j)]] <- (z - m) / s
  scaling$center[j] <- m; scaling$sd[j] <- s
}
log_n <- log(ifelse(n > 0, n, NA_real_))
m <- mean(log_n, na.rm = TRUE); s <- sd(log_n, na.rm = TRUE)
if (!is.finite(s) || s <= 0) stop("Cannot standardize log necropsy counts.")
dat$.zi <- (log_n - m) / s
scaling$center[nrow(scaling)] <- m; scaling$sd[nrow(scaling)] <- s

x_cols <- paste0(".x", seq_along(predictors))
keep <- complete.cases(dat[, c(cases_col, trials_col, ".zi", x_cols), drop = FALSE])
df <- dat[keep, , drop = FALSE]
sp <- intersect(tree$tip.label, df$.species)
df <- df[match(sp, df$.species), , drop = FALSE]
if (nrow(df) < 20L) stop("Fewer than 20 complete-case species.")
message("Matched species: ", length(common), "; analyzed species: ", nrow(df))

stan_data <- list(N = nrow(df),
  y = as.integer(df[[cases_col]]), trials = as.integer(df[[trials_col]]),
  K_mu = length(predictors), X_mu = as.matrix(df[, x_cols, drop = FALSE]),
  K_zi = 1L, X_zi = matrix(df$.zi, ncol = 1L),
  rho_prior_a = 1, rho_prior_b = 9,
  zi_intercept_prior_mean = -2, zi_intercept_prior_sd = 1.5)

if (phylo_model %in% c("BM", "lambda")) {
  A <- cov2cor(ape::vcv(tree, corr = FALSE))[sp, sp, drop = FALSE]
  diag(A) <- 1
  if (min(eigen(A, symmetric = TRUE, only.values = TRUE)$values) <= 0)
    stop("BM correlation matrix is not positive definite.")
  if (phylo_model == "BM") {
    stan_data$L_A <- t(chol(A))
  } else {
    stan_data$A_BM <- A
  }
} else {
  height <- max(ape::node.depth.edgelength(tree)[seq_along(tree$tip.label)])
  if (!is.finite(height) || height <= 0) stop("Invalid tree height.")
  scaled_tree <- tree; scaled_tree$edge.length <- scaled_tree$edge.length / height
  shared <- ape::vcv(scaled_tree, corr = FALSE)
  stan_data$T_shared <- shared[sp, sp, drop = FALSE]
  stan_data$tip_heights <- unname(diag(shared)[sp])
  if (any(stan_data$tip_heights <= 0) ||
      min(eigen(stan_data$T_shared, symmetric = TRUE, only.values = TRUE)$values) < -1e-8)
    stop("Invalid OU shared-time matrix or tip heights.")
}

# Fit and save: no Slurm dependency; honors SLURM_CPUS_PER_TASK if provided.
cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "0"))
if (is.na(cpus) || cpus < 1L) cpus <- parallel::detectCores()
if (is.na(cpus) || cpus < 1L) cpus <- 1L
cores <- min(chains, cpus)
run_name <- paste(c(phylo_model, cases_col, paste(predictors, transforms, sep = "-")),
                  collapse = "_")
result_dir <- file.path(output_dir, gsub("[^A-Za-z0-9_-]", "_", run_name))
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
rstan::rstan_options(auto_write = FALSE)
fit <- rstan::sampling(rstan::stan_model(file = stan_file), data = stan_data,
  chains = chains, iter = warmup + samples, warmup = warmup, cores = cores,
  seed = seed,
  control = list(adapt_delta = if (phylo_model == "lambda") 0.99 else 0.999,
                 max_treedepth = 14))
saveRDS(fit, file.path(result_dir, "fit.rds"))
write.csv(data.frame(Species = df$.species), file.path(result_dir, "analysis_species.csv"), row.names = FALSE)
write.csv(scaling, file.path(result_dir, "scaling_parameters.csv"), row.names = FALSE)
saveRDS(list(data_csv = data_csv, tree_nwk = tree_nwk, aliases_csv = aliases_csv,
             stan_md5 = unname(tools::md5sum(stan_file)),
             phylo_model = phylo_model, cases_col = cases_col, trials_col = trials_col,
             predictors = predictors, transforms = transforms, min_trials = min_trials,
             chains = chains, warmup = warmup, samples = samples, seed = seed),
        file.path(result_dir, "run_settings.rds"))

pars <- c("Intercept_mu_original", "b_mu", "Intercept_zi_original", "b_zi",
          "rho", "phi", "sd_phylo", if (phylo_model == "lambda") "lambda",
          if (phylo_model == "OU") "alpha")
s <- as.data.frame(rstan::summary(fit, pars = pars)$summary, optional = TRUE)
s$parameter <- rownames(s)
for (j in seq_along(predictors))
  s$parameter[s$parameter == paste0("b_mu[", j, "]")] <- paste0("b_", predictors[j])
write.csv(s[, c("parameter", setdiff(names(s), "parameter"))],
          file.path(result_dir, "posterior_summary.csv"), row.names = FALSE)

sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
diag_pars <- c("Intercept_mu", "b_mu", "Intercept_zi", "b_zi",
               "rho", "sd_phylo", "V_species", "w_gauss", "z_phylo",
               if (phylo_model == "lambda") "lambda",
               if (phylo_model == "OU") "alpha")
diag_summary <- rstan::summary(fit, pars = diag_pars)$summary
bfmi <- vapply(sampler, function(p) mean(diff(p[, "energy__"])^2) /
                 var(p[, "energy__"]), numeric(1))
checks <- data.frame(phylo_model = phylo_model, n_species = nrow(df),
  divergences = sum(vapply(sampler, function(p) sum(p[, "divergent__"]), numeric(1))),
  treedepth_hits = sum(vapply(sampler, function(p) sum(p[, "treedepth__"] >= 14), numeric(1))),
  minimum_bfmi = min(bfmi),
  maximum_rhat = max(diag_summary[, "Rhat"], na.rm = TRUE),
  minimum_effective_samples = min(diag_summary[, "n_eff"], na.rm = TRUE))
write.csv(checks, file.path(result_dir, "sampler_diagnostics.csv"), row.names = FALSE)
if (checks$divergences > 0 || checks$treedepth_hits > 0 ||
    checks$minimum_bfmi < 0.3 || checks$maximum_rhat > 1.01 ||
    checks$minimum_effective_samples < 100 * chains)
  warning("Sampling diagnostics need attention; see sampler_diagnostics.csv.")
# One compact Excel workbook: model setup, effects, zero inflation and diagnostics.
e <- s[match(paste0("b_", predictors), s$parameter), , drop = FALSE]
effects <- data.frame(
  Predictor = predictors, Transform = transforms,
  Median_log_odds = e[["50%"]], Lower_95 = e[["2.5%"]], Upper_95 = e[["97.5%"]],
  Median_odds_ratio = exp(e[["50%"]]), OR_lower_95 = exp(e[["2.5%"]]),
  OR_upper_95 = exp(e[["97.5%"]]), Rhat = e$Rhat, ESS = e$n_eff
)
zi_coef <- s[match(c("Intercept_zi_original", "b_zi[1]"), s$parameter), , drop = FALSE]
zi_prob <- plogis(rstan::extract(fit, pars = "Intercept_zi_original")$Intercept_zi_original)
zi_q <- unname(quantile(zi_prob, c(0.025, 0.5, 0.975)))
zero_inflation <- data.frame(
  Quantity = c("Intercept (log-odds)", "Standardized log-necropsy slope",
               "Excess-zero probability at reference necropsy effort"),
  Median = c(zi_coef[["50%"]], zi_q[2]),
  Lower_95 = c(zi_coef[["2.5%"]], zi_q[1]),
  Upper_95 = c(zi_coef[["97.5%"]], zi_q[3]),
  Interpretation = c("Log-odds of a structural zero at standardized log necropsies = 0",
    "Change in structural-zero log-odds per 1 SD higher log necropsies",
    "Structural-zero probability at standardized log necropsies = 0; not the total chance of observing zero cases")
)
model_info <- data.frame(
  Item = c("Model", "Evolutionary model", "Response (cases)", "Denominator", "Predictors",
    "Number of predictors", "Matched species", "Analyzed species", "Overdispersion prior",
    "Zero-inflation intercept prior", "Zero-inflation covariate"),
  Value = c(regression_type, phylo_model, cases_col, trials_col,
    paste(predictors, collapse = ", "), length(predictors), length(common), nrow(df),
    "rho ~ Beta(1,9)", "Normal(-2,1.5)", "Standardized log(necropsies)")
)
writexl::write_xlsx(
  list("Model summary" = model_info, "Regression effects" = effects,
       "Zero inflation" = zero_inflation, "All parameters" = s,
       "Diagnostics" = checks),
  path = file.path(result_dir, "model_summary.xlsx")
)
message("*** Finished *** :) \n Outputs: ", normalizePath(result_dir))

