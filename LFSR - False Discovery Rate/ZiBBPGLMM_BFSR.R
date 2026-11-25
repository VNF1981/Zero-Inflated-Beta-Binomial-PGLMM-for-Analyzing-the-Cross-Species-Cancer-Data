# VNF1981 - ZiBBPGLMM V2 downstream Bayesian false sign rate analysis
# Run: Rscript ZiBBPGLMM_BFSR.R
# Test calculation logic: Rscript ZiBBPGLMM_BFSR.R --self-test

# -------------------- USER SETTINGS --------------------
results_dir <- "results"
output_dir <- "BFSR_results"
phylo_model <- "BM"                 # BM, lambda, OU; select before analysis
responses <- NULL                   # NULL = all discovered responses
predictors <- NULL                  # NULL = all discovered single-predictor models
# Optional named list of intended predictors for each response, e.g.:
# expected_predictors <- list(NeoplasiaWithDenominators = c("adult_weight.g.", "residLongevity"))
expected_predictors <- NULL
family_scope <- "response"          # response or pooled
bfsr_threshold <- 0.05
require_ci_excludes_zero <- TRUE
strict_diagnostics <- TRUE          # FALSE permits exploratory output, flagged in reports
# -------------------------------------------------------

fail <- function(...) stop(..., call. = FALSE)
check <- function(ok, ...) if (!isTRUE(ok)) fail(...)
text_vector <- function(x) is.character(x) && length(x) > 0L && !anyNA(x) && all(nzchar(trimws(x))) && !anyDuplicated(x)
scalar_text <- function(x) text_vector(x) && length(x) == 1L
scalar_number <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
read_csv <- function(path) read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)

posterior_stats <- function(beta) {
  check(is.numeric(beta) && length(beta) >= 2L && all(is.finite(beta)), "Invalid posterior draws.")
  ppos <- mean(beta > 0); pneg <- mean(beta < 0); pzero <- mean(beta == 0)
  q <- unname(quantile(beta, c(.025, .5, .975), type = 7))
  data.frame(N_draws = length(beta), Beta_mean = mean(beta), Beta_median = q[2],
             Lower_95 = q[1], Upper_95 = q[3], P_beta_gt_0 = ppos,
             P_beta_lt_0 = pneg, P_beta_eq_0 = pzero,
             Direction = if (ppos >= pneg) "positive" else "negative",
             lfsr = 1 - max(ppos, pneg),
             CI_pass = q[1] > 0 || q[3] < 0)
}

rank_family <- function(d, alpha, require_ci) {
  check(nrow(d) > 0L && all(is.finite(d$lfsr)) && all(d$lfsr >= 0 & d$lfsr <= 1), "Invalid family error probabilities.")
  check(!anyNA(d$CI_pass), "Missing credible interval decision.")
  d <- d[order(d$lfsr, d$Response, d$Predictor, method = "radix"), , drop = FALSE]
  rownames(d) <- NULL
  d$Rank <- seq_len(nrow(d))
  d$Cumulative_BFSR <- cumsum(d$lfsr) / d$Rank
  # Keep equal lfsr values together: evaluate the cumulative mean at the end of each tie.
  tie_end <- rep(cumsum(rle(d$lfsr)$lengths), rle(d$lfsr)$lengths)
  d$Tie_group_BFSR <- d$Cumulative_BFSR[tie_end]
  d$BFSR_pass <- d$Tie_group_BFSR <= alpha + 1e-12
  d$Final_discovery <- d$BFSR_pass & (!require_ci | d$CI_pass)
  # The optional CI screen changes the selected set; verify its actual posterior risk.
  # If removing a member raises the average above alpha, drop entire worst-lfsr groups.
  while (any(d$Final_discovery) && mean(d$lfsr[d$Final_discovery]) > alpha + 1e-12) {
    worst <- max(d$lfsr[d$Final_discovery])
    d$Final_discovery[d$lfsr == worst] <- FALSE
  }
  d
}

self_test <- function() {
  d <- data.frame(Response = "N", Predictor = letters[1:7],
                  lfsr = c(.0025,.0063,.057,.067,.072,.077,.106),
                  CI_pass = c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE))
  a <- rank_family(d, .05, FALSE); b <- rank_family(d, .05, TRUE)
  stopifnot(sum(a$BFSR_pass) == 6L, sum(a$Final_discovery) == 6L,
            sum(b$Final_discovery) == 2L, abs(a$Cumulative_BFSR[6] - .0469666666666667) < 1e-12)
  z <- posterior_stats(c(rep(-1, 10), rep(1, 990)))
  stopifnot(abs(z$lfsr - .01) < 1e-12, z$CI_pass, z$Direction == "positive")
  z <- posterior_stats(c(rep(-1, 990), rep(1, 10)))
  stopifnot(abs(z$lfsr - .01) < 1e-12, z$CI_pass, z$Direction == "negative")
  z <- posterior_stats(c(-1, 0, 1)); stopifnot(abs(z$lfsr - 2/3) < 1e-12)
  d <- data.frame(Response = "N", Predictor = letters[1:3], lfsr = c(.01,.08,.08), CI_pass = TRUE)
  stopifnot(sum(rank_family(d, .05, FALSE)$Final_discovery) == 1L)
  d$lfsr <- c(.1,.2,.3); stopifnot(!any(rank_family(d, .05, FALSE)$Final_discovery))
  d$lfsr <- c(.001,.06,.07); d$CI_pass <- c(FALSE,TRUE,TRUE)
  stopifnot(!any(rank_family(d, .05, TRUE)$Final_discovery))
  d <- d[1,,drop=FALSE]; d$lfsr <- .05; d$CI_pass <- TRUE
  stopifnot(rank_family(d,.05,TRUE)$Final_discovery)
  message("All built-in calculation self-tests passed.")
}

read_metadata <- function(path) {
  s <- readRDS(path)
  needed <- c("data_csv","tree_nwk","aliases_csv","stan_md5","phylo_model","cases_col","trials_col","predictors","transforms","min_trials","chains","warmup","samples","seed")
  check(is.list(s) && all(needed %in% names(s)), "Incompatible metadata: ", path)
  for (key in c("data_csv","tree_nwk","stan_md5","phylo_model","cases_col","trials_col"))
    check(scalar_text(s[[key]]), "Invalid metadata field ", key, ": ", path)
  check(text_vector(s$predictors) && is.character(s$transforms) && !anyNA(s$transforms) &&
          length(s$predictors) == length(s$transforms) && all(s$transforms %in% c("identity","log10")), "Invalid predictors/transforms: ", path)
  check(s$phylo_model %in% c("BM","lambda","OU"), "Unknown evolutionary model: ", path)
  for (key in c("min_trials","chains","warmup","samples","seed"))
    check(scalar_number(s[[key]]) && s[[key]] >= 0 && s[[key]] == floor(s[[key]]), "Invalid metadata field ", key, ": ", path)
  check(s$chains >= 1 && s$samples >= 2, "Invalid sampling settings: ", path)
  check(is.null(s$aliases_csv) || scalar_text(s$aliases_csv), "Invalid alias metadata: ", path)
  s
}

read_model <- function(path, s) {
  folder <- dirname(path)
  files <- file.path(folder, c("fit.rds","analysis_species.csv","scaling_parameters.csv","posterior_summary.csv","sampler_diagnostics.csv","model_summary.xlsx"))
  check(all(file.exists(files)), "Incomplete model: ", folder, "\nMissing: ", paste(basename(files[!file.exists(files)]), collapse = ", "))
  check(all(file.info(files)$size > 0), "Empty output file: ", folder)
  species <- read_csv(files[2]); scaling <- read_csv(files[3]); post <- read_csv(files[4]); diag <- read_csv(files[5])
  check("Species" %in% names(species) && text_vector(species$Species) && nrow(species) >= 20L, "Invalid species file: ", folder)
  check(all(c("predictor","transform","center","sd") %in% names(scaling)) && nrow(scaling) == 2L &&
          identical(as.character(scaling$predictor), c(s$predictors,"log_trials")) &&
          identical(as.character(scaling$transform), c(s$transforms,"ln")) &&
          all(is.finite(scaling$center)) && all(is.finite(scaling$sd) & scaling$sd > 0), "Inconsistent scaling file: ", folder)
  fields <- c("n_species","divergences","treedepth_hits","minimum_bfmi","maximum_rhat","minimum_effective_samples")
  check(nrow(diag) == 1L && all(c("phylo_model",fields) %in% names(diag)), "Invalid diagnostics structure: ", folder)
  check(all(vapply(diag[fields], scalar_number, logical(1))) && diag$phylo_model == s$phylo_model && diag$n_species == nrow(species), "Inconsistent diagnostics: ", folder)
  check(all(c("parameter","mean","2.5%","50%","97.5%","Rhat","n_eff") %in% names(post)), "Invalid posterior summary: ", folder)
  row <- post[post$parameter == paste0("b_",s$predictors), , drop = FALSE]
  check(nrow(row) == 1L && all(vapply(row[c("mean","2.5%","50%","97.5%","Rhat","n_eff")], scalar_number, logical(1))), "Invalid coefficient summary: ", folder)
  fit <- readRDS(files[1])
  check(inherits(fit,"stanfit"), "fit.rds is not an RStan stanfit: ", folder)
  draws <- rstan::extract(fit, pars = "b_mu", permuted = TRUE, inc_warmup = FALSE)$b_mu
  check(is.matrix(draws) && ncol(draws) == 1L && nrow(draws) == s$chains * s$samples, "Unexpected b_mu posterior shape or draw count: ", folder)
  beta <- draws[,1L]; stats <- posterior_stats(beta)
  check(isTRUE(all.equal(unname(c(stats$Beta_mean,stats$Lower_95,stats$Beta_median,stats$Upper_95)),
                          unname(as.numeric(unlist(row[1,c("mean","2.5%","50%","97.5%"), drop = FALSE], use.names = FALSE))), tolerance = 1e-6)),
        "Posterior draws and saved summary disagree: ", folder)
  # Recompute coefficient diagnostics directly from the fitted object.
  sm <- rstan::summary(fit, pars = "b_mu")$summary
  check(nrow(sm) == 1L && all(c("Rhat","n_eff") %in% colnames(sm)) && all(is.finite(sm[1,c("Rhat","n_eff")])), "Invalid fitted coefficient diagnostics: ", folder)
  check(isTRUE(all.equal(as.numeric(sm[1,c("Rhat","n_eff")]), c(row$Rhat,row$n_eff), tolerance = 1e-6)), "Fitted and saved coefficient diagnostics disagree: ", folder)
  problems <- c(if (s$chains < 2) "fewer than two chains",
                if (diag$divergences != 0) "divergences",
                if (diag$treedepth_hits != 0) "treedepth hits",
                if (diag$minimum_bfmi < .3) "BFMI < 0.3",
                if (max(diag$maximum_rhat,row$Rhat) > 1.01) "Rhat > 1.01",
                if (min(diag$minimum_effective_samples,row$n_eff) < 100*s$chains) "ESS < 100 per chain")
  data.frame(Predictor = s$predictors, Response = s$cases_col, Denominator = s$trials_col,
             Phylogenetic_model = s$phylo_model, Transform = s$transforms, N_species = nrow(species),
             stats, Rhat = row$Rhat, ESS = row$n_eff, Chains = s$chains,
             Divergences = diag$divergences, Treedepth_hits = diag$treedepth_hits,
             Minimum_BFMI = diag$minimum_bfmi, Maximum_Rhat = diag$maximum_rhat,
             Minimum_ESS = diag$minimum_effective_samples, QC_pass = !length(problems),
             QC_issues = paste(problems,collapse = "; "),
             Data_file = s$data_csv, Tree_file = s$tree_nwk,
             Alias_file = if (is.null(s$aliases_csv)) "<none>" else s$aliases_csv,
             Stan_MD5 = s$stan_md5, Min_trials = s$min_trials,
             Model_directory = folder, stringsAsFactors = FALSE)
}

main <- function() {
  check(dir.exists(results_dir), "Results directory does not exist: ", results_dir)
  check(scalar_text(output_dir), "Invalid output directory.")
  check(scalar_text(phylo_model) && phylo_model %in% c("BM","lambda","OU"), "Select BM, lambda, or OU.")
  check(scalar_text(family_scope) && family_scope %in% c("response","pooled"), "Invalid family_scope.")
  check(scalar_number(bfsr_threshold) && bfsr_threshold > 0 && bfsr_threshold <= .5, "Invalid BFSR threshold.")
  for (x in list(require_ci_excludes_zero, strict_diagnostics)) check(is.logical(x) && length(x) == 1L && !is.na(x), "Invalid logical setting.")
  check(is.null(responses) || text_vector(responses), "Invalid responses setting.")
  check(is.null(predictors) || text_vector(predictors), "Invalid predictors setting.")
  check(is.null(expected_predictors) || (is.list(expected_predictors) && text_vector(names(expected_predictors)) &&
          all(vapply(expected_predictors,text_vector,logical(1)))), "Invalid expected_predictors list.")
  for (pkg in c("rstan","writexl")) check(requireNamespace(pkg,quietly=TRUE), "Install package: ", pkg)
  check(packageVersion("rstan") >= "2.26.0", "rstan >= 2.26 is required.")
  root <- normalizePath(results_dir,mustWork=TRUE)
  paths <- list.files(root, pattern = "^run_settings\\.rds$", recursive = TRUE, full.names = TRUE)
  fit_paths <- list.files(root, pattern = "^fit\\.rds$", recursive = TRUE, full.names = TRUE)
  orphan <- setdiff(dirname(fit_paths),dirname(paths))
  check(!length(orphan), "Fits without run_settings.rds found: ", paste(orphan,collapse="; "))
  check(length(paths) > 0L, "No run_settings.rds found under ", root)
  metadata <- lapply(paths,read_metadata)
  eligible <- vapply(metadata,function(s) length(s$predictors)==1L && s$phylo_model==phylo_model &&
                      (is.null(responses) || s$cases_col %in% responses) &&
                      (is.null(predictors) || s$predictors %in% predictors) &&
                      (is.null(expected_predictors) || (s$cases_col %in% names(expected_predictors) &&
                        s$predictors %in% expected_predictors[[s$cases_col]])), logical(1))
  inventory <- do.call(rbind,lapply(seq_along(paths),function(i) {
    s <- metadata[[i]]
    data.frame(Response=s$cases_col, Phylogenetic_model=s$phylo_model, Predictors=paste(s$predictors,collapse="; "),
               Transforms=paste(s$transforms,collapse="; "), Eligible=eligible[i], Model_directory=dirname(paths[i]))
  }))
  paths <- paths[eligible]; metadata <- metadata[eligible]
  check(length(paths)>0L,"No eligible single-predictor models.")
  keys <- do.call(rbind,lapply(metadata,function(s) data.frame(Response=s$cases_col,Predictor=s$predictors)))
  check(!anyDuplicated(keys),"Duplicate response/predictor models found (including alternative transformations). Keep one prespecified model per predictor.")
  actual_responses <- sort(unique(keys$Response))
  if (!is.null(responses)) check(setequal(responses,actual_responses),"Requested response missing.")
  if (!is.null(expected_predictors)) check(setequal(names(expected_predictors),actual_responses),"Expected response missing.")
  for (response in actual_responses) {
    expected <- if (!is.null(expected_predictors)) expected_predictors[[response]] else predictors
    actual <- keys$Predictor[keys$Response==response]
    if (!is.null(expected)) check(setequal(expected,actual),"Incomplete family for ",response,". Missing: ",paste(setdiff(expected,actual),collapse=", "),"; unexpected: ",paste(setdiff(actual,expected),collapse=", "))
  }
  if (is.null(expected_predictors) && is.null(predictors)) warning("Automatic discovery cannot detect intended models that never produced files. Verify the family inventory.",call.=FALSE)
  d <- do.call(rbind,lapply(seq_along(paths),function(i) {
    message("Reading ", i,"/",length(paths),": ",dirname(paths[i])); read_model(paths[i],metadata[[i]])
  }))
  # Do not combine incompatible datasets or model specifications within a response.
  for (response in actual_responses) {
    sub <- d[d$Response==response,,drop=FALSE]
    for (field in c("Denominator","Data_file","Tree_file","Alias_file","Stan_MD5","Min_trials"))
      check(length(unique(sub[[field]]))==1L,"Mixed ",field," in response ",response,"; use a homogeneous results directory.")
  }
  check(!file.exists(output_dir),"Output directory already exists. Choose a new path to avoid overwriting: ",output_dir)
  check(dir.create(output_dir,recursive=TRUE),"Cannot create output directory.")
  write.csv(inventory,file.path(output_dir,"BFSR_inventory.csv"),row.names=FALSE)
  qc_cols <- c("Response","Predictor","Phylogenetic_model","N_species","N_draws","Chains","Rhat","ESS","Divergences","Treedepth_hits","Minimum_BFMI","Maximum_Rhat","Minimum_ESS","QC_pass","QC_issues","Model_directory")
  write.csv(d[qc_cols],file.path(output_dir,"BFSR_diagnostics.csv"),row.names=FALSE)
  if (strict_diagnostics && any(!d$QC_pass)) fail("Sampling QC failed. See BFSR_diagnostics.csv. No discovery report generated. Resolve diagnostics and use a new output directory.")
  if (any(!d$QC_pass)) warning("QC failures retained; results are exploratory only.",call.=FALSE)
  d$Family <- if (family_scope=="pooled") rep("All_responses",nrow(d)) else d$Response
  groups <- split(d,d$Family)
  ranked <- lapply(groups,rank_family,alpha=bfsr_threshold,require_ci=require_ci_excludes_zero)
  res <- do.call(rbind,unname(ranked)); rownames(res) <- NULL
  res$Family_QC_pass <- ave(as.integer(res$QC_pass),res$Family,FUN=min)==1L
  res$Reportable_discovery <- res$Final_discovery & res$Family_QC_pass
  summary <- do.call(rbind,lapply(ranked,function(x) {
    selected <- x$Final_discovery; raw <- x$BFSR_pass
    data.frame(Family=x$Family[1],N_predictors=length(unique(x$Predictor)),N_associations=nrow(x),
               N_BFSR_pass=sum(raw),N_CI_pass=sum(x$CI_pass),N_final=sum(selected),
               N_reportable=if (all(x$QC_pass)) sum(selected) else 0L,
               Achieved_BFSR=if(any(raw)) mean(x$lfsr[raw]) else NA_real_,
               Achieved_final_BFSR=if(any(selected)) mean(x$lfsr[selected]) else NA_real_,
               All_QC_pass=all(x$QC_pass))
  })); rownames(summary) <- NULL
  discoveries <- res[res$Reportable_discovery,,drop=FALSE]
  settings <- data.frame(Setting=c("Method","Threshold","Family_scope","Phylogenetic_model","Require_CI","Strict_diagnostics","Expected_family_verified","Results_directory","R_version","RStan_version","Interpretation"),
                         Value=c("Bayesian false sign rate (not point-null FDR)",bfsr_threshold,family_scope,phylo_model,require_ci_excludes_zero,strict_diagnostics,
                                 !is.null(expected_predictors)||!is.null(predictors),root,R.version.string,as.character(packageVersion("rstan")),"Posterior expected proportion of incorrect effect signs, conditional on fitted models and priors"))
  write.csv(res,file.path(output_dir,"BFSR_results.csv"),row.names=FALSE)
  write.csv(discoveries,file.path(output_dir,"BFSR_discoveries.csv"),row.names=FALSE)
  write.csv(summary,file.path(output_dir,"BFSR_summary.csv"),row.names=FALSE)
  write.csv(settings,file.path(output_dir,"BFSR_settings.csv"),row.names=FALSE)
  writexl::write_xlsx(list(Results=res,Discoveries=discoveries,Summary=summary,Diagnostics=d[qc_cols],Inventory=inventory,Settings=settings),file.path(output_dir,"BFSR_results.xlsx"))
  grDevices::pdf(file.path(output_dir,"BFSR_cumulative_plot.pdf"),width=8,height=5)
  tryCatch(for (x in ranked) {
    plot(x$Rank,x$Cumulative_BFSR,type="b",pch=19,xlab="Rank (ascending lfsr)",ylab="Cumulative Bayesian false sign rate",
         main=x$Family[1],ylim=c(0,max(bfsr_threshold,x$Cumulative_BFSR)*1.08))
    abline(h=bfsr_threshold,lty=2,col="red")
    points(x$Rank[x$Final_discovery],x$Cumulative_BFSR[x$Final_discovery],pch=21,bg="darkgreen",col="black",cex=1.4)
    legend("topleft",legend=c("Cumulative BFSR","Threshold","Final discovery"),pch=c(19,NA,21),lty=c(1,2,NA),col=c("black","red","black"),pt.bg=c(NA,NA,"darkgreen"),bty="n")
  },finally=grDevices::dev.off())
  print(summary,row.names=FALSE)
  message("Finished. Outputs: ",normalizePath(output_dir))
  invisible(res)
}

if ("--self-test" %in% commandArgs(trailingOnly=TRUE)) self_test() else main()
