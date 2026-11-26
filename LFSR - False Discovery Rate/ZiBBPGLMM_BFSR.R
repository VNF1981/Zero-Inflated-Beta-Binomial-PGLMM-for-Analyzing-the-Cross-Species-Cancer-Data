# VNF1981
# ZiBBPGLMM Bayesian false sign rate: minimal-input version
# Run from the directory containing results/:
# Rscript ZiBBPGLMM_BFSR.R

# USER SETTINGS
results_dir <- "results"
output_dir <- "BFSR_results"
manifest_file <- "model_manifest.csv"
phylo_model <- "BM"
responses <- NULL
predictors <- NULL
expected_predictors <- NULL # optional named list of intended predictors by response
family_scope <- "response" # response or pooled
bfsr_threshold <- 0.05
require_ci_excludes_zero <- TRUE
strict_diagnostics <- TRUE

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


# Only fit RDS files and labels are needed. No summary or report files are required.
read_auto <- function(root) {
  paths <- list.files(root, "^run_settings\\.rds$", recursive=TRUE, full.names=TRUE)
  fits <- list.files(root, "^fit\\.rds$", recursive=TRUE, full.names=TRUE)
  check(length(paths)>0L, "No run_settings.rds files found. Run with --make-manifest, fill the CSV, then rerun.")
  check(!length(setdiff(dirname(fits),dirname(paths))), "Some fit.rds files lack metadata. Use --make-manifest to label all intended models.")
  rows <- lapply(paths,function(p) {
    s <- readRDS(p)
    check(is.list(s) && text_vector(s$predictors) && scalar_text(s$cases_col) && scalar_text(s$phylo_model), "Invalid metadata: ",p)
    if(length(s$predictors)!=1L) return(NULL)
    data.frame(File=file.path(dirname(p),"fit.rds"), Predictor=s$predictors,
      Response=s$cases_col, Phylogenetic_model=s$phylo_model, Coefficient="b_mu[1]",
      stringsAsFactors=FALSE)
  })
  do.call(rbind,rows)
}

make_manifest <- function() {
  check(dir.exists(results_dir),"Results directory not found: ",results_dir)
  check(!file.exists(manifest_file),"Manifest already exists; edit it rather than overwriting it.")
  root <- normalizePath(results_dir,winslash="/",mustWork=TRUE)
  paths <- list.files(root,"\\.rds$",recursive=TRUE,full.names=TRUE,ignore.case=TRUE)
  paths <- paths[tolower(basename(paths))!="run_settings.rds"]
  check(length(paths)>0L,"No candidate RDS files found.")
  d <- data.frame(File=substring(paths,nchar(root)+2L),Predictor="",Response="",Phylogenetic_model="",Coefficient="",stringsAsFactors=FALSE)
  # Prepopulate labels when current-version metadata exists.
  for(i in seq_len(nrow(d))) {
    p <- file.path(dirname(paths[i]),"run_settings.rds")
    if(basename(paths[i])=="fit.rds" && file.exists(p)) {
      s <- readRDS(p)
      if(is.list(s) && length(s$predictors)==1L && scalar_text(s$cases_col) && scalar_text(s$phylo_model)) {
        d$Predictor[i] <- s$predictors; d$Response[i] <- s$cases_col
        d$Phylogenetic_model[i] <- s$phylo_model; d$Coefficient[i] <- "b_mu[1]"
      }
    }
  }
  write.csv(d,manifest_file,row.names=FALSE,na="")
  message("Created ",manifest_file,". Remove non-model rows and fill missing labels. Leave Coefficient blank for automatic matching.")
}

read_manifest <- function(root) {
  d <- read.csv(manifest_file,colClasses="character",check.names=FALSE,na.strings=NULL)
  required <- c("File","Predictor","Response","Phylogenetic_model")
  check(all(required %in% names(d)) && nrow(d)>0L,"Manifest must include File, Predictor, Response, Phylogenetic_model.")
  if(!"Coefficient" %in% names(d)) d$Coefficient <- ""
  for(n in c(required,"Coefficient")) d[[n]] <- trimws(d[[n]])
  for(n in required) check(!anyNA(d[[n]]) && all(nzchar(d[[n]])),"Fill every ",n," cell in the manifest.")
  d$Coefficient[is.na(d$Coefficient)] <- ""
  # Paths in the manifest are relative to results_dir, not the script directory.
  d$File <- file.path(root,d$File)
  d
}

read_model <- function(row) {
  path <- row$File
  obj <- readRDS(path)
  if(inherits(obj,"brmsfit")) {
    check(requireNamespace("brms",quietly=TRUE),"Install brms to read older brms models.")
    fit <- obj$fit
    kind <- "brmsfit"
  } else { fit <- obj; kind <- "stanfit" }
  check(inherits(fit,"stanfit"),"Unsupported fitted object in ",path," (requires RStan stanfit or brmsfit with RStan backend).")
  pn <- fit@sim$fnames_oi
  check(length(pn)>0L,"No posterior draws: ",path)
  coeff <- row$Coefficient
  if(!nzchar(coeff)) {
    if(kind=="stanfit") {
      candidates <- intersect(c("b_mu[1]","b_mu"),pn)
      # Automatic selection is allowed only for a scalar or one-element b_mu.
      if(sum(grepl("^b_mu\\[",pn))>1L) candidates <- character()
    } else {
      candidates <- intersect(c(paste0("b_",row$Predictor),paste0("b_mu_",row$Predictor)),pn)
    }
    check(length(candidates)==1L,"Cannot uniquely identify coefficient for ",row$Predictor," in ",path,
      ". Set Coefficient explicitly in the manifest. Available b parameters: ",paste(grep("^b_",pn,value=TRUE),collapse=", "))
    coeff <- candidates[1L]
  }
  check(coeff %in% pn,"Coefficient ",coeff," not found in ",path,". Available b parameters: ",paste(grep("^b_",pn,value=TRUE),collapse=", "))
  if(kind=="stanfit" && grepl("^b_mu\\[",coeff)) check(sum(grepl("^b_mu\\[",pn))==1L,"Multiple b_mu coefficients: use single-predictor models only: ",path)
  a <- rstan::extract(fit,pars=coeff,permuted=FALSE,inc_warmup=FALSE)
  check(length(dim(a))==3L && dim(a)[3L]==1L,"Unexpected coefficient draw shape: ",path)
  beta <- as.numeric(a); stats <- posterior_stats(beta)
  sm <- rstan::summary(fit,pars=coeff)$summary
  check(nrow(sm)==1L && all(c("Rhat","n_eff") %in% colnames(sm)),"Unable to extract coefficient diagnostics: ",path)
  rh <- unname(sm[1,"Rhat"]); ess <- unname(sm[1,"n_eff"]); chains <- dim(a)[2L]
  sampler <- tryCatch(rstan::get_sampler_params(fit,inc_warmup=FALSE),error=function(e) NULL)
  div <- NA_real_; bfmi <- NA_real_; depth_hits <- NA_real_
  if(length(sampler)>0L && all(vapply(sampler,function(x) all(c("divergent__","energy__","treedepth__") %in% colnames(x)),logical(1)))) {
    div <- sum(vapply(sampler,function(x) sum(x[,"divergent__"]),numeric(1)))
    bfmi <- min(vapply(sampler,function(x) mean(diff(x[,"energy__"])^2)/var(x[,"energy__"]),numeric(1)))
    # RStan NUTS defaults to 10 when no explicit maximum is supplied.
    limits <- vapply(fit@stan_args,function(x) {
      z <- x$control$max_treedepth
      if(is.null(z)) 10 else as.numeric(z)
    },numeric(1))
    if(length(limits)==length(sampler)) depth_hits <- sum(vapply(seq_along(sampler),function(i) sum(sampler[[i]][,"treedepth__"]>=limits[i]),numeric(1)))
  }
  problems <- c(if(chains<2) "fewer than two chains",
    if(!is.finite(rh) || rh>1.01) "Rhat missing or > 1.01",
    if(!is.finite(ess) || ess<100*chains) "ESS missing or < 100 per chain",
    if(!is.finite(div) || div>0) "divergences present or unavailable",
    if(!is.finite(bfmi) || bfmi<.3) "BFMI missing or < 0.3",
    if(!is.finite(depth_hits) || depth_hits>0) "treedepth hits present or unavailable")
  data.frame(Predictor=row$Predictor,Response=row$Response,Phylogenetic_model=row$Phylogenetic_model,
    Coefficient=coeff,Model_type=kind,stats,Chains=chains,Rhat=rh,ESS=ess,Divergences=div,
    Minimum_BFMI=bfmi,Treedepth_hits=depth_hits,QC_pass=!length(problems),QC_issues=paste(problems,collapse="; "),File=path)
}

main <- function() {
  check(dir.exists(results_dir),"Results directory not found: ",results_dir)
  check(scalar_text(phylo_model) && phylo_model %in% c("BM","lambda","OU"),"Select BM, lambda or OU.")
  check(family_scope %in% c("response","pooled"),"Invalid family_scope.")
  check(scalar_number(bfsr_threshold) && bfsr_threshold>0 && bfsr_threshold<=.5,"Invalid BFSR threshold.")
  for(x in list(require_ci_excludes_zero,strict_diagnostics)) check(is.logical(x)&&length(x)==1L&&!is.na(x),"Invalid logical setting.")
  check(requireNamespace("rstan",quietly=TRUE),"Install rstan: install.packages('rstan')")
  root <- normalizePath(results_dir,winslash="/",mustWork=TRUE)
  manual <- file.exists(manifest_file)
  d <- if(manual) read_manifest(root) else read_auto(root)
  check(!is.null(d) && nrow(d)>0L,"No models found.")
  d <- d[d$Phylogenetic_model==phylo_model & (if(is.null(responses)) TRUE else d$Response %in% responses) & (if(is.null(predictors)) TRUE else d$Predictor %in% predictors),,drop=FALSE]
  check(nrow(d)>0L,"No models match your settings.")
  check(all(file.exists(d$File)),"Missing fitted files: ",paste(d$File[!file.exists(d$File)],collapse="; "))
  d$File <- normalizePath(d$File,winslash="/",mustWork=TRUE)
  check(!anyDuplicated(d$File),"Repeated fitted model file in selected family.")
  check(!anyDuplicated(d[c("Response","Predictor")]),"Duplicate predictor/response models. Select one prespecified model per association.")
  if(!is.null(responses)) check(setequal(responses,unique(d$Response)),"A requested response is missing.")
  if(!is.null(expected_predictors)) {
    check(setequal(names(expected_predictors),unique(d$Response)),"Expected responses do not match discovered responses.")
    for(r in names(expected_predictors)) check(setequal(expected_predictors[[r]],d$Predictor[d$Response==r]),"Incomplete or unexpected predictor family: ",r)
  } else if(!is.null(predictors)) {
    for(r in unique(d$Response)) check(setequal(predictors,d$Predictor[d$Response==r]),"Incomplete predictor family: ",r)
  } else warning("Verify that all intended models are represented. Missing unlisted analyses cannot be detected automatically.",call.=FALSE)
  check(!file.exists(output_dir),"Output directory already exists; choose a new output_dir.")
  res <- do.call(rbind,lapply(seq_len(nrow(d)),function(i) {
    message("Reading ",i,"/",nrow(d),": ",d$Predictor[i]); read_model(d[i,,drop=FALSE])
  }))
  res$Family <- if(family_scope=="pooled") "All_responses" else res$Response
  groups <- lapply(split(res,res$Family),rank_family,alpha=bfsr_threshold,require_ci=require_ci_excludes_zero)
  res <- do.call(rbind,unname(groups)); rownames(res) <- NULL
  res$Family_QC_pass <- ave(as.integer(res$QC_pass),res$Family,FUN=min)==1L
  res$Reportable_discovery <- res$Final_discovery & res$Family_QC_pass
  discoveries <- res[res$Reportable_discovery,,drop=FALSE]
  overview <- do.call(rbind,lapply(groups,function(x) {
    selected <- x$Final_discovery; raw <- x$BFSR_pass
    data.frame(Family=x$Family[1],N_predictors=length(unique(x$Predictor)),N_associations=nrow(x),N_BFSR_pass=sum(raw),N_final=sum(selected),
      N_reportable=if(all(x$QC_pass)) sum(selected) else 0L,Achieved_BFSR=if(any(raw)) mean(x$lfsr[raw]) else NA_real_,
      Achieved_final_BFSR=if(any(selected)) mean(x$lfsr[selected]) else NA_real_,All_QC_pass=all(x$QC_pass))
  }))
  qc <- res[c("Response","Predictor","Coefficient","Chains","Rhat","ESS","Divergences","Minimum_BFMI","Treedepth_hits","QC_pass","QC_issues","File")]
  settings <- data.frame(Setting=c("Method","Threshold","Family_scope","Phylogenetic_model","Require_CI","Metadata_source","Strict_diagnostics"),
    Value=c("Bayesian false sign rate (not point-null FDR)",bfsr_threshold,family_scope,phylo_model,require_ci_excludes_zero,if(manual) "manifest" else "run_settings.rds",strict_diagnostics))
  check(dir.create(output_dir,recursive=TRUE),"Cannot create output directory.")
  tables <- list(Results=res,Discoveries=discoveries,Summary=overview,Diagnostics=qc,Inventory=d,Settings=settings)
  for(n in names(tables)) write.csv(tables[[n]],file.path(output_dir,paste0("BFSR_",tolower(n),".csv")),row.names=FALSE)
  if(requireNamespace("writexl",quietly=TRUE)) writexl::write_xlsx(tables,file.path(output_dir,"BFSR_results.xlsx"))
  else message("Optional writexl package unavailable; CSV reports were generated.")
  grDevices::pdf(file.path(output_dir,"BFSR_cumulative_plot.pdf"),width=8,height=5)
  tryCatch(for(x in groups) {
    plot(x$Rank,x$Cumulative_BFSR,type="b",pch=19,xlab="Predictor rank",ylab="Cumulative Bayesian false sign rate",main=x$Family[1],ylim=c(0,max(bfsr_threshold,x$Cumulative_BFSR)*1.08))
    abline(h=bfsr_threshold,lty=2); points(x$Rank[x$Final_discovery],x$Cumulative_BFSR[x$Final_discovery],pch=21,bg="grey")
  },finally=grDevices::dev.off())
  print(overview,row.names=FALSE)
  if(any(!res$QC_pass)) {
    if(strict_diagnostics) fail("Sampling QC failed. Reports saved, but affected families have no reportable discoveries. See BFSR_diagnostics.csv.")
    warning("Sampling QC failed; affected families have no reportable discoveries.",call.=FALSE)
  }
  message("Finished: ",normalizePath(output_dir))
}
args <- commandArgs(trailingOnly=TRUE)
if("--self-test" %in% args) self_test() else if("--make-manifest" %in% args) make_manifest() else main()
