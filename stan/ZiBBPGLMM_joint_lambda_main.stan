// ZiBB-PGLMM, Pagel's lambda -- exact reparameterisation (October 2026).
//
// SAME generative model, SAME likelihood, SAME priors and SAME data block as
// ZiBB_lambda_ORIGINAL.stan / ZiBB_lambda_EIGEN_CANDIDATE.stan.
// Only the coordinates in which HMC moves are changed:
//
// (1) Variance coordinates.  Two species-level variance terms enter the model:
//     the Gaussian latent variance sd_phylo^2 (split by lambda into a
//     phylogenetic part and an i.i.d. "nugget" part) and the beta-binomial
//     nugget, whose approximate logit-scale variance is kappa / phi, with
//     kappa = 1 / (pbar (1 - pbar)) a fixed data constant.  HMC samples
//         V_species = sd_phylo^2 + kappa / phi     (total; well identified)
//         w_gauss   = sd_phylo^2 / V_species       (Gaussian share; weakly identified)
//         lambda                                   (unchanged)
//     The priors are still stated on sd_phylo, lambda and rho, and the
//     log-Jacobian of (V_species, w_gauss, lambda) -> (sd_phylo, lambda, rho)
//     is added, so the joint prior on (sd_phylo, lambda, rho) is unchanged.
//     kappa is only a coordinate scale; any positive constant gives the same
//     posterior.
//
// (2) Latent field.  Fixed eigenbasis A_BM = Q diag(d) Q' (as in the eigen
//     candidate) with information-scaled partial non-centring:
//         e_k  = tau_k * c_k * z_phylo[k],   z_phylo[k] ~ N(0, 1 / c_k^2)
//     which is identical in distribution to e_k ~ N(0, tau_k^2), where
//         tau_k^2 = sd_phylo^2 (1 - lambda + lambda d_k),
//         c_k     = (1 + tau_k^2 * I_k)^(-1/2),
//         I_k     = sum_i Q[i,k]^2 * trials_i pbar (1 - pbar) / (1 + (trials_i - 1) rho)
//     (I_k approximates the beta-binomial information about eigen coordinate k).
//     c_k -> 1 gives the usual non-centred map; large I_k moves toward centring.
//     Purpose: keep the conditional scale of z_phylo near 1 wherever
//     (sd_phylo, rho) sit, which removes the energy (BFMI) pathology.
//
// Remaining, by design: w_gauss (how species-level overdispersion is split
// between the Gaussian nugget sd_phylo^2 (1 - lambda) and the beta-binomial
// nugget) is only weakly identified with one count per species.  Its posterior
// is partly prior-driven and may need more draws to reach ESS >= 400, but it no
// longer produces divergences or low BFMI.
//
// Output names used by the runner are preserved: Intercept_mu, b_mu,
// Intercept_zi, b_zi, rho, phi, sd_phylo, lambda, z_phylo, u_phylo, log_lik,
// y_rep, Intercept_mu_original, Intercept_zi_original.
// NOTE: z_phylo is now the partially non-centred latent vector; it is still
// the sampled latent block and is still N(0, I)-like a posteriori.

functions {
  real zibeta_binomial2_lpmf(int y, int trials, real mu, real phi, real zi) {
    real bb_lp = beta_binomial_lpmf(y | trials, mu * phi, (1 - mu) * phi);
    if (y == 0) {
      return log_sum_exp(
        bernoulli_lpmf(1 | zi),
        bernoulli_lpmf(0 | zi) + bb_lp
      );
    }
    return bernoulli_lpmf(0 | zi) + bb_lp;
  }
}

data {
  int<lower=1> N;
  array[N] int<lower=0> y;
  array[N] int<lower=1> trials;

  int<lower=1> K_mu;
  matrix[N, K_mu] X_mu;

  int<lower=1> K_zi;
  matrix[N, K_zi] X_zi;

  matrix[N, N] A_BM;

  real<lower=0> rho_prior_a;
  real<lower=0> rho_prior_b;
  real zi_intercept_prior_mean;
  real<lower=0> zi_intercept_prior_sd;
}

transformed data {
  vector[K_mu] X_mu_means;
  vector[K_zi] X_zi_means;
  matrix[N, K_mu] X_mu_c;
  matrix[N, K_zi] X_zi_c;

  matrix[N, N] Q_BM = eigenvectors_sym(A_BM);
  vector[N] eig_BM = eigenvalues_sym(A_BM);
  matrix[N, N] Q_BM_sq_t = (Q_BM .* Q_BM)';   // transposed once; used every evaluation

  // Fixed coordinate constants (data only; they do not change the model).
  real pbar = (sum(y) + 0.5) / (sum(trials) + 1.0);
  real kappa = 1.0 / (pbar * (1.0 - pbar));
  vector[N] info_num = to_vector(trials) * (pbar * (1.0 - pbar));
  vector[N] trials_m1 = to_vector(trials) - 1.0;

  for (k in 1:K_mu) {
    X_mu_means[k] = mean(col(X_mu, k));
    X_mu_c[, k] = X_mu[, k] - X_mu_means[k];
  }
  for (k in 1:K_zi) {
    X_zi_means[k] = mean(col(X_zi, k));
    X_zi_c[, k] = X_zi[, k] - X_zi_means[k];
  }

  if (min(eig_BM) < -1e-8) reject("A_BM must be positive semidefinite");
}

parameters {
  real Intercept_mu;
  vector[K_mu] b_mu;

  real Intercept_zi;
  vector[K_zi] b_zi;

  real<lower=0> V_species;               // sd_phylo^2 + kappa / phi
  real<lower=0, upper=1> w_gauss;        // sd_phylo^2 / V_species
  real<lower=0, upper=1> lambda;

  vector[N] z_phylo;
}

transformed parameters {
  real<lower=0> sd_phylo = sqrt(w_gauss * V_species);
  real<lower=0> phi = kappa / ((1 - w_gauss) * V_species);
  real<lower=0, upper=1> rho = 1 / (1 + phi);

  vector[N] c_pnc;
  vector[N] u_phylo;
  vector[N] eta_mu;
  vector[N] eta_zi;
  vector[N] mu;
  vector[N] zi;

  {
    vector[N] tau_sq = square(sd_phylo)
                       * (rep_vector(1.0 - lambda + 1e-10, N) + lambda * eig_BM);
    vector[N] info_eig = Q_BM_sq_t * (info_num ./ (1.0 + trials_m1 * rho));
    c_pnc = inv_sqrt(1.0 + tau_sq .* info_eig);
    u_phylo = Q_BM * (sqrt(tau_sq) .* c_pnc .* z_phylo);
  }

  eta_mu = Intercept_mu + X_mu_c * b_mu + u_phylo;
  eta_zi = Intercept_zi + X_zi_c * b_zi;
  mu = inv_logit(eta_mu);
  zi = inv_logit(eta_zi);
}

model {
  Intercept_mu ~ normal(0, 1.5);
  b_mu ~ normal(0, 1);

  Intercept_zi ~ normal(zi_intercept_prior_mean, zi_intercept_prior_sd);
  b_zi ~ normal(0, 1);

  // Original priors, on the original parameters (unchanged).
  rho ~ beta(rho_prior_a, rho_prior_b);
  sd_phylo ~ normal(0, 1);
  lambda ~ beta(1, 1);

  // log |d(sd_phylo, lambda, rho) / d(V_species, w_gauss, lambda)|
  //   = log(V_species) - log(sd_phylo) + 2 * log1m(rho) + const.
  // (Stan adds the Jacobians of the <lower/upper> transforms on V_species,
  //  w_gauss and lambda itself.)
  target += log(V_species) - log(sd_phylo) + 2 * log1m(rho);

  // Equivalent to (Q' u_phylo)[k] ~ N(0, tau_k^2); the -log(1/c_k) term is
  // kept by '~' because the scale depends on parameters.
  z_phylo ~ normal(0, inv(c_pnc));

  for (n in 1:N) {
    target += zibeta_binomial2_lpmf(y[n] | trials[n], mu[n], phi, zi[n]);
  }
}

generated quantities {
  vector[N] log_lik;
  array[N] int y_rep;

  real Intercept_mu_original = Intercept_mu - dot_product(X_mu_means, b_mu);
  real Intercept_zi_original = Intercept_zi - dot_product(X_zi_means, b_zi);

  // Approximate share of species-level latent (logit-scale) variance that is
  // phylogenetically structured: sd_phylo^2 * lambda / V_species.
  real phylo_share_approx = w_gauss * lambda;

  for (n in 1:N) {
    log_lik[n] = zibeta_binomial2_lpmf(y[n] | trials[n], mu[n], phi, zi[n]);

    if (bernoulli_rng(zi[n]) == 1) {
      y_rep[n] = 0;
    } else {
      y_rep[n] = beta_binomial_rng(
        trials[n],
        mu[n] * phi,
        (1 - mu[n]) * phi
      );
    }
  }
}
