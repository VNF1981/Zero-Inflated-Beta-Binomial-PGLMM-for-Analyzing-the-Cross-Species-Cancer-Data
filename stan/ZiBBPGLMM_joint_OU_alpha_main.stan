// ZiBB-PGLMM, Ornstein-Uhlenbeck (alpha) -- exact reparameterisation CANDIDATE
// (October 2026).  SAME model, likelihood, priors, data block and OU
// correlation function as the previous file.  Only sampling coordinates change:
//
// (1) Variance coordinates, as in the lambda/BM files:
//         V_species = sd_phylo^2 + kappa / phi,  w_gauss = sd_phylo^2 / V_species,
//     priors kept on sd_phylo and rho, log-Jacobian
//     log(V) - log(sd_phylo) + 2 log(1 - rho) added.
//
// (2) Latent field: the OU correlation depends on alpha, so no fixed
//     eigenbasis exists; the Cholesky factor L(alpha) is kept.  Partial
//     non-centring is applied in Cholesky coordinates:
//         u = sd_phylo * L(alpha) * (c .* z_phylo),  z_phylo[i] ~ N(0, 1/c_i^2)
//     <=> u ~ N(0, sd_phylo^2 A_ou(alpha)) exactly, for ANY positive c.
//     c_i = (1 + sd_phylo^2 * J_i)^(-1/2), with J_i = sum_j L[j,i]^2 D_j the
//     approximate data information about Cholesky coordinate i and
//     D_j = trials_j pbar(1-pbar) / (1 + (trials_j - 1) rho).
//
// UNTESTED on real data: unlike lambda/BM, this does not remove the
// alpha-dependent rotation of the Cholesky coordinates. Validate before use.

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

  // OU correlation for a tree expressed in root-to-tip/shared-root time.
  // alpha is dimensionless because the R driver scales tree height to 1.
  // The alpha -> 0 limit is the Brownian correlation matrix.
  matrix ou_corr_from_shared_time(matrix T_shared, vector tip_heights, real alpha) {
    int N = rows(T_shared);
    matrix[N, N] C;

    if (alpha < 1e-6) {
      for (i in 1:N) {
        for (j in i:N) {
          real val = T_shared[i, j] / sqrt(tip_heights[i] * tip_heights[j]);
          C[i, j] = val;
          C[j, i] = val;
        }
      }
    } else {
      for (i in 1:N) {
        for (j in i:N) {
          real tij = T_shared[i, j];
          real ti = tip_heights[i];
          real tj = tip_heights[j];

          real numerator =
            exp(-alpha * (ti + tj - 2.0 * tij)) *
            (-expm1(-2.0 * alpha * tij));

          real denominator = sqrt(
            (-expm1(-2.0 * alpha * ti)) *
            (-expm1(-2.0 * alpha * tj))
          );

          real val = numerator / denominator;
          C[i, j] = val;
          C[j, i] = val;
        }
      }
    }

    for (i in 1:N) {
      C[i, i] = 1.0;
    }
    return C;
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

  matrix[N, N] T_shared;
  vector<lower=0>[N] tip_heights;

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
}

parameters {
  real Intercept_mu;
  vector[K_mu] b_mu;

  real Intercept_zi;
  vector[K_zi] b_zi;

  real<lower=0> V_species;            // sd_phylo^2 + kappa / phi
  real<lower=0, upper=1> w_gauss;     // sd_phylo^2 / V_species

  vector[N] z_phylo;

  real<lower=0> alpha;
}

transformed parameters {
  real<lower=0> sd_phylo = sqrt(w_gauss * V_species);
  real<lower=0> phi = kappa / ((1 - w_gauss) * V_species);
  real<lower=0, upper=1> rho = 1 / (1 + phi);

  matrix[N, N] A_ou;
  matrix[N, N] L_A_ou;
  vector[N] c_pnc;
  vector[N] u_phylo;
  vector[N] eta_mu;
  vector[N] eta_zi;
  vector[N] mu;
  vector[N] zi;

  A_ou = ou_corr_from_shared_time(T_shared, tip_heights, alpha);
  for (i in 1:N) {
    A_ou[i, i] += 1e-10;
  }

  L_A_ou = cholesky_decompose(A_ou);
  {
    vector[N] D = info_num ./ (1.0 + trials_m1 * rho);
    vector[N] J = (L_A_ou .* L_A_ou)' * D;
    c_pnc = inv_sqrt(1.0 + square(sd_phylo) * J);
  }
  u_phylo = sd_phylo * (L_A_ou * (c_pnc .* z_phylo));

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

  // Original priors on the original parameters (unchanged).
  rho ~ beta(rho_prior_a, rho_prior_b);
  sd_phylo ~ normal(0, 1);

  // log |d(sd_phylo, rho) / d(V_species, w_gauss)| + const.
  target += log(V_species) - log(sd_phylo) + 2 * log1m(rho);

  // Equivalent to the previous z ~ N(0, I) with u = sd_phylo * L * z.
  z_phylo ~ normal(0, inv(c_pnc));

  // Tree height is scaled to 1 in R, so alpha is on a comparable scale
  // across all model-specific species subsets.
  alpha ~ exponential(1);

  for (n in 1:N) {
    target += zibeta_binomial2_lpmf(y[n] | trials[n], mu[n], phi, zi[n]);
  }
}

generated quantities {
  vector[N] log_lik;
  array[N] int y_rep;

  real Intercept_mu_original = Intercept_mu - dot_product(X_mu_means, b_mu);
  real Intercept_zi_original = Intercept_zi - dot_product(X_zi_means, b_zi);
  real phylo_share_approx = w_gauss;

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
