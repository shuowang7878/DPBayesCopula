data {
  int<lower=1> p;   // Number of variables
  int<lower=1> N;   // Upper bound on the feasible true statistic
  real<lower=0> epsilon;  // Privacy budget for differential privacy
  real<lower=0> Delta_f;  // Sensitivity of the function
  matrix[p, p] noisy_counts; // Noisy contingency table counts
}

parameters {
  cholesky_factor_corr[p] L;  // Cholesky decomposition to ensure a valid correlation matrix
}

transformed parameters {
  matrix[p, p] R;  
  R = multiply_lower_tri_self_transpose(L);  // Compute R = L * L^T to guarantee a valid correlation matrix
}

model {
  // Prior for the correlation matrix R
  L ~ lkj_corr_cholesky(1);

  // Likelihood computation
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      real rho_ij = R[i, j];  // Extract the (i, j) correlation from R
      real log_theta_ij = log(pi() + 2 * asin(rho_ij)) - log(pi() - 2 * asin(rho_ij));

      // Compute the posterior distribution of X_ij
      vector[N+1] log_probs;
      for (X_ij in 0:N) {
        real log_geom = -abs(noisy_counts[i, j] - X_ij) * epsilon / Delta_f;  // Geometric noise
        real log_prob = 2 * X_ij * log_theta_ij;  // power term
        real log_binom = 2 * lchoose(N, X_ij);  // squared binomial coefficient
        log_probs[X_ij + 1] = log_geom + log_prob + log_binom;
      }
      
      // Use log-sum-exp trick to marginalize over X_ij
      target += log_sum_exp(log_probs);
    }
  }
}
