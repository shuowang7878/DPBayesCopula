data {
  int<lower=1> p;   // Number of variables
  int<lower=1> N;   // Upper bound on the feasible true statistic
  real<lower=0> epsilon;  // Privacy budget for differential privacy
  real<lower=0> Delta_f;  // Sensitivity of the function
  matrix[p, p] noisy_counts; // Noisy contingency table counts
  matrix[p, p] count_ub;
  matrix[p, p] count_lb;
}

parameters {
  cholesky_factor_corr[p] L;  // Cholesky decomposition to ensure a valid correlation matrix
}

transformed parameters {
  matrix[p, p] R;  
  R = multiply_lower_tri_self_transpose(L);  // Compute R = L * L^T to guarantee a valid correlation matrix
  real eps_sens = epsilon / Delta_f;
  real coeff = log1m_exp(-eps_sens) - log1p_exp(-eps_sens);
}

model {
  // Prior for the correlation matrix R
  L ~ lkj_corr_cholesky(1);

  // Likelihood computation
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      real rho_ij = R[i, j];
      real log_theta_ij = log(pi() + 2 * asin(rho_ij)) - log(pi() - 2 * asin(rho_ij));
      
      int noisy = to_int(noisy_counts[i, j]);
  
      // Compute count bounds
      int ub = to_int(count_ub[i, j]);
      int lb = to_int(count_lb[i, j]);
  
      vector[ub - lb + 1] log_probs;
      vector[N + 1] hypergeo_loglik;
      for (x in 0:N) {
        hypergeo_loglik[x + 1] = 2 * lchoose(N, x) + 2 * x * log_theta_ij;
      }
      real norm_const = log_sum_exp(hypergeo_loglik);
      
      for (x in lb:ub) {
        real log_geom = -fabs(noisy - x) * eps_sens;
        log_probs[x - lb + 1] = coeff + log_geom + hypergeo_loglik[x + 1] - norm_const;
      }

  
      target += log_sum_exp(log_probs);
    }
  }

}
