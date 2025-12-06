#' Bayesian noise-aware differentially private Gaussian copula (Bayes-NA)
#'
#' This function implements the Bayesian noise-aware (Bayes-NA) estimator for
#' a Gaussian copula correlation matrix under differential privacy. It first
#' constructs differentially private pairwise high–high counts based on
#' median thresholds, then fits a Bayesian latent Gaussian model in Stan to
#' recover the correlation matrix.
#'
#' For each pair \eqn{(i, j)}, the sufficient statistic is
#' \deqn{
#'   T_{ij} = \sum_{s=1}^n \mathbf{1}\{X_{s,i} \ge \mathrm{med}(X_i),
#'                                    X_{s,j} \ge \mathrm{med}(X_j)\},
#' }
#' where \eqn{\mathrm{med}(X_i)} is the sample median of column \eqn{i}.
#' Ties at the median are included in the "\eqn{\ge} median" group.
#'
#' @param data A numeric data.frame or matrix (rows = observations, columns = variables).
#'   No missing values allowed.
#' @param epsilon Global privacy budget. The per-pair budget is
#'   \code{epsilon / choose(p, 2)}.
#' @param iter_sampling Number of post-warmup sampling iterations. Defaults to 100.
#' @param iter_warmup Number of warmup iterations. Defaults to 100.
#' @param seed Optional random seed for reproducibility (affects DP noise).
#' @param chains Number of Stan MCMC chains (default 1).
#' @param cores Number of CPU cores for Stan (default 1).
#'
#' @return A list containing:
#'   \describe{
#'     \item{mean}{Posterior mean correlation matrix, i.e. estimated 
#'     differentially private \eqn{p \times p} correlation matrix.}
#'     \item{sd}{Posterior standard deviation matrix.}
#'     \item{ci_lower}{2.5\% elementwise posterior quantiles.}
#'     \item{ci_upper}{97.5\% elementwise posterior quantiles.}
#'     \item{samples}{Posterior samples of the correlation matrix,
#'     an array of dimension \code{(iter_sampling * chains) x p x p}.}
#'     \item{dp}{A list with differential privacy metadata, including
#'       \code{epsilon}, \code{epsilon_per_pair}, \code{sensitivity}[].}
#'   }
#'
#' @examples
#' \dontrun{
#'   set.seed(123)
#'   X <- matrix(rnorm(300), ncol = 3)
#'   fit <- bayes_noiseaware(X, epsilon = 1, iter_sampling = 500, iter_warmup = 500)
#'   fit$mean
#' }
#'
#' @export
bayes_noiseaware <- function(
    data,
    epsilon,
    iter_sampling = 100,
    iter_warmup   = 100,
    seed          = NULL,
    chains        = 1,
    cores         = 1
) {
  ## Basic input validation
  if (!is.data.frame(data) && !is.matrix(data))
    stop("'data' must be a data.frame or matrix.", call. = FALSE)
  data <- as.data.frame(data)
  
  if (any(!sapply(data, is.numeric)))
    stop("All columns of 'data' must be numeric.", call. = FALSE)
  if (anyNA(data))
    stop("Missing values are not allowed.", call. = FALSE)
  
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(data)
  p <- ncol(data)
  if (p < 2) stop("At least two variables are required.", call. = FALSE)
  if (n < 2) stop("At least two observations are required.", call. = FALSE)
  if (!is.numeric(epsilon) || length(epsilon) != 1L || epsilon <= 0) {
    stop("'epsilon' must be a positive scalar numeric.", call. = FALSE)
  }
  
  ## Privacy budget per pair
  eps_pair <- epsilon / choose(p, 2)
  ## Global sensitivity is 1 (proof can be found in the paper)
  sensitivity <- 1
  
  ## Load the package's built-in Stan model
  model <- .get_bayesna_model()
  
  ## Pairwise high–high counts based on medians
  pair_counts <- .pairwise_median_counts(data)
  
  ## Apply geometric DP noise to each T_ij
  noisy_pairvals <- pair_counts + vapply(
    seq_along(pair_counts),
    function(k) .geometric_mechanism(eps_pair, sensitivity),
    numeric(1L)
  )
  
  ## Reconstruct p×p matrix for Stan
  noisy_mat <- .reconstruct_matrix_from_pairs(noisy_pairvals, diag_value = 1)
  
  ## Stan uses -1 to mark unavailable (lower-triangular entries)
  noisy_mat[lower.tri(noisy_mat, diag = TRUE)] <- -1
  
  ## Compute k for the symmetric geometric mechanism such that
  ## P(|noise| > k) ≤ 10^-8 when noise ~ Geom(eps_pair / sensitivity).
  k <- ceiling(
    (log(2) - log1p(exp(-eps_pair / sensitivity)) - (-8 * log(10))) /
      (eps_pair / sensitivity) - 1
  )
  
  ## Upper/lower bounds for the unobserved true statistic
  count_ub <- pmin(noisy_mat + k, floor((n + 1) / 2))
  count_lb <- pmax(noisy_mat - k, 0)
  
  count_ub[lower.tri(count_ub, diag = TRUE)] <- -1
  count_lb[lower.tri(count_lb, diag = TRUE)] <- -1
  
  ## Stan data list
  stan_data <- list(
    p            = p,
    N            = floor((n + 1) / 2),
    epsilon      = eps_pair,
    Delta_f      = sensitivity,
    noisy_counts = noisy_mat,
    count_ub     = count_ub,
    count_lb     = count_lb
  )
  
  ## Fit Bayesian model
  fit <- rstan::sampling(
    model,
    data   = stan_data,
    warmup = iter_warmup,
    iter   = iter_warmup + iter_sampling,
    chains = chains,
    cores  = cores,
    refresh = 0
  )
  
  ## Extract posterior draws
  draws <- rstan::extract(fit, pars = "R")$R
  dimnames(draws) <- list(NULL, colnames(data), colnames(data))
  
  ## Posterior summaries
  post_mean  <- apply(draws, c(2, 3), mean)
  post_sd    <- apply(draws, c(2, 3), sd)
  ci_lower   <- apply(draws, c(2, 3), quantile, probs = 0.025)
  ci_upper   <- apply(draws, c(2, 3), quantile, probs = 0.975)
  
  dimnames(post_mean) <- list(colnames(data), colnames(data))
  dimnames(post_sd)   <- list(colnames(data), colnames(data))
  dimnames(ci_lower)  <- list(colnames(data), colnames(data))
  dimnames(ci_upper)  <- list(colnames(data), colnames(data))
  
  ## Output
  list(
    samples  = draws,
    mean     = post_mean,
    sd       = post_sd,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    dp       = list(
      epsilon          = epsilon,
      epsilon_per_pair = eps_pair,
      sensitivity      = sensitivity
    )
  )
}




#' MLE noise-naive differentially private Gaussian copula
#'
#' Implements the MLE noise-naive (MLE-NN) estimator for a Gaussian copula
#' correlation matrix under differential privacy. Each variable is binarized
#' into "low" and "high" groups using the sample median with random
#' tie-breaking so that exactly \code{ceil(n/2)} observations fall into the
#' "high" group. Pairwise high–high counts are then perturbed using one of three
#' range-preserving mechanisms (RTG, TGM, BTGM). The resulting noisy statistics
#' are treated as if they were true counts, yielding a noise-naive likelihood
#' estimator of the copula correlations.
#'
#' @param data Numeric data.frame or matrix. No missing values allowed.
#' @param epsilon Total privacy budget. The per-pair budget is
#'   \code{epsilon / choose(p, 2)}.
#' @param method One of \code{"RTG"}, \code{"TGM"}, or \code{"BTGM"}.
#' @param seed Optional random seed (affects DP noise and tie-breaking).
#'
#' @return A list with components:
#'   \describe{
#'     \item{estimation}{Estimated \eqn{p \times p} correlation matrix
#'     (positive-semidefinite).}
#'     \item{dp}{List containing \code{epsilon}, \code{epsilon_per_pair},
#'       and \code{sensitivity}.}
#'   }
#'
#' @export
mle_noise_naive <- function(
    data,
    epsilon,
    method = c("RTG", "TGM", "BTGM"),
    seed   = NULL
) {
  ## Data validation
  if (!is.data.frame(data) && !is.matrix(data))
    stop("'data' must be a data.frame or matrix.", call. = FALSE)
  
  data <- as.data.frame(data)
  
  if (any(!vapply(data, is.numeric, logical(1L))))
    stop("All columns of 'data' must be numeric.", call. = FALSE)
  
  if (anyNA(data))
    stop("Missing values are not allowed in 'data'.", call. = FALSE)
  
  if (!is.numeric(epsilon) || length(epsilon) != 1L || epsilon <= 0)
    stop("'epsilon' must be a positive scalar.", call. = FALSE)
  
  method <- match.arg(method)
  if (!is.null(seed)) set.seed(seed)
  
  ## Basic dimensions
  n <- nrow(data)
  p <- ncol(data)
  
  if (p < 2L)
    stop("At least two variables are required.", call. = FALSE)
  
  if (n < 2L)
    stop("At least two observations are required.", call. = FALSE)
  
  ## Privacy budget for each pair of variables
  eps_pair <- epsilon / choose(p, 2)
  ## Global sensitivity of each high–high count
  sensitivity <- 1
  
  ## Compute pairwise high–high counts with lexicographic tie-breaking
  pair_counts <- .pairwise_median_counts(data)
  
  ## Number of observations assigned to the “high” group
  n_half <- floor((n + 1L) / 2L)
  
  ## Apply range-preserving DP noise to the vector of counts
  noisy_counts <- range_preserving_mechanism(
    counts       = pair_counts,
    n            = n_half,
    epsilon      = eps_pair,
    sensitivity  = sensitivity,
    method       = method
  )
  
  k_grid <- 0:n_half
  
  ## Computes the MLE correlation for a single noisy count
  single_mle <- function(k) {
    if (k <= 0)  return(-1)
    if (k >= n_half) return(1)
    
    obj <- function(r) {
      log_theta <- log(pi + 2 * asin(r)) - log(pi - 2 * asin(r))
      vals <- 2 * lchoose(n_half, k_grid) + 2 * k_grid * log_theta
      norm_vals <- vals - .logSumExp(vals)
      expected_k <- sum(k_grid * exp(norm_vals))
      expected_k - k
    }
    
    out <- try(
      stats::uniroot(obj, lower = -1 + 1e-8, upper = 1 - 1e-8),
      silent = TRUE
    )
    
    if (inherits(out, "try-error")) NA_real_ else out$root
  }
  
  ## Vector of pairwise correlations
  raw_corr <- vapply(noisy_counts, single_mle, numeric(1L))
  names(raw_corr) <- names(noisy_counts)
  
  ## Reconstruct full correlation matrix
  corr_mat <- .reconstruct_matrix_from_pairs(raw_corr, diag_value = 1)
  diag(corr_mat) <- 1
  
  ## Enforce positive semidefiniteness
  eigvals <- eigen(corr_mat, only.values = TRUE)$values
  if (any(eigvals < 0)) {
    corr_mat <- as.matrix(Matrix::nearPD(corr_mat, corr = TRUE)$mat)
  }
  
  ## Final output
  list(
    estimation = corr_mat,
    dp = list(
      epsilon          = epsilon,
      epsilon_per_pair = eps_pair,
      sensitivity      = sensitivity
    )
  )
}
