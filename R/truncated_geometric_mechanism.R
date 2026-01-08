#' Truncated geometric mechanisms for counts
#'
#' Applies one of three truncated geometric mechanisms to a numeric vector of
#' counts. Each coordinate \code{i} has its own valid integer range
#' \code{[lower_i, upper_i]}, sensitivity \code{Delta_i}, and privacy budget
#' \code{epsilon_i}. 
#' 
#' The following mechanisms are implemented:
#' \describe{
#'   \item{\code{"RTG"}}{Renormalized truncated geometric mechanism. A renormalized
#'     geometric distribution on \code{lower_i,...,upper_i} is used, with the
#'     privacy-adjusted geometric parameter obtained by solving a root
#'     equation.}
#'
#'   \item{\code{"TGM"}}{Truncated geometric mechanism. Noise is first
#'     added to each count using a geometric mechanism, and
#'     the resulting value is truncated to the valid range
#'     \code{lower_i,...,upper_i}.}
#'
#'   \item{\code{"BTGM"}}{Bayesian-inspired truncated geometric mechanism.
#'     Noise is first added to each count using a geometric mechanism,
#'     followed by a Bayesian-type posterior expectation
#'     estimator under a uniform prior on \code{lower_i,...,upper_i}.
#'     The user may choose whether to round the posterior expectation.}
#' }
#'
#' @param counts Numeric vector of counts (length \code{d}). Each element must
#'   lie within its corresponding \code{[lower_i, upper_i]}.
#' @param n Optional integer; total sample size used only to set default bounds
#'   when \code{lower} or \code{upper} is missing.
#' @param total_epsilon Optional positive scalar specifying the total privacy
#'   budget for all coordinates combined. If supplied together with
#'   \code{epsilon}, the sum of \code{epsilon} must equal \code{total_epsilon}.
#' @param epsilon Optional positive numeric specifying per-coordinate privacy
#'   budgets. May be a scalar (replicated to length \code{d}) or a
#'   length-\code{d} vector. If \code{epsilon} is not supplied but
#'   \code{total_epsilon} is provided, then the default per-coordinate
#'   allocation is \code{total_epsilon / d}.
#' @param sensitivity Positive numeric; scalar or length-\code{d} vector;
#'   elementwise \eqn{\ell_1}-sensitivity.
#' @param method Character; one of \code{"RTG"}, \code{"TGM"}, \code{"BTGM"}.
#' @param lower,upper Integer bounds; scalar or length-\code{d} vectors. If both
#'   are missing, they default to \code{0} and \code{n}, respectively.
#' @param round_output Logical; if \code{TRUE} and \code{method = "BTGM"},
#'   rounds the posterior expectation before truncation.
#'
#' @return A numeric vector of perturbed counts with the same length and names
#'   as \code{counts}.
#'
#' @examples
#' counts <- c(10, 15, 20)
#' noisy  <- truncated_geometric_mechanism(
#'   counts        = counts,
#'   n             = 50,
#'   total_epsilon = 1,
#'   sensitivity   = 1,
#'   method        = "BTGM",
#'   round_output  = TRUE
#' )
#' noisy
#'
#' @export
truncated_geometric_mechanism <- function(counts,
                                          n = NULL,
                                          total_epsilon = NULL,
                                          epsilon = NULL,
                                          sensitivity,
                                          method = c("RTG", "TGM", "BTGM"),
                                          lower = NULL,
                                          upper = NULL,
                                          round_output = FALSE) {
  ## Basic checks on counts
  if (!is.numeric(counts)) {
    stop("`counts` must be a numeric vector.", call. = FALSE)
  }
  d <- length(counts)
  if (d == 0L) {
    return(counts)
  }
  
  method <- match.arg(method)
  
  ## Sensitivity checks and broadcasting
  if (!is.numeric(sensitivity) || any(sensitivity <= 0) || length(sensitivity) < 1L) {
    stop("`sensitivity` must be positive numeric (scalar or vector).", call. = FALSE)
  }
  sens_vec <- if (length(sensitivity) == 1L) rep(sensitivity, d) else sensitivity
  if (length(sens_vec) != d) {
    stop("`sensitivity` must be scalar or length equal to `counts`.", call. = FALSE)
  }
  
  ## Determine per-coordinate bounds
  if (is.null(lower) || is.null(upper)) {
    if (is.null(n)) {
      stop("`n` must be supplied when `lower` or `upper` is missing.", call. = FALSE)
    }
    if (!is.numeric(n) || length(n) != 1L || n <= 0) {
      stop("`n` must be a positive scalar numeric.", call. = FALSE)
    }
    lower <- rep.int(0L, d)
    upper <- rep.int(as.integer(n), d)
  } else {
    lower <- as.integer(lower)
    upper <- as.integer(upper)
    if (length(lower) == 1L) lower <- rep.int(lower, d)
    if (length(upper) == 1L) upper <- rep.int(upper, d)
    if (length(lower) != d || length(upper) != d) {
      stop("`lower` and `upper` must be scalar or length equal to `counts`.", call. = FALSE)
    }
  }
  
  if (any(lower > upper)) {
    stop("Each `lower[i]` must be <= `upper[i]`.", call. = FALSE)
  }
  
  ## Counts must lie in their corresponding ranges
  if (any(counts < lower | counts > upper, na.rm = TRUE)) {
    stop("Each `counts[i]` must lie inside [lower[i], upper[i]].", call. = FALSE)
  }
  
  ## Epsilon / total_epsilon handling
  if (is.null(epsilon) && is.null(total_epsilon)) {
    stop("Either `total_epsilon` or `epsilon` must be supplied.", call. = FALSE)
  }
  
  if (!is.null(total_epsilon)) {
    if (!is.numeric(total_epsilon) ||
        length(total_epsilon) != 1L ||
        total_epsilon <= 0) {
      stop("`total_epsilon` must be a positive scalar numeric.", call. = FALSE)
    }
  }
  
  eps_vec <- NULL
  
  if (!is.null(epsilon)) {
    if (!is.numeric(epsilon) || any(epsilon <= 0) || length(epsilon) < 1L) {
      stop("`epsilon` must be positive numeric (scalar or vector).", call. = FALSE)
    }
    eps_vec <- if (length(epsilon) == 1L) rep(epsilon, d) else epsilon
    if (length(eps_vec) != d) {
      stop("`epsilon` must be scalar or length equal to `counts`.", call. = FALSE)
    }
    
    if (!is.null(total_epsilon)) {
      total_from_vec <- sum(eps_vec)
      tol <- 1e-8 * max(1, abs(total_epsilon))
      if (abs(total_from_vec - total_epsilon) > tol) {
        stop("The sum of `epsilon` does not match `total_epsilon` (within tolerance).",
             call. = FALSE)
      }
    }
  } else {
    ## epsilon is NULL, derive per-coordinate budgets from total_epsilon
    eps_vec <- rep(total_epsilon / d, d)
  }
  
  out <- counts
  
  ## RTG: precompute per-coordinate epsilon_geom_i
  eps_geom_vec <- NULL
  if (method == "RTG") {
    eps_geom_vec <- numeric(d)
    
    for (i in seq_len(d)) {
      L_i <- upper[i] - lower[i]
      
      if (L_i <= 0) {
        # degenerate interval (single point): no noise needed
        eps_geom_vec[i] <- eps_vec[i]
        next
      }
      
      g_i <- function(eps_geom) {
        alpha  <- exp(-eps_geom / sens_vec[i])
        d_star <- min(sens_vec[i], ceiling(L_i / 2))
        num <- 1 + alpha - alpha^(d_star + 1) - alpha^(L_i + 1 - d_star)
        den <- 1 - alpha^(L_i + 1)
        num / den
      }
      
      f_i <- function(eps_geom) {
        eps_geom + log(g_i(eps_geom)) - eps_vec[i]
      }
      
      root <- stats::uniroot(f_i,
                             interval = c(1e-6, eps_vec[i]),
                             tol = 1e-10)$root
      eps_geom_vec[i] <- root
    }
  }
  
  ## Per-coordinate noisy count generation
  for (k in seq_len(d)) {
    x       <- counts[k]
    low     <- lower[k]
    up      <- upper[k]
    Delta_k <- sens_vec[k]
    eps_k   <- eps_vec[k]
    vals_k  <- low:up
    
    if (method == "RTG") {
      logp <- -(eps_geom_vec[k] / Delta_k) * abs(vals_k - x)
      p    <- exp(logp - .logSumExp(logp))
      out[k] <- sample(vals_k, size = 1L, prob = p)
      
    } else if (method == "TGM") {
      noisy <- x + .geometric_mechanism(eps_k, Delta_k)
      out[k] <- max(low, min(up, noisy))
      
    } else if (method == "BTGM") {
      noisy <- x + .geometric_mechanism(eps_k, Delta_k)
      logw  <- -(eps_k / Delta_k) * abs(vals_k - noisy)
      a1    <- .logSumExp(logw)
      w_norm <- exp(logw - a1)
      m <- sum(vals_k * w_norm)
      if (isTRUE(round_output)) {
        m <- round(m)
      }
      out[k] <- m
    }
  }
  
  out
}

