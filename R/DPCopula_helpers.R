#' Pairwise high–high counts using lexicographic tie-breaking
#'
#' For each variable, an \eqn{n \times p} matrix of i.i.d. N(0,1) keys is drawn,
#' independent of the data. Observations are ordered lexicographically by
#' \code{(X[, j], key[, j])}. The first \code{ceiling(n/2)} rows are assigned
#' to the `"low"` group and the remainder to `"high"`. For every variable pair
#' \code{i < j}, the function counts the number of observations that fall into
#' the `"high"` group for both variables. This produces a set of pairwise
#' high–high counts used as the statistics for the DP Bayesian model.
#'
#' @param data Numeric data.frame or matrix.
#' @return Named numeric vector of pairwise high–high counts.
#' @keywords internal
.pairwise_median_counts <- function(data) {
  data <- as.data.frame(data)
  n <- nrow(data)
  p <- ncol(data)
  
  ## Generate lexicographic keys
  keys <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  colnames(keys) <- colnames(data)
  
  ## Determine high/low groups
  group_df <- as.data.frame(matrix(NA_character_, nrow = n, ncol = p))
  colnames(group_df) <- colnames(data)
  
  for (j in seq_len(p)) {
    ord <- order(data[[j]], keys[, j])
    n_low <- floor(n / 2)
    
    grp <- character(n)
    grp[ord[seq_len(n_low)]] <- "low"
    grp[ord[(n_low + 1L):n]] <- "high"
    
    group_df[[j]] <- grp
  }
  
  ## Compute high–high counts for all pairs
  pairs <- utils::combn(seq_len(p), 2, simplify = FALSE)
  coln  <- colnames(data)
  
  counts <- vapply(
    pairs,
    function(idx) {
      i <- idx[1L]
      j <- idx[2L]
      sum(group_df[[i]] == "high" & group_df[[j]] == "high")
    },
    numeric(1L)
  )
  
  names(counts) <- vapply(
    pairs,
    function(idx) paste(coln[idx[1L]], "vs", coln[idx[2L]]),
    character(1L)
  )
  
  counts
}



#' Two-sided geometric mechanism
#'
#' Draws discrete DP noise as the difference of two independent geometric random
#' variables. Used to perturb count data.
#'
#' @param epsilon Privacy budget.
#' @param sensitivity Global L1 sensitivity.
#' @return Integer DP noise.
#' @keywords internal
.geometric_mechanism <- function(epsilon, sensitivity) {
  lambda <- exp(-epsilon / sensitivity)
  stats::rgeom(1, 1 - lambda) - stats::rgeom(1, 1 - lambda)
}

#' Numerically stable log-sum-exp
#'
#' Computes \eqn{\log\left(\sum_i e^{x_i}\right)} in a numerically stable way.
#' This avoids overflow when \code{x} contains large positive values and avoids
#' underflow when values are very negative.
#'
#' @param x A numeric vector.
#'
#' @return A numeric scalar equal to \code{log(sum(exp(x)))} computed
#'   stably.
#'
#' @keywords internal
.logSumExp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}
