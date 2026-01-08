#' Pairwise high–high counts using lexicographic tie-breaking
#'
#' For each variable, an \eqn{n \times p} matrix of i.i.d. N(0,1) keys is drawn,
#' independent of the data. Observations are ordered lexicographically by
#' \code{(X[, j], key[, j])}. The first \code{ceiling(n/2)} rows are assigned
#' to the `"low"` group and the remainder to `"high"`. For every variable pair
#' \code{i < j}, the function counts the number of observations that fall into
#' the `"high"` group for both variables.
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


#' Reconstruct a symmetric matrix from named pairwise values
#'
#' Given a named numeric vector whose names encode variable pairs in the form
#' \code{"var1 vs var2"}, this function reconstructs the corresponding
#' symmetric \eqn{p \times p} matrix. The diagonal entries are set to
#' \code{diag_value}.
#'
#' @param pair_values Named numeric vector. Each name must be of the form
#'   \code{"var1 vs var2"}, where \code{var1} and \code{var2} are variable
#'   names. The order of entries is arbitrary.
#' @param diag_value Numeric scalar used to fill the diagonal. Defaults to
#'   \code{NA_real_}. For correlation matrices, \code{diag_value = 1} is
#'   typical.
#'
#' @return A symmetric numeric matrix with row and column names inferred
#'   from \code{pair_values}, and diagonal equal to \code{diag_value}.
#'
#' @keywords internal
.reconstruct_matrix_from_pairs <- function(pair_values,
                                           diag_value = NA_real_) {
  pair_values <- unlist(pair_values)
  
  if (is.null(names(pair_values))) {
    stop("`pair_values` must be a named vector with names like 'var1 vs var2'.",
         call. = FALSE)
  }
  
  ## Split names "var1 vs var2" into variable names
  name_split <- strsplit(names(pair_values), " vs ", fixed = TRUE)
  var_names  <- unique(unlist(name_split))
  p          <- length(var_names)
  
  ## Initialize matrix
  mat <- matrix(NA_real_, nrow = p, ncol = p,
                dimnames = list(var_names, var_names))
  
  ## Fill upper and lower triangles
  for (i in seq_along(pair_values)) {
    pair <- name_split[[i]]
    if (length(pair) != 2L) {
      stop("Each name in `pair_values` must contain exactly one ' vs ' separator.",
           call. = FALSE)
    }
    v1 <- pair[1L]
    v2 <- pair[2L]
    mat[v1, v2] <- pair_values[i]
    mat[v2, v1] <- pair_values[i]
  }
  
  ## Set diagonal
  diag(mat) <- diag_value
  
  mat
}

