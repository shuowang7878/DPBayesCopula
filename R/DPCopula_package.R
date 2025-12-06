#' DPCopula: Differentially Private Copula Correlation Estimation
#'
#' The `DPCopula` package provides differentially private estimators for
#' Gaussian copula correlation matrices. It implements:
#'
#' \enumerate{
#'   \item A Bayesian noise-aware (Bayes-NA) estimator that explicitly accounts
#'   for differential privacy noise through a Bayesian model fitted in Stan.
#'   \item A maximum-likelihood noise-naive (MLE-NN) estimator that treats
#'   range-preserving noisy counts as if they were true statistics.
#' }
#'
#' Both methods rely on pairwise high--high counts based on median
#' binarization with lexicographic tie-breaking, together with
#' differentially private geometric or range-preserving mechanisms.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{bayes_noiseaware()}}{Bayesian noise-aware DP Gaussian copula
#'   estimator via Stan.}
#'   \item{\code{mle_noise_naive()}}{MLE noise-naive DP Gaussian copula
#'   estimator using range-preserving mechanisms.}
#'   \item{\code{range_preserving_mechanism()}}{Range-preserving DP mechanisms
#'   (RTG, TGM, BTGM) for counts.}
#' }
#'
#' @importFrom stats rnorm uniroot
#' @importFrom utils combn
#' @importFrom Matrix nearPD
#' @importFrom rstan sampling extract stan_model
#' @keywords internal
"_PACKAGE"
