.bayesna_stan_env <- new.env(parent = emptyenv())

.get_bayesna_model <- function() {
  if (!is.null(.bayesna_stan_env$model)) {
    return(.bayesna_stan_env$model)
  }
  
  stan_file <- system.file("stan", "BayesNA.stan", package = "DPBayesCopula")
  if (stan_file == "")
    stop("Stan model file 'BayesNA.stan' not found.")
  
  message("Compiling Stan model BayesNA.stan ...")
  .bayesna_stan_env$model <- rstan::stan_model(stan_file)
  
  .bayesna_stan_env$model
}
