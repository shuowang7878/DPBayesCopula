.bayesna_stan_env <- new.env(parent = emptyenv())

.get_bayesna_model <- function() {
  if (!is.null(.bayesna_stan_env$model)) {
    return(.bayesna_stan_env$model)
  }
  
  stan_file <- system.file("stan", "DPCopula_BayesNA.stan", package = "DPCopula")
  if (stan_file == "")
    stop("Stan model file 'DPCopula_BayesNA.stan' not found.")
  
  message("Compiling Stan model DPCopula_BayesNA.stan ...")
  .bayesna_stan_env$model <- rstan::stan_model(stan_file)
  
  .bayesna_stan_env$model
}
