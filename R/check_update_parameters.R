check_update_parameters <- function(parameters){
  if (any(parameters$phi < 1e-4)) {
    warning("Parameter phi too small or negative, set to 1e-4 instead.")
    parameters$phi[parameters$phi < 1e-4] <- 1e-4
  }
  # parameters <- recenter(parameters, 1) # DO NOT RESCALE! BIASES BINOMIAl p=4, q=1, n=10000
  # warning("recentered parameters")
  parameters
}
