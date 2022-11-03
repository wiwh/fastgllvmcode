update_parameters <- function(parameters, gradients, alpha, learning_rate_i, ...){
  args <- list(...)
  if(!is.null(args[["median"]])) median <- args[["median"]] else median <- F

  # update step
  for (par in names(parameters)) {
    if(is.null(gradients[[par]])) next()
    if(par %in% c("Z")) {
      # we replace the previous value with the new one  THIS IS IMPORTANT
      parameters[[par]] <- parameters[[par]] - gradients[[par]]
    } else {
      # we do one step
      if (median) {
        step <- alpha * learning_rate_i * gradients[[par]]
        id_too_big <- abs(step) > median
        step[id_too_big] <- median * sign(step[id_too_big])
        parameters[[par]] <- parameters[[par]] - step
      } else {
        parameters[[par]] <- parameters[[par]] - alpha * learning_rate_i * gradients[[par]]
      }
    }
  }
  parameters <- check_update_parameters(parameters)
}

check_update_parameters <- function(parameters){
  if (any(parameters$phi < 1e-2)) {
    warning("Parameter phi too small or negative, set to 1e-4 instead.")
    parameters$phi[parameters$phi < 1e-2] <- 1e-2
  }
  if(any(abs(parameters$A) > 10)) {
    parameters$A[parameters$A > 10] <- 10
    parameters$A[parameters$A < -10] <- -10
  }
  if(!is.null(parameters$B) && any(abs(parameters$B) > 10)) {
    parameters$B[parameters$B > 10] <- 10
    parameters$B[parameters$B < -10] <- -10
  }
  # parameters <- recenter(parameters, 1) # DO NOT RESCALE! BIASES BINOMIAl p=4, q=1, n=10000
  # warning("recentered parameters")
  parameters
}


# Check that the parameters wiggle around some value
check_convergence <- function(history){
  hist_len <- nrow(history[[1]])
  # the last 30% must be wiggly
  id_plus  <- floor(hist_len/3):hist_len
  id_minus <- id_plus - 1


  diff <- sapply(names(history), function(par){
    if(par != "Z"){
      colMeans(sign(history[[par]][id_plus,, drop=F] - history[[par]][id_minus,, drop=F]))
      # re-scale
    } else {
      0
    }
  })
  diff
}
