compute_parameters_update <- function(fg, step_size, hessian, grad_batch, grad_simul, controls=list()){
  # update for AB
  if (!is.null(hessian)) {
    dAB <- mult_invHessian_dAB(grad_batch$AB - grad_simul$AB, hessian_AB = hessian)
  } else {
    dAB <-  grad_simul$AB - grad_batch$AB
  }
  dAB <- AB_separate(dAB, fg$dimensions)
  dA <- dAB$A
  dB <- dAB$B

  # update for phi
  dphi <- grad_simul$phi - grad_batch$phi

  # update for covZ
  dcovZ <- fg$parameters$covZ - grad_simul$covZ

  # updating all parameters
  fg$parameters$A <- fg$parameters$A - step_size * trim(dA, controls$trim)
  if (!is.null(fg$parameters$B)) {
    fg$parameters$B <- fg$parameters$B - step_size * trim(dB, controls$trim)
  }

  fg$parameters$phi <- fg$parameters$phi - step_size * .1 * trim(dphi, controls$trim)

  fg$parameters$covZ <- fg$parameters$covZ - min(step_size, .5) * (dcovZ) # TODO: check if it's ok to do that... I think so... but this may add some dependence

  fg$parameters <- check_update_parameters(fg$parameters)
  fg$parameters
}

