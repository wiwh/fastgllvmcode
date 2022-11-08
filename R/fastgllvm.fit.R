#' Fit a fastgllvm object for fitting a gllvm model.
#'
fastgllvm.fit <- function(fg, parameters.init = NULL, controls) {

  if(is.null(fg$parameters)) {
    # compute initial values for A, B phi, Z
    initial_values <- initialize_parameters(fg)
    fg$Z <- initial_values$Z
    fg$parameters <- list(A = initial_values$A,
                          B = initial_values$B,
                          phi = initial_values$phi)
    rm("initial_values")
    fg$parameters <- initialize_additional_parameters(fg, controls$method)
  } else {
    cat("\nInitial parameters values set to those provided.")
  }

  fg$mean <- fg$linpar <- matrix(0, nrow=fg$dimensions$n, ncol=fg$dimensions$p)

  gradients <- initialize_gradients(fg$parameters)
  if (controls$hessian) {
    hessian <- simulate_hessian_AB(fg)
  } else {
    hessian <- NULL
  }

  params_hist <- list()
  if(controls$hist) params_hist <- c(params_hist, list(fg$parameters))

  # Beginning of the iterations
  # ---------------------------

  # impute for Y
  if(!is.null(fg$Miss)) {
    fg$Y[fg$Miss] <- fg$mean[fg$Miss]
  }
  signs <- fg$parameters$A * 0 + 1
  for (i in 1:controls$maxit) {
    # compute criterion after each pass...
    # only update the hessian once per pass...
    # TODO: compute_Z must accept linpar as argument! for starting values or whatever

    if (i < 20){
      step_size = controls$alpha
    } else {
      step_size = controls$alpha*20/i
    }

    # Compute the Hessian
    # fg$hessian <- update_hessian_AB(fg$hessian, compute_hessian_AB(fg$X, fg$dimensions, fg$parameters, fg$families), weight=.9)
    # Importantly, the hessian must be computed independently from the rest
    fg <- compute_mean(fg)
    fg$deviance <- mean(compute_deviance(subset(fg, 1:20)))
    cat("\ni: ", i, "dev:", fg$deviance)

    warning("the hessian must be recomputed at every batch....")



    # re-draw random batches every pass
    batches <- initialize_batches(fg$dimensions$n, controls$batch_size)

    for (batch in batches) {
      fg_batch <- subset(fg, batch)
      dAB <- compute_dAB_centered(fg_batch, method = controls$method, hessian=hessian)

      # update_a <-  trim(signs * step_size *dAB$dA, controls$trim)
      # sign_a <- sign(update_a)
      # signs[sign(signs) == sign_a] <- signs[sign(signs) == sign_a] + sign_a[sign(signs) == sign_a]
      # signs[sign(signs) != sign_a] <- -sign_a[sign(signs) !=sign_a]

      fg$parameters$A <- fg$parameters$A - trim(step_size * dAB$dA, controls$trim)
      if (!is.null(fg$parameters$B)) {
        fg$parameters$B <- fg$parameters$B - trim(step_size * dAB$dB, controls$trim)
      }
      # by the (inverse of ) Hessian
      if(controls$hist) params_hist <- c(params_hist, list(fg$parameters))
      # cat("\nSign:", signs)
    }
  }
  if(controls$hist){
    history <- sapply(names(fg$parameters), function(par_name) {
      do.call(rbind, lapply(params_hist, function(parameters_i) as.vector(parameters_i[[par_name]])))
    }, simplify=F)
  }
  fg$fit$hist <- if(controls$hist) history else NULL
  fg$Z <- compute_Z(fg, start=fg$Z)$Z
  fg
}

#
#   moving_average <- parameters
#   crit <- Inf
#
#   if(is.null(controls$learning_rate.args)) controls$learning_rate.args <- list(method="spall", rate=2, end=.1)
#   if(is.null(controls$learning_rate.args$method)) controls$learning_rate.args$method <- "spall"
#   if(is.null(controls[["learning_rate"]])) controls$learning_rate <- initialize_learning_rate(maxit=controls$maxit, learning_rate.args = controls$learning_rate.args)
#
#   for(i in 1:controls$maxit){
#     # if(i < controls$maxit/2) median <- .5 else median <- F
#     moving_average_old <- moving_average
#     gradients_old <- gradients
#     # get the gradient
#
#     gradients <- compute_gradients(fg$Y, fg$X, parameters, fg$families, fg$Miss, debiase=T)
#
#     # exponential smoothing step
#     # for (par in names(parameters)) {
#     #   gradients[[par]] <- (controls$beta * gradients[[par]] + (1-controls$beta) * gradients[[par]])
#     # }
#     parameters <- update_parameters(parameters, gradients, alpha=controls$alpha, learning_rate_i = controls$learning_rate(i), median=median)
#
#     # We track a moving average of the last 1/controls$ma iterates
#     for(k in seq_along(parameters)){
#       moving_average[[k]] <- controls$ma * moving_average[[k]] + (1 - controls$ma) * parameters[[k]]
#     }
#
#     if(i == 1 || i%%10 == 0 || i == controls$maxit){
#       crit <- compute_error(moving_average$A, moving_average_old$A, rotate=F)
#       cat("\n Iteration: ", i, " - crit: ", crit)
#     }
#     if(hist) params_hist <- c(params_hist, list(parameters))
#     if(i >= controls$minit && crit < controls$eps) break()
#   }
#
#   if(hist){
#     history <- sapply(names(parameters), function(par_name) {
#       do.call(rbind, lapply(params_hist, function(parameters_i) as.vector(parameters_i[[par_name]])))
#     }, simplify=F)
#   }
#
#   # Update the fastgllvm object
#   fg$parameters <- moving_average
#   fg$fit <- list(
#     crit = crit,
#     controls = controls,
#     hist = if(hist) history else NULL
#   )
#   fg$converged <- ifelse(crit < controls$eps, T, F)
#   fg
# }


if(0) {

  library(gmf)
  # TODO: GROS PROBLEME AVEC LA HESSIENNE?
  devtools::load_all()
  poisson  <- 0
  gaussian <- 0
  binomial <- 50
  nobs <- 1000
  q <- 2
  p <- poisson + gaussian + binomial

  intercept <- T
  k <- 1
  if(k==0 & intercept) k <- 1

  family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial))
  # set.seed(1234)
  fg <- gen_fastgllvm(nobs=nobs, p=p, q=q, k=k, family=family, intercept=intercept, miss.prob = 0, scale=1)

  # set.seed(1304)
  fit1 <- fastgllvm(fg$Y, q = q, family=family, hist=T, method="full", batch_size=nobs, trim=.5, intercept=intercept, alpha=1, hessian=T, maxit=50)
  plot(fit1)

  fit.gmf <- gmf(fg$Y, fg$X, family=binomial(), p=1)

  fitm <- compute_mean(fit1)
  plot(fg$Y, fitm$mean)
  points(fg$Y, fit.gmf$fit, col=2)
  abline(0,1,col=2)

  fit2 <- fitm
  fit2$mean <- fit.gmf$fit

  mean(compute_deviance(fitm))
  mean(compute_deviance(fit2))




  # fit <- fit2

  fit <- fit1
  plot(fit)
  plot(fg$parameters$A, psych::Procrustes(fit$parameters$A, fg$parameters$A)$loadings, col=3)
  points(fg$parameters$A, psych::Procrustes(fit.gmf$v, fg$parameters$A)$loadings, col=4)
  points(fg$parameters$B, t(fit.gmf$beta), col=4)
  points(fg$parameters$B, fit$parameters$B, pch=2, col=3)
  abline(0,1,col=2)
  plot(fg$Z, fit$Z); abline(0,-1,col=2); abline(0,1,col=2)


  library(ltm)
  if(q==1) fit.ltm <- ltm(fg$Y ~ z1)
  if(q==2) fit.ltm <- ltm(fg$Y ~ z1 + z2)

  plot(fg$parameters$A, psych::Procrustes(fit$parameters$A, fg$parameters$A)$loadings, col=3)
  points(fg$parameters$B, fit$parameters$B, pch=2, col=3)
  points(fg$parameters$A, psych::Procrustes(fit.ltm$coefficients[,2:(1+q)], fg$parameters$A)$loadings)#, xlim=c(-5,5), ylim=c(-5,5))
  points(fg$parameters$B, fit.ltm$coefficients[,1], pch=2)
  # points(fg$parameters$A, psych::Procrustes(fit.full$parameters$A, fg$parameters$A)$loadings, col=2)
  # points(fg$parameters$B, fit.full$parameters$B, pch=2, col=2)
  legend(x="bottomright",pch=1,  legend=c("ltm","full", "simple"), col=1:3)
  abline(0,1,col=2)
  abline(0,-1,col=2)
  compute_error(fit.ltm$coefficients[,2:(1+q), drop=F], fg$parameters$A, rotate=T)
  # compute_error(fit$parameters$A, fg$parameters$A, rotate=T)
  compute_error(fit$parameters$A, fg$parameters$A, rotate=T)
  # compute_error(fit.full$parameters$A, fg$parameters$A, rotate=T)

  # plot(fg$Z, fit$Z)
  # set.seed(1304)
  plot(fit)
}
