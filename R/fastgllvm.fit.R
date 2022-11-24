#' Fit a fastgllvm object for fitting a gllvm model.
#'
fastgllvm.fit <- function(fg, parameters.init = NULL, controls) {
  # Initialization
  # --------------

  if (!is.null(parameters.init)) {
    fg$parameters <- parameters.init
    if(is.null(fg$Z)) fg$Z <- compute_Z(fg)$Z
  } else if(is.null(fg$parameters)) {
    # compute initial values for A, B phi, Z, covZ
    fg <- initialize_parameters(fg, return_object = T)
  } else {
    cat("\nInitial parameters set to the values from the supplied fit.")
  }

  fg <- compute_mean(fg, return_object = T)
  fg$deviance <- mean(compute_deviance(fg))

  params_hist <- list()
  if(!is.null(controls$hist)) params_hist <- c(params_hist, list(c(fg$parameters, deviance=fg$deviance)))

  ma <- fg$parameters

  # Impute for Y
  if(!is.null(fg$Miss)) {
    fg$Y[fg$Miss] <- fg$mean[fg$Miss]
  }

  # Beginning of the iterations
  # ---------------------------
  for (i in 1:controls$maxit) {

    if (controls$hessian) hessian <- simulate_hessian_AB(fg)

    if (i < 20){
      step_size = controls$alpha
    } else {
      step_size = controls$alpha*20/i
    }


    warning("the hessian must be recomputed at every batch....")
    if(!is.null(controls$H.seed)) warning("The seed is beeing reset!")


    # re-draw random batches every pass
    batches <- initialize_batches(fg$dimensions$n, controls$batch_size)

    for (batch in batches) {
      fg_batch <- subset(fg, batch)

      if (!is.null(controls$H.seed)) {
        set.seed(controls$H.seed)
      }
      if (controls$H > 1) {
        sims <- lapply(1:controls$H, function(na) {
          compute_dAB_centered(fg_batch, controls=controls, hessian=hessian)
        })
        dA <- Reduce("+", lapply(sims, function(sim)sim$dA)) / length(sims)
        dB <- Reduce("+", lapply(sims, function(sim)sim$dB)) / length(sims)
        dphi <- Reduce("+", lapply(sims, function(sim)sim$dphi)) / length(sims)
        covZ <- Reduce("+", lapply(sims, function(sim)sim$covZ)) / length(sims)
      } else {
        dAB <- compute_dAB_centered(fg_batch, controls=controls, hessian=hessian)
        dA <- dAB$dA
        dB <- dAB$dB
        dphi <- dAB$dphi
        covZ <- dAB$covZ
      }

      # updating all parameters
      fg$parameters$A <- fg$parameters$A - trim(step_size * dA, controls$trim)
      if (!is.null(fg$parameters$B)) {
        fg$parameters$B <- fg$parameters$B - trim(step_size * dB, controls$trim)
      }

      fg$parameters$phi <- fg$parameters$phi - trim(step_size * dphi, controls$trim)
      fg$parameters <- check_update_parameters(fg$parameters)

      fg$parameters$covZ <- covZ

      # cat("\nSign:", signs)
    }

    fg <- compute_Z(fg, start=fg$Z, return_object = T)
    fg <- compute_mean(fg, return_object = T)  # TODO: compute mean based on ma mayhaps?

    # rescale the model
    if (controls$rescale) {
      fg_rescaled <- rescale(fg, rescale.A = T, rescale.B = T, target.cov = fg$parameters$covZ)
      fg$parameters$A <- fg$parameters$A * .9 + fg_rescaled$parameters$A * .1
      fg$parameters$B <- fg_rescaled$parameters$B * .9 + fg_rescaled$parameters$B * .1
    }


    # fg$parameters$phi <- compute_phi(fg)
    # impute for Y
    if(!is.null(fg$Miss)) {
      fg$Y[fg$Miss] <- fg$mean[fg$Miss]
    }
    ma <- update_moving_average(ma, fg$parameters, 0.9)

    fg$deviance <- mean(compute_deviance(fg))
    cat("\ni: ", i, "dev:", fg$deviance)

    if(!is.null(controls$hist)){
      if (length(params_hist) > controls$hist) params_hist[[1]] <- NULL
      params_hist <- c(params_hist, list(c(fg$parameters, deviance=fg$deviance)))
    }
  }

  fg$parameters <- ma
  fg <- compute_mean(fg, return_object = T)

  if(!is.null(controls$hist)){
    history <- sapply(names(params_hist[[1]]), function(par_name) {
      do.call(rbind, lapply(params_hist, function(parameters_i) as.vector(parameters_i[[par_name]])))
    }, simplify=F)
  }

  fg$fit$hist <- if(controls$hist) history else NULL
  fg$controls <- controls
  fg$controls$alpha <- step_size
  fg$Z <- compute_Z(fg, start=fg$Z)$Z
  fg
}

update_signs_count <- function(signs_count, signs_old, signs_new) {
  signs_equal <- signs_new == signs_old
  signs_count[signs_equal] <- signs_count[signs_equal] + 1
  signs_count[!signs_equal] <- 0
  signs_count
}

# update moving average

update_moving_average <- function(moving_average, parameters, weight) {
  sapply(names(moving_average), function(par) {
    moving_average[[par]] <- weight * moving_average[[par]] + (1 - weight) *parameters[[par]]
  }, simplify=F)
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

  devtools::load_all()
  poisson  <- 0
  gaussian <- 100
  binomial <- 0
  nobs <- 100
  q <- 2
  p <- poisson + gaussian + binomial



  intercept <- T
  k <- 1
  if(k==0 & intercept) k <- 1


  family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial))
  set.seed(14240)
  fg <- gen_fastgllvm(nobs=nobs, p=p, q=q, k=k, family=family, intercept=intercept, miss.prob = 0.4, scale=1, phi=rep(1, p))

  if(0) {
    fg$Y <- scale(fg$Y, scale=F)
    sds <- sqrt(colMeans(fg$Y^2))

    fit.fad <- fad::fad(fg$Y, factors = 10)
    fit.ffa <- ffa(fg$Y, 10, maxiter = 100, eps=1e-10, savepath = T)

    compute_FFA_error(fg$Y, fit.fad$loadings * sds, fit.fad$uniquenesses*sds^2)
    compute_FFA_error(fg$Y, fit.ffa$A, fit.ffa$Psi)

    ts.plot(fit.ffa$path$Psi)

    compute_error(fit1$loadings  * sds, fg$parameters$A) - compute_error(fit2$A, fg$parameters$A)
    compute_error(fit1$loadings  * sds, fg$parameters$A)
    compute_error(fit2$A, fg$parameters$A)

    sims <- sapply(1:10, function(na) {
      set.seed(213131+na)
      fg <- gen_fastgllvm(nobs=nobs, p=p, q=q, k=k, family=family, intercept=intercept, miss.prob = 0, scale=1, phi = runif(p, .2, .8))
      sds <- apply(fg$Y, 2, sd)
      fit1 <- fad::fad(fg$Y, factors = q)
      fit2 <- ffa(fg$Y, q, maxiter = 20, eps=1e-5, savepath = T)
      cat("\n", fit2$niter)
      compute_error(fit1$loadings  * sds, fg$parameters$A) - compute_error(fit2$A, fg$parameters$A)
    })

  hist(sims)

  mbm <- microbenchmark::microbenchmark(fad::fad(fg$Y, factors = 10), ffa(fg$Y, 10, maxiter = 100, eps=1e-5, savepath = T), times=10)

  plot(fg$parameters$A, psych::Procrustes(fit1$loadings*sds, fg$parameters$A)$loadings, col=1)
  points(fg$parameters$A, psych::Procrustes(fit2$A, fg$parameters$A)$loadings, col=2)
  abline(0,1,col=3)


  mean((fit1$uniquenesses*sds - 1)^2)
  mean((fit2$Psi - 1)^2)

  # # rescaling inside the loop: needs to add the effect of the rescaling.... anyway it only affects the main since there is no rescaling on the simulation since we take the scaling factor there!
  # # ok for rescaling!
  # set.seed(13342)
  # fit0 <- fastgllvm(fg$Y, q = q, family=family, hist=T, method="rescaled", batch_size=1000, trim=.1, intercept=intercept, alpha=.5, hessian=T, maxit=100, use_signs = F, H=1, rescale=F)
  #
  # plot(fit0)

  }

  # full, with hessian, no rescaling
  set.seed(13342)
  fit1 <- fastgllvm(fg$Y, q = q, family=family, hist=100, method="full", batch_size=100, trim=.1, intercept=intercept, alpha=.2, hessian=T, maxit=100, use_signs = F, H=1, rescale=F)
  fit1 <- update(fit1, H=10, maxit=10)
  plot(fit1)
  # clear winner in poisson
  # full, with rescaling outside the loop
  set.seed(13342)
  fit2 <- fastgllvm(fg$Y, X=fg$X, q = q, family=family, hist=100, method="full", batch_size=100, trim=.3, intercept=intercept, alpha=.3, hessian=T, maxit=100, use_signs = F, H=1, rescale=T)
  plot(fit2)
  fit2 <- update(fit2, H=10, alpha=fit2$controls$alpha*10, maxit=10)
  plot(fit2)

  # simple, without rescaling
  set.seed(13342)
  fit3 <- fastgllvm(fg$Y, q = q, family=family, hist=100, method="simple", batch_size=500, trim=.2, intercept=intercept, alpha=.1, hessian=T, maxit=100, use_signs = F, H=1, rescale=T)
  plot(fit3)
  fit3 <- update(fit3, H=10, maxit=10)
  plot(fit3)
  # for(i in 1:5) {
  #   fit3 <- update(fit3, H=i, alpha=fit3$controls$alpha/2, trim=fit3$controls$trim/2, H=i, maxit=10)
  #   plot(fit3)
  # }

  # Approx for large p
  set.seed(1334)
  fit4 <- fastgllvm(fg$Y, q = q, family=family, hist=100, method="approx", batch_size=1000, trim=.1, intercept=intercept, alpha=.1, hessian=F, maxit=100, use_signs = T, H=1, rescale=F)
  fit4 <- update(fit4, H=10, alpha=fit4$controls$alpha/10, maxit=10)
  plot(fit4)
  # for(i in 1:5) {
  #   fit4 <- update(fit4, trim=fit4$controls$trim/1.5, H=i, maxit=10)
  #   plot(fit4)
  # }

  # fit <- update(fit, H=100, H.seed=1231, alpha=.1, maxit=10)
  # plot(fit)
  library(gllvm)
  fit.gllvm <- gllvm(y = fg$Y, X= fg$X, num.lv = q, formula = ~0 + ., family=binomial(), method = "EVA", sd.errors=F)

  library(mirtjml)
  fit.mirtjml <- mirtjml::mirtjml_expr(fg$Y, K=q, tol = .01)

  library(gmf)
  fit.gmf <- gmf(fg$Y,X = fg$X, family=poisson(), p=q, intercept = F)

  library(ltm)
  if(q==1) fit.ltm <- ltm(fg$Y ~ z1)
  if(q==2) fit.ltm <- ltm(fg$Y ~ z1 + z2)


  # plot(fg$parameters$A, psych::Procrustes(fit1$parameters$A, fg$parameters$A)$loadings, col=1)
  # abline(0,1,col=2)
  # points(fg$parameters$A, psych::Procrustes(fit2$parameters$A, fg$parameters$A)$loadings, col=2)
  # points(fg$parameters$A, psych::Procrustes(fit3$parameters$A, fg$parameters$A)$loadings, col=3)
  # points(fg$parameters$B, fit1$parameters$B, pch=2, col=1)

  plot(fg$parameters$A, psych::Procrustes(fit1$parameters$A, fg$parameters$A)$loadings, col=1)
  points(fg$parameters$A, psych::Procrustes(fit2$parameters$A, fg$parameters$A)$loadings, col=2)
  # points(fg$parameters$B, fit2$parameters$B, col=3)
  abline(0,1,col=2)
  points(fg$parameters$A, psych::Procrustes(fit3$parameters$A, fg$parameters$A)$loadings, col=3)
  points(fg$parameters$A, psych::Procrustes(fit4$parameters$A, fg$parameters$A)$loadings, col=4)
  points(fg$parameters$A, psych::Procrustes(fit.gmf$v, fg$parameters$A)$loadings, col=5)
  points(fg$parameters$A, psych::Procrustes(fit.ltm$coefficients[,2:(1+q)], fg$parameters$A)$loadings, col=6)#, xlim=c(-5,5), ylim=c(-5,5))
  # points(fg$parameters$B, fit$parameters$B, pch=2, col=1)
  points(fg$parameters$A, psych::Procrustes(fit.mirtjml$A_hat, fg$parameters$A)$loadings, col=3)
  # points(fg$parameters$B, t(fit.gmf$beta), col=4)
  # compute_error(fit0$parameters$A, fg$parameters$A)
  compute_error(fit1$parameters$A, fg$parameters$A)
  compute_error(fit2$parameters$A, fg$parameters$A)
  compute_error(fit3$parameters$A, fg$parameters$A)
  compute_error(fit4$parameters$A, fg$parameters$A)
  compute_error(fit.gllvm$params$theta, fg$parameters$A)
  compute_error(fit.gmf$v, fg$parameters$A)
  compute_error(fit.mirtjml$A_hat, fg$parameters$A)
  compute_error(fit.ltm$coefficients[,2:(1+q), drop=F], fg$parameters$A, rotate=T)

  plot(fg$Z, fit$Z); abline(0,-1,col=2); abline(0,1,col=2)
  points(fg$Z, -fit.gmf$u, col=2)

  plot(fg$parameters$A, psych::Procrustes(fit$parameters$A, fg$parameters$A)$loadings)
  points(fg$parameters$B, fit$parameters$B, pch=2)
  abline(0,1,col=2)
  points(fg$parameters$A, psych::Procrustes(fit.ltm$coefficients[,2:(1+q)], fg$parameters$A)$loadings, col=6)#, xlim=c(-5,5), ylim=c(-5,5))
  points(fg$parameters$B, fit.ltm$coefficients[,1], pch=2, col=2)
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
