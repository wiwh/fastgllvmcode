#' Title
#'
#' @param fg: a gllvmprime model
#' @param parameters.init: a list of initial values for the parameters
#' @param controls: a list of controls for the fitting process
#'
#' @return a gllvmprime object fitted to the data
#' @export
#'
gllvmprime.fit <- function(fg, parameters.init = NULL, controls) {
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

  if (controls$method == "full") {
    compute_gradient <- gradient_full
  } else if (controls$method == "simple") {
    compute_gradient <- gradient_simple
  } else {
    stop("Unknown gradient method 'method'.")
  }

  hessian <- NULL

  fg <- compute_mean(fg, return_object = T)
  fg$deviance <- mean(compute_deviance(fg))

  params_hist <- list()
  if(!is.null(controls$hist)) params_hist <- c(params_hist, list(c(fg$parameters, deviance=fg$deviance)))

  # Initialize the moving average
  ma <- fg$parameters

  # Impute for Y
  if(!is.null(fg$Miss)) {
    fg$Y[fg$Miss] <- fg$mean[fg$Miss]
  }

  # Beginning of the iterations
  # ---------------------------
  for (i in 1:controls$maxit) {

    if (controls$hessian) {
      hessian_old <- hessian
      hessian_new <- simulate_hessian_AB(fg)
      hessian <- update_hessian(hessian_old, hessian_new, weight_old=.9)
    }

    step_size = (controls$alpha + 5)/ (5+i)

    if(!is.null(controls$H.seed)) warning("The seed is beeing reset! This corresponds to the sample path method, not the SA method.")

    # re-draw random batches every pass
    batches <- initialize_batches(fg$dimensions$n, controls$batch_size)

    for (batch in batches) {
      fg_batch <- subset(fg, batch)

      if (!is.null(controls$H.seed)) {
        set.seed(controls$H.seed)
      }
      # Compute the simulated gradients
      if (controls$H > 1) {
        sims <- lapply(1:controls$H, function(na) {
          fg_simul <- simulate(fg_batch, return_object=TRUE)
          compute_gradient(fg_simul)
        })
        grad_simul <- list(
          AB = Reduce("+", lapply(sims, function(sim)sim$AB)) / length(sims),
          phi = Reduce("+", lapply(sims, function(sim)sim$phi)) / length(sims),
          covZ = Reduce("+", lapply(sims, function(sim)sim$covZ)) / length(sims)
        )
        if(is.null(fg$parametersB)) grad_simul$B <- NULL
      } else {
        fg_simul <- simulate(fg_batch, return_object=TRUE)

        grad_simul <- compute_gradient(fg_simul)
      }

      grad_batch <- compute_gradient(fg_batch)
      fg$parameters <- compute_parameters_update(fg_batch, step_size, hessian, grad_batch, grad_simul, controls)
    }

    # Update Z and the estimated means after the parameter update
    fg <- compute_Z(fg, start=fg$Z, return_object = T)
    fg <- compute_mean(fg, return_object = T)  # TODO: compute mean based on ma mayhaps?


    # Rescale the model: this helps the training, but does not affect the result
    if (controls$rescale) {
      fg_rescaled <- rescale(fg, rescale.A = T, rescale.B = T, target.cov = fg$parameters$covZ)
      fg$parameters$A <- fg$parameters$A * .9 + fg_rescaled$parameters$A * .1
      fg$parameters$B <- fg$parameters$B * .95 + fg_rescaled$parameters$B * .05
      # fg$parameters$A <- fg_rescaled$parameters$A
      # fg$parameters$B <- fg_rescaled$parameters$B
      fg$Z <- fg_rescaled$Z
    }

    # Impute values for Y
    if(!is.null(fg$Miss)) {
      fg$Y[fg$Miss] <- fg$mean[fg$Miss]
    }
    ma <- update_moving_average(ma, fg$parameters, 0.9)

    fg$deviance <- mean(compute_deviance(fg))

    cat("\nIteration: ", i, "dev:", fg$deviance, "range: ", range(fg$parameters$A))

    if(!is.null(controls$hist)){
      if (length(params_hist) > controls$hist) params_hist[[1]] <- NULL
      params_hist <- c(params_hist, list(c(fg$parameters, deviance=fg$deviance)))
    }
  }

  fg$parameters <- ma
  fg$Z <- compute_Z(fg, start=fg$Z)$Z
  fg <- compute_mean(fg, return_object = T)

  if(!is.null(controls$hist)){
    history <- sapply(names(params_hist[[1]]), function(par_name) {
      do.call(rbind, lapply(params_hist, function(parameters_i) as.vector(parameters_i[[par_name]])))
    }, simplify=F)
  }

  fg$hist <- if(controls$hist) history else NULL
  fg$controls <- controls
  fg$controls$alpha <- step_size
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

# update hessian
update_hessian <- function(hessian_old, hessian_new, weight_old) {
  if (is.null(hessian_old)) {
    return(hessian_new)
  }
  lapply(seq_along(hessian_old), function(h_i){
    weight_old * hessian_old[[h_i]]  + (1-weight_old) * hessian_new[[h_i]]
  })
}

if(0) {
  devtools::load_all()
  poisson  <- 0
  gaussian <- 0
  binomial <- 4
  nobs <- 200
  q <- 1

  p <- poisson + gaussian + binomial

  intercept <- T
  k <- 1
  if(k==0 & intercept) k <- 1


  family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial))

  set.seed(14240)
  fg <- gen_gllvmprime(nobs=nobs, p=p, q=q, k=k, family=family, intercept=intercept, miss.prob = 0, scale=1, phi=rep(1, p))


  # full, with hessian, no rescaling
  set.seed(13342)
  fit1 <- gllvmprime(fg$Y, q = q, family=family, maxit=100, hist=200, method="full", trim=.1, intercept=intercept, alpha=.5, hessian=T, rescale=T)
  plot(fit1)
  fit1 <- update(fit1, rescale=F)
  plot(fit1)

  # fit1 <- update(fit1, H=10, maxit=10)

  plot(fit1)
  MPE(fit1$parameters$A, fg$parameters$A)

  # clear winner in poisson
  # full, with rescaling outside the loopa-
  set.seed(13342)
  fit2 <- gllvmprime(fg$Y, X=fg$X, q = q, family=family, hist=100, method="full", trim=.5, intercept=intercept, alpha=.1, hessian=T, maxit=20, use_signs = F, H=1, rescale=T)
  fit2 <- update(fit2, H=1, alpha=.05, maxit=50)
  plot(fit2)
  fit2 <- update(fit2, H=10, alpha=fit2$controls$alpha*10, maxit=10)
  fit2 <- update(fit2, H=10, maxit=10)


  # simple, with rescaling
  set.seed(13342)
  fit3 <- gllvmprime(fg$Y, X=fg$X, q = q, family=family, hist=100, method="simple", batch_size=100, trim=.1, intercept=intercept, alpha=.1, hessian=T, maxit=100, use_signs = F, H=1, rescale=F)
  plot(fit3)
  fit3 <- update(fit2, H=10, alpha=fit2$controls$alpha*10, maxit=10)
  plot(fit3)

  # simple, with rescaling
  set.seed(13342)
  fit4 <- gllvmprime(fg$Y, X=fg$X, q = q, family=family, hist=100, method="simple", batch_size=500, trim=.3, intercept=intercept, alpha=.3, hessian=T, maxit=100, use_signs = F, H=1, rescale=F)
  plot(fit4)
  fit4 <- update(fit2, H=10, alpha=fit2$controls$alpha*10, maxit=10)
  plot(fit4)

  # Approx for large p
  set.seed(1334)
  fit5 <- gllvmprime(fg$Y, q = q, family=family, hist=100, method="approx", batch_size=1000, trim=.1, intercept=intercept, alpha=.1, hessian=F, maxit=100, use_signs = T, H=1, rescale=F)
  fit5 <- update(fit4, H=10, alpha=fit4$controls$alpha/10, maxit=10)
  plot(fit5)

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
  fit.gmf <- gmf(fg$Y,X = fg$X, family=binomial(), p=q, intercept = F, method="quasi")


  library(ltm)
  if(q==1) fit.ltm <- ltm(fg$Y ~ z1)
  if(q==2) fit.ltm <- ltm(fg$Y ~ z1 + z2)

  # plot(fg$parameters$A, psych::Procrustes(fit1$parameters$A, fg$parameters$A)$loadings, col=1)
  # abline(0,1,col=2)
  # points(fg$parameters$A, psych::Procrustes(fit2$parameters$A, fg$parameters$A)$loadings, col=2)
  # points(fg$parameters$A, psych::Procrustes(fit3$parameters$A, fg$parameters$A)$loadings, col=3)
  # points(fg$parameters$B, fit1$parameters$B, pch=2, col=1)

  plot(fg$parameters$A, psych::Procrustes(fit1$parameters$A, fg$parameters$A)$loadings, col=1)
  abline(0,1,col=2)
  points(fg$parameters$A, psych::Procrustes(fit2$parameters$A, fg$parameters$A)$loadings, col=2)
  # points(fg$parameters$B, fit2$parameters$B, col=3)
  points(fg$parameters$A, psych::Procrustes(fit3$parameters$A, fg$parameters$A)$loadings, col=3)
  points(fg$parameters$A, psych::Procrustes(fit4$parameters$A, fg$parameters$A)$loadings, col=4)
  points(fg$parameters$A, psych::Procrustes(fit5$parameters$A, fg$parameters$A)$loadings, col=5)

  points(fg$parameters$A, psych::Procrustes(fit.gmf$v, fg$parameters$A)$loadings, col=5)
  points(fg$parameters$A, psych::Procrustes(fit.ltm$coefficients[,2:(1+q)], fg$parameters$A)$loadings, col=7)#, xlim=c(-5,5), ylim=c(-5,5))
  # points(fg$parameters$B, fit$parameters$B, pch=2, col=1)
  points(fg$parameters$A, psych::Procrustes(fit.mirtjml$A_hat, fg$parameters$A)$loadings, col=3)
  # points(fg$parameters$B, t(fit.gmf$beta), col=4)
  # compute_error(fit0$parameters$A, fg$parameters$A)
  MPE(fit1$parameters$A, fg$parameters$A)
  MPE(fit2$parameters$A, fg$parameters$A)
  MPE(fit3$parameters$A, fg$parameters$A)
  MPE(fit4$parameters$A, fg$parameters$A)
  MPE(fit.gllvm$params$theta, fg$parameters$A)
  MPE(fit.gmf$v, fg$parameters$A)
  MPE(fit.mirtjml$A_hat, fg$parameters$A)
  MPE(fit.ltm$coefficients[,2:(1+q), drop=F], fg$parameters$A, rotate=T)

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
