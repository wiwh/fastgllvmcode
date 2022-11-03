#' Fit a fastgllvm object for fitting a gllvm model.
#'
fastgllvm.fit <- function(fg, method, parameters.init = NULL, controls=NULL, verbose=F, hist=F, median=F) {

  if(is.null(fg$parameters)) {
    # compute initial values for A, B phi, Z
    initial_values <- initialize_parameters(fg)
    fg$Z <- initial_values$Z
    fg$parameters <- list(A = initial_values$A,
                          B = initial_values$B,
                          phi = initial_values$phi)
    rm("initial_values")
    fg$parameters <- initialize_additional_parameters(fg, method)
    # fg$hessian <- initialize_hessian_AB(fg)
  } else {
    cat("\nInitial parameters values set to those provided.")
  }

  fg$mean <- fg$linpar <- matrix(0, nrow=fg$dimensions$n, ncol=fg$dimensions$p)

  gradients <- initialize_gradients(fg$parameters)


  # Beginning of the iterations
  # ---------------------------

  # impute for Y
  if(!is.null(fg$Miss)) {
    fg$Y[fg$Miss] <- fg$mean[fg$Miss]
  }

  # compute criterion after each pass...
  # only update the hessian once per pass...
  # TODO: compute_Z must accept linpar as argument! for starting values or whatever

  # Compute the Hessian
  # fg$hessian <- update_hessian_AB(fg$hessian, compute_hessian_AB(fg$X, fg$dimensions, fg$parameters, fg$families), weight=.9)
  # Importantly, the hessian should be computed independently from the rest
  fg$hessian <- compute_hessian_AB(fg$X, fg$dimensions, fg$parameters, fg$families)
  browser()

  # re-draw random batches every pass
  batches <- initialize_batches(fg$dimensions$n, controls$batch_size)
  for (batch in batches) {
    compute_dAB_centered(fg, method="simple")
    H_old <- hessian_old(fg)
    set.seed(123)
    compute_gradients(subset(fg, batch))
    # by the (inverse of ) Hessian

  }



  zstar <- compute_Z(Y, X, fg$parameters, fg$families, start=Z)
  params_hist <- list()
  if(hist) params_hist <- c(params_hist, list(parameters))

  moving_average <- parameters
  crit <- Inf

  if(is.null(controls$learning_rate.args)) controls$learning_rate.args <- list(method="spall", rate=2, end=.1)
  if(is.null(controls$learning_rate.args$method)) controls$learning_rate.args$method <- "spall"
  if(is.null(controls[["learning_rate"]])) controls$learning_rate <- initialize_learning_rate(maxit=controls$maxit, learning_rate.args = controls$learning_rate.args)

  for(i in 1:controls$maxit){
    # if(i < controls$maxit/2) median <- .5 else median <- F
    moving_average_old <- moving_average
    gradients_old <- gradients
    # get the gradient

    gradients <- compute_gradients(fg$Y, fg$X, parameters, fg$families, fg$Miss, debiase=T)

    # exponential smoothing step
    # for (par in names(parameters)) {
    #   gradients[[par]] <- (controls$beta * gradients[[par]] + (1-controls$beta) * gradients[[par]])
    # }
    parameters <- update_parameters(parameters, gradients, alpha=controls$alpha, learning_rate_i = controls$learning_rate(i), median=median)

    # We track a moving average of the last 1/controls$ma iterates
    for(k in seq_along(parameters)){
      moving_average[[k]] <- controls$ma * moving_average[[k]] + (1 - controls$ma) * parameters[[k]]
    }

    if(i == 1 || i%%10 == 0 || i == controls$maxit){
      crit <- compute_error(moving_average$A, moving_average_old$A, rotate=F)
      cat("\n Iteration: ", i, " - crit: ", crit)
    }
    if(hist) params_hist <- c(params_hist, list(parameters))
    if(i >= controls$minit && crit < controls$eps) break()
  }

  if(hist){
    history <- sapply(names(parameters), function(par_name) {
      do.call(rbind, lapply(params_hist, function(parameters_i) as.vector(parameters_i[[par_name]])))
    }, simplify=F)
  }

  # Update the fastgllvm object
  fg$parameters <- moving_average
  fg$fit <- list(
    crit = crit,
    controls = controls,
    hist = if(hist) history else NULL
  )
  fg$converged <- ifelse(crit < controls$eps, T, F)
  fg
}


if(0) {
  devtools::load_all()
  set.seed(1234)
  poisson  <- 10
  gaussian <- 10
  binomial <- 10
  nobs <- 100
  q <- 2
  p <- poisson + gaussian + binomial

  intercept <- T
  k <- 1
  if(k==0 & intercept) k <- 1
  family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial))
  set.seed(10030)
  fg <- gen_fastgllvm(nobs=nobs, p=p, q=q, k=k, family=family, intercept=intercept, phi=runif(p) + 0.5, miss.prob = 0, scale=1)

  set.seed(1304)
  fit.simple <- fastgllvm(fg$Y, X= fg$X, q = q, family=family,  intercept = T, hist=T, controls = list(maxit=200, alpha=5, beta=0, eps=1e-10, learning_rate.args=list(end=0.01, method="spall", rate=2)), method="simple", median=.5, batch_size=12)
  # set.seed(1304)
}
