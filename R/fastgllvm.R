#' Fit a fastgllvm
#' @param Y: a `n` times `p` matrix of observations
#' @param q: the number of factors
#' @param family: specifies the family. See below for more information.
#' @param X: either `NULL` or a `n` timeis `k` matrix of covariates.
#' @param intercept: a boolean (default:TRUE) indicating whether an intercept should be included in the model. If `X` is supplied and an intercept is desired, the first column of `X` must be a vector of ones to model the intercept.
#' @description
#'
#' @details
#' The implemented families currently are "gaussian", "binomial", or "poisson.
#' The `family` argument can be one of the following:
#'
#' * a string, the name of one of the implemented families.
#' * a vector of strings of size `p`, each of which being the name of the implemented families.
#' * a list: `list(objects, id, vec)`, where `objects` is a named list of objects of type "family",
#' `id` is a named list of integer vectors specifying the column numbers corresponding, whose names correspond to those of `objects`, and
#' `vec` is a named list of string vectors whose elements are the names corresponding to the `objects`$family name. This can be useful to specify
#' different link functions for the families.
#'
#' @examples
#' # Setting: declare the types of the response (in the order in which they appear)
#' family=c(rep("poisson", 50), rep("gaussian", 50), rep("binomial", 50))
#' q <- 2
#' p <- length(family)
#' set.seed(1234)
#' # Simulate data
#' data <- gen_fastgllvm(nobs=500, p=p, q=q, family=family, k=1, intercept=1, miss.prob = 0)
#' # Fit the data
#' fit <- fastgllvm(data$Y, q = q, X=data$X, family=family,  intercept = 1, controls = list(minit=20, maxit=100,alpha=1, eps=1e-3))
#' # Evaluate the fit
#' plot(fit)
#' If necessary, update the fit.
#' fit <- update(fit)
#' # Evaluate the fit
#' plot(fit)
#' @export

fastgllvm <- function(Y,
                      q=1,
                      family="gaussian",
                      X=NULL,
                      intercept=T,
                      Z.init = NULL,
                      parameters.init=NULL,
                      controls=list(),
                      verbose = F,
                      hist = T) {

  stopifnot(is.matrix(Y))

  families <- generate_families(family, ncol(Y))

  if (is.null(X)) {
    if (intercept) {
      X <- matrix(1, nrow(Y), 1)
      k <- 1
    } else {
      X <- matrix(0, nrow(Y), 0)
      k <- 0
    }
  } else {
    if (intercept && any(X[,1] != 1)) warning("When X is supplied, the intercept argument is ignored. Add a column of ones as the first column if you want to have an intercept.")
    k <- ncol(X)
  }

  dimensions <- list(
    p = ncol(Y),
    q = q,
    k = k
  )

  if (!is.null(Z.init)) {
    stopifnot(is.matrix(Z.init))
    Z <- Z.init
  } else {
    Z <- matrix(0, dimensions$p, dimensions$q)
  }

  # TODO: this is useless at this stage...
  if (!is.null(parameters.init)) {
    if(!is.list(parameters.init)) stop("The supplied 'parameters.init' must be a list.")
    if (any(!(names(parameters.init) %in% c("A", "B", "phi", "covZ")))) stop("Unrecognized name of parameter in supplied 'parameters.init'.")
    parameters <- parameters.init
  } else {
    parameters <- initialize_parameters(parameters.init$A, parameters.init$B, parameters.init$phi, dimensions$p, dimensions$q, dimensions$k)
  }


  stopifnot(is.list(controls))
  controls <- initialize_controls(controls)

  fg <- new_fastgllvm(Y, X, Z, parameters, families, dimensions)
  fastgllvm.fit(fg, controls, verbose, hist)

}


#' Workhorse function to fit a fastgllvm model
#' @param fg: an object of type "fastgllvm"
#' @param controls: a list of parameters controlling the fitting process
fastgllvm.fit <- function(fg, controls, verbose, hist, parameters.init=NULL) {
  #TODO: problem with missing value when q=1
  values <- list()
  if (!is.null(parameters.init)) {
    values$parameters <- parameters.init
  } else {
    values$parameters <- compute_parameters_initial_values(fg, rescale=F)
  }
  values$gradients = initialize_gradients(values$parameters)

  # One may define other functions to compute gradients
  compute_gradients_function <- compute_gradients

  params_hist <- list()
  moving_average <- values$parameters
  crit <- Inf

  learning_rate <- initialize_learning_rate(method="spall", maxit=controls$maxit, learning_rate.args = list(rate=3, end=.01))

  for(i in 1:controls$maxit){
    moving_average_old <- moving_average
    if (i < 10) {
      # safe start
      values <- update_parameters(fg$Y, fg$X, values$parameters, values$gradients, fg$families, fg$Miss, controls$alpha * learning_rate(i), controls$beta, debiase=ifelse(controls$safe.start, F, T), compute_gradients = compute_gradients_function)
    } else {
      # now we converge
      values <- update_parameters(fg$Y, fg$X, values$parameters, values$gradients, fg$families, fg$Miss, controls$alpha * learning_rate(i), controls$beta, debiase=T,  compute_gradients = compute_gradients_function)
    }
    for(k in seq_along(values$parameters)){
      moving_average[[k]] <- controls$ma * moving_average[[k]] + (1 - controls$ma) * values$parameters[[k]]
    }
    if(i%%5 == 1 || i == controls$maxit){
      crit <- compute_error(moving_average$A, moving_average_old$A, rotate=F)
      cat("\n Iteration: ", i, " - crit: ", crit)
    }
    if(hist) params_hist <- c(params_hist, list(values$parameters))
    if(i >= controls$minit && crit < controls$eps) break()
  }
  if(hist){
    history <- list()
    history$A <- do.call(rbind, lapply(params_hist, function(parameters_i) as.vector(parameters_i$A)))
    history$B <- do.call(rbind, lapply(params_hist, function(parameters_i) as.vector(parameters_i$B)))
    history$phi <- do.call(rbind, lapply(params_hist, function(parameters_i) as.vector(parameters_i$phi)))
    history$Z.cov <- do.call(rbind, lapply(params_hist, function(parameters_i) as.vector(parameters_i$covZ)))
  }

  # Update the fastgllvm object
  fg$Z <- moving_average[["Z"]]
  fg$parameters <- moving_average[c("A", "B", "phi", "covZ")]
  fg$fit <- list(
    crit = crit,
    controls = controls,
    hist = if(hist) history else NULL,
    learning_rate = learning_rate(1:i)
  )
  fg$converged <- ifelse(crit < controls$eps, T, F)
  fg
}


# Constructor
# -----------

#' Generates a fastgllvm object
#'
#' @param A the matrix of loadings.
#' @param B the matrix of fixed effect coefficients, of dimensions p * k
#' @param phi a vector of scale parameters.
#' @param X either 0 (no covariates, no intercept), 1 (an intercept), or a matrix of n * k covariates (with, possibly, the first column of 1s being an intercept)
#'
#' @return a list corresponding to the model
new_fastgllvm <- function(Y, X, Z, parameters, families, dimensions, linpar=NULL, fit=list(), Miss=NULL) {
  stopifnot(is.matrix(Y))
  if(!is.null(parameters$B)) {
    stopifnot(is.matrix(parameters$B))
    stopifnot(is.matrix(X))
    if(ncol(X) != ncol(parameters$B)) stop("Dimension mismatch between X and B.")
  }

  if (is.null(Z)) {
    Z <- matrix(0, nrow(Y), dimensons$q)
  }

  if (is.null(parameters)) {
    parameters <- list()
    parameters$A <- matrix(0, dimensions$p, dimensions$q)
    if (dimensions$k > 0) {
      parameters$B <- matrix(0, dimensions$p, dimensions$k)
    } else {
      parameters$B <- NULL
    }
    parameters$phi <- rep(1, dimensions$p)
  }


  stopifnot(is.matrix(parameters$A))
  stopifnot(is.vector(parameters$phi))
  stopifnot(is.list(fit))

  if (dimensions$q >= dimensions$p) {
    stop("The number of latent variables (q) must be strictly smaller than the number of observed variables (p).")
  }
  if (!is.vector(parameters$phi) || length(parameters$phi) != dimensions$p) {
    stop("phi must be a vector of length p.")
  }

  # check that the families have been correctly specified
  if(!all(sort(as.vector(do.call(c, families$id))) == 1:dimensions$p)) {
    stop("Family incorrectly specified: check the indices.")
  }

  if (is.null(Miss)) {
    Miss <- is.na(Y)
    if (!any(Miss)) {
      Miss <- NULL
    }
  }

  if (!is.null(Miss)) {
    if(any(rowcheck <- rowSums(!Miss) < (dimensions$q+1))) stop(paste0("Rows ", paste0(which(rowcheck), collapse = ","), " do not have enough observations."))
    if(any(colcheck <- colSums(!Miss) < (dimensions$q+1))) stop(paste0("Columns ", paste0(which(colcheck), collapse = ","), "do not have enough observations."))
  }

  fastgllvm <- structure(
    list(Y=Y,
         X=X,
         Z=Z,
         parameters=parameters,
         families=families,
         dimensions=dimensions,
         fit=fit,
         Miss=Miss),
    class="fastgllvm")
  fastgllvm
}


# validate a fastgllvm object
# TODO: families$id must total n!!!
validate_fastgllvm <- function(fastgllvm){

  # TODO: check that controls have the required stuff... or do that at control
  # TODO: test that the Miss matrix has no row with all ones...
  with(fastgllvm, {
  })
}

#' Generates a GLLVM model.
#'
#' Returns an (unfitted) gllvm model of class "fastglllvm" with simulated data.
#' Because it is a fastgllvm object, it has components such as "hist" that are NULL
#' since no fit actually took place.
#'
#' This function is useful for experimentation.
#'
#' @param nobs the number of observations to draw
#' @param p the number of manifest variables
#' @param q the number of latent variables
#' @param k the number of covariates common to all responses, k=0 yields an intercept of 0, while k=1 yields a random intercept for each response, and k>1 generates covariates.
#' @param family one of "normal", "bernoulli", or a vector of both of size p, specifying individual responses' distribution
#' @param par optional parameters in a list
#' @param sigma sigma
#'
#' @return a list corresponding to the model
#' @export
gen_fastgllvm <- function(nobs=100,
                          p=5,
                          q=1,
                          k=1,
                          family="gaussian",
                          A=NULL,
                          B=NULL,
                          phi=NULL,
                          X=NULL,
                          Z=NULL,
                          intercept=T,
                          scale=1,
                          miss.prob=0){

  # Generate unspecified parameters
  parameters <- generate_parameters(A, B, phi, p, q, k)
  if(scale!=1) {
    parameters$A <- parameters$A * scale
    parameters$B <- parameters$B * scale
  }
  dimensions <- generate_dimensions(parameters)


  # Generate unspecified variates
  if (is.null(Z)) Z <- gen_Z(nobs, q)
  if (is.null(X)) {
    if (k==0) {
      X <- matrix(NA, nobs, 0)
    } else {
      X <- gen_X(nobs, dimensions$k, intercept)
    }
  } else {
    if (intercept && any(X[,1] != 1)) warning("When X is supplied, the intercept argument is ignored. Add a column of ones as the first column if you want to have an intercept.")
  }
  if (intercept && k==0) warning("Intercept could not be added because k is set to 0. Set k>=1 if intercept=True.")

  # Construct families
  families <- generate_families(family, dimensions$p)

  # Generate data
  linpar <- compute_linpar(Z, parameters$A, X, parameters$B)
  Y <- generate_y(linpar, parameters$phi, families)$Y

  if (miss.prob != 0) {
    Y[runif(prod(dim(Y))) < miss.prob] <- NA
    Miss <- is.na(Y)
  } else {
    Miss <- NULL
  }

  fastgllvm <- new_fastgllvm(
    Y=Y,
    X=X,
    Z=Z,
    parameters=parameters,
    families=families,
    dimensions=dimensions,
    linpar=linpar,
    Miss = Miss
  )

  validate_fastgllvm(fastgllvm)
  fastgllvm
}


generate_parameters <- function(A, B, phi, p, q, k){
  if (is.null(A)) {
    if (is.null(p) | is.null(q)) stop("Either A, or p and q must be supplied.")
    A <- matrix(runif(p*q,-2, 2), p, q)
  } else {
    stopifnot(is.matrix(A))
    p <- nrow(A)
    q <- ncol(A)
  }
  if (is.null(B)) {
    if (is.null(p) | is.null(k)) stop("Either B, or p and k must be supplied")
    B <- matrix(runif(p*k, -1, 1), p, k)  # this can be a matrix with 0 columns.
  } else {
    stopifnot(is.matrix(B))
    stopifnot(nrow(B) == p)
    k <- ncol(B)
  }
  # if B has 0 col, then set it to NULL
  if(ncol(B)==0) B <- NULL


  if (is.null(phi)) {
    if(is.null(p)) stop("Either phi, or p must be supplied")
    phi <- rep(1, p)
  } else {
    stopifnot(length(phi) == p)
  }
  list(A=A, B=B, phi=phi)
}

generate_dimensions <- function(parameters){
  if (is.null(parameters$B)) {
    k <- 0
  } else {
    k <- ncol(parameters$B)
  }

  list(p=nrow(parameters$A), q=ncol(parameters$A), k=k)
}

#' Generate the families list as used in the other functions
generate_families <- function(family, p){
  families <- list(
    objects = list(),
    id = list(),
    vec = vector(length=p)
  )
  if (!is.list(family)) {
    if(!is.vector(family)) stop("Supplied 'family' must be either a list or a vector.")
    if(length(family) > 1) {
      fam_unique <- unique(family)
      for(i in 1:length(fam_unique)) {
        family_obj <- get(fam_unique[i], mode = "function", envir = parent.frame())()
        families$objects[[family_obj$family]] <- family_obj
        families$id[[family_obj$family]] <- which(family==fam_unique[i])
        if (is.null(family_obj$family)) {
          stop("'family' not recognized")
        }
      }
    } else {
      family_obj <- get(family, mode = "function", envir = parent.frame())()
      families$objects[[family_obj$family]] <- family_obj
      families$id[[family_obj$family]] <- 1:p
      if (is.null(family_obj$family)) {
        stop("'family' not recognized")
      }
    }
  } else {
    # TODO: test that the user entered the family correctly
    # families$id <- family
    # families$objects <- lapply(names(families$id), function(family_name)
    #   get(family_name, mode = "function", envir = parent.frame())())
    # names(families$objects) <- names(families$id)
  }

  for (i in seq_along(families$id)) {
    families$vec[families$id[[i]]] <- rep(families$objects[[i]]$family, length(families$id[[i]]))
  }
  families
}

# SOME TESTS
if(0) {
    devtools::load_all()
    set.seed(121234)
    poisson  <- 50
    gaussian <- 50
    binomial <- 50
    q <- 2
    p <- poisson + gaussian + binomial
    family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial))
    set.seed(1203)
    fg <- gen_fastgllvm(nobs=500, p=p, q=q, family=family, k=1, intercept=1, miss.prob = 0, scale=1)
    fit <- fastgllvm(fg$Y, q = q, X=fg$X, family=family,  intercept = 1, hist=T, controls = list(minit=20, maxit=100,alpha=1, eps=1e-3))
    plot(fit)
    fit <- update(fit)
    plot(fit)


  library(mirtjml)
  fit.m <- mirtjml_expr(fg$Y, q, tol = 1e-2)
  compute_error(fit.m$A_hat, fg$parameters$A, rotate = T)

  plot(fg$Z, fit.fg$Z); abline(0,-1,col=2)
  plot(fg$par$A, psych::Procrustes(fit.fg$parameters$A, fg$parameters$A)$loadings); abline(0,1,col=2)
  points(fg$par$A, psych::Procrustes(fit.m$A_hat, fg$parameters$A)$loadings, col=2)
  plot(fg$par$B, fit.fg$parameters$B, ylim=range(fg$parameters$B*1.5)); abline(0,1,col=2)
  points(fg$par$B, fit.m$d_hat, col=2)

  ts.plot(fit.fg$fit$hist$A[,1:min(100, p*q)])
  ts.plot(fit.fg$fit$hist$B[,1:p])
  ts.plot(fit.fg$fit$hist$Z.cov[,1:q])
  plot(fit.fg$fit$learning_rate)

  compute_error(fit.fg$parameters$B, fg$parameters$B)
  compute_error(fit.m$d_hat, fg$parameters$B)

  # load a simulated dataset
  attach(data_sim)
  # run the exploratory analysis
  res <- mirtjml_expr(data_sim$response, data_sim$K, tol = 1e-1)
  fit.fg <- fastgllvm(data_sim$response, data_sim$K, family="binomial", controls=list(eps=1e-4))

  plot(data_sim$A, psych::Procrustes(fit.fg$parameters$A, data_sim$A)$loadings, xlim=c(0,1.5), ylim=c(-.1, 2), pch=19); abline(h=0, col=2); abline(0,1,col=2)
  points(data_sim$A, psych::Procrustes(res$A_hat, data_sim$A)$loadings, col=2, pch=1)
  compute_error(fit.fg$parameters$A, data_sim$A, rotate = T)
  compute_error(res$A_hat, data_sim$A, rotate = T)
  plot(data_sim$d, fit.fg$parameters$B); abline(0,1,col=2)
  points(data_sim$d, res$d_hat, col=2)
  compute_error(fit.fg$parameters$B, matrix(data_sim$d))
  compute_error(matrix(res$d_hat), matrix(data_sim$d))
}
if(0){
  fastgllvm.fit <- function(Y, X, A.init=NULL, B.init=NULL, phi.init=NULL, Z.init=NULL, H=1, maxit=250 , tol=1e-5, learning_rate = NULL,  learning_rate.args = NULL, verbose = T ){
    if(is.null(learning_rate)) learning_rate <- ifelse(method=="SA", "exp", "constant")
    if(is.character(learning_rate)){
      learning_rate <- ff_learning_rate(method=learning_rate, maxit=maxit, learning_rate.args = learning_rate.args)
    }
    learning_rate.seq <- learning_rate(1:maxit)
    hist <- list(
      A = list(),
      B = list(),
      phi = list(),
      crit = list()
    )

    generate_Z <- generate_Z_functionfactory(n, q, method=method, H=H)

    A <- A.init
    B <- B.init
    phi <- phi.init
    Z <- generate_Z()[[1]]

    sol <- compute_pi(Y=Y, Z=Z, X=X, A=A, B=B, phi=phi, families=families, maxit=100)
    for(i in 1:100){
    sol <- compute_pi(Y=Y, Z=Z, X=X, A=A, B=B, phi=phi, families=families, maxit=1)
    sol$A <- psych::Procrustes(sol$A, f$A)$loadings
    plot(f$A, sol$A)


    A <- sol$A
    B <- sol$B
    phi <- sol$phi

    Z <- compute_zstar(Y, A, phi, X %*% t(B), families, start=Z, dims, scale=F, maxit=1)$Zstar
    print(diag(var(Z)))
    }
    plot(f$Z, Z)

    i <- 0
                            converged <- FALSE
                          while(i < maxit & !converged){
                            i <- i+1
                            Psi <- get_Psi(Y, Y.c, A, B, phi, X, family, generate_Z)

                            # udate A
                            A <- A + learning_rate.seq[i] * Psi$A
                            broken.A <- abs(A)>10
                            A[broken.A] <- 10 * sign(A[broken.A])

                            # update B
                            B  <- B + learning_rate.seq[i] * Psi$B

                            # update phi
                            phi <- phi + learning_rate.seq[i] * Psi$phi

                            # save
                            hist$A[hist.i + i, ] <- as.vector(A)
                            hist$B[hist.i + i, ] <- as.vector(B)
                            hist$phi[hist.i + i, ] <- as.vector(phi)

                            hist$crit[hist.i + i] <- learning_rate.seq[i] * norm(Psi$A)/(p*q)
                            if(hist$crit[hist.i + i] < tol) converged <- TRUE

                            if(verbose)cat("\ni: ", i, " - norm:", hist$crit[hist.i + i], " - learning rate:", learning_rate.seq[i])
                            # check if the criterion is small enough to jump to the next "repetition", where the learning rate increases again
                            if(converged){
                              # fill in the histories
                              hist$A <- hist$A[1:(hist.i + i),]
                              hist$B <- hist$B[1:(hist.i + i),]
                              hist$phi <- hist$phi[1:(hist.i + i),]
                              hist$crit <- hist$crit[1:(hist.i + i)]
                            }
                          }
    class(fastgllvm) <- "fastgllvm"
    fastgllvm
  }
}

