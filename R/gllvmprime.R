#' Fit a GLLVM using an M-estimator a stochatic approximation algorithm.
#' @param Y: a `n` times `p` matrix, organized with the observations in the rows and the responses in the columns
#' @param q: the number of factors
#' @param family: specifies the families. See below for more information.
#' @param X: either `NULL` or a `n` timeis `k` matrix of covariates.
#' @param intercept: a boolean (default:TRUE) indicating whether an intercept should be included in the model. If `X` is supplied and an intercept is desired, the first column of `X` must be a vector of ones to model the intercept.
#' @param gradient_function: one of: "simple", "full", or a user-defined function.
#' @param parameters.init: a list containing "A", "B", and "phi".
#' @param hessian: one of "median", "F", or "T"
#' @param hist: the last `hist` iterations are kept for the graph. Do not keep any if set to NULL.
#' @description Fits a gllvm model.
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
#' data <- gen_gllvmprime(nobs=100, p=p, q=q, family=family, k=1, intercept=1, miss.prob = 0)
#' # Fit the data
#' fit <- gllvmprime(data$Y, q = q, X=data$X, family=family, alpha=.5, maxit=20)
#' # Evaluate the fit
#' plot(fit)
#' # If necessary, continue the fitting process:
#' fit <- update(fit)
#' # Evaluate the fit
#' plot(fit)
#' @export

gllvmprime <- function(Y,
                      q=1,
                      family="gaussian",
                      X=NULL,
                      intercept=T,
                      Z.init = NULL,
                      parameters.init=NULL,
                      method="full",
                      verbose = F,
                      hist = 100,
                      hessian = T,
                      use_signs = F,
                      trim = 0.2,
                      alpha = 0.2,
                      batch_size=nrow(Y),
                      maxit=100,
                      H=1,
                      H.seed = NULL,
                      rescale=T) {

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
    if (intercept && ncol(X) > 0 && any(X[,1] != 1)) warning("When X is supplied, the intercept argument is ignored. Add a column of ones as the first column if you want to have an intercept.")
    k <- ncol(X)
  }

  dimensions <- list(
    n = nrow(Y),
    p = ncol(Y),
    q = q,
    k = k
  )

  controls <- list(
    batch_size = batch_size,
    trim = trim,
    method = method,
    hessian = hessian,
    use_signs = use_signs,
    hist = hist,
    alpha = alpha,
    verbose = verbose,
    maxit = maxit,
    H = H,
    H.seed = H.seed,
    rescale = rescale
  )

  # get mising values
  Miss <- is.na(Y)
  if (!any(Miss)) {
    Miss <- NULL
  }

  fg <- new_gllvmprime(Y, X, Z.init, parameters.init, families, dimensions, Miss)

  gllvmprime.fit(fg, parameters.init=parameters.init, controls=controls)
}

# Constructor
# -----------

#' Generates a gllvmprime object
#'
#' @param A the matrix of loadings.
#' @param B the matrix of fixed effect coefficients, of dimensions p * k
#' @param phi a vector of scale parameters.
#' @param X either 0 (no covariates, no intercept), 1 (an intercept), or a matrix of n * k covariates (with, possibly, the first column of 1s being an intercept)
#'
#' @return a list corresponding to the model
new_gllvmprime <- function(Y, X, Z, parameters, families, dimensions, Miss, fit=list()) {
  gllvmprime <- structure(
    list(Y=Y,
         X=X,
         Z=Z,
         parameters=parameters,
         families=families,
         dimensions=dimensions,
         Miss=Miss,
         fit=fit),
    class="gllvmprime")

  validate_gllvmprime(gllvmprime)

  gllvmprime
}



#' Validates a gllvmprime object
#'
#' @param gllvmprime: a gllvmprime object
#'
#' Strives to find errors in the gllvmprime object to simplify debugging.
validate_gllvmprime <- function(gllvmprime) {
  with(gllvmprime, {
    stopifnot(is.matrix(Y))
    #' Check whether the parameters are correctly specified
    if(!is.null(parameters)) {
      stopifnot(is.matrix(parameters$A))
      stopifnot(all(dim(parameters$A) == c(dimensions$p, dimensions$q)))

      if (!is.vector(parameters$phi) || length(parameters$phi) != dimensions$p) {
        stop("phi must be a vector of length p.")
      }
    }

    if(!is.null(parameters$B)) {
      stopifnot(is.matrix(parameters$B))
      stopifnot(is.matrix(X))
      if(ncol(X) != ncol(parameters$B)) stop("Dimension mismatch between X and B.")
    }

    stopifnot(is.list(fit))


    if (dimensions$q >= dimensions$p) {
      stop("The number of latent variables (q) must be strictly smaller than the number of observed variables (p).")
    }

    # check that the families have been correctly specified
    if(!all(sort(as.vector(do.call(c, families$id))) == 1:dimensions$p)) {
      stop("Family incorrectly specified: check the indices.")
    }

    if (!is.null(Miss)) {
      if(any(rowcheck <- rowSums(!Miss) < (dimensions$q+1))) stop(paste0("Rows ", paste0(which(rowcheck), collapse = ","), " do not have enough observations."))
      if(any(colcheck <- colSums(!Miss) < (dimensions$q+1))) stop(paste0("Columns ", paste0(which(colcheck), collapse = ","), "do not have enough observations."))
    }
  })
}






#' Generates a GLLVM model.
#'
#' Returns an (unfitted) gllvm model of class "fastglllvm" with simulated data.
#' Because it is a gllvmprime object, it has components such as "hist" that are NULL
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
gen_gllvmprime <- function(nobs=100,
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

  dimensions <- list(n=nobs,
                     p=nrow(parameters$A),
                     q=ncol(parameters$A),
                     k=if(is.null(parameters$B)) 0 else ncol(parameters$B))



  if (is.null(X)) {
    if (dimensions$k==0) {
      X <- matrix(NA, dimensions$n, 0)
    } else {
      X <- gen_X(dimensions$n, dimensions$k, intercept)
    }
  } else {
    if (intercept && any(X[,1] != 1)) warning("When X is supplied, the intercept argument is ignored. Add a column of ones as the first column if you want to have an intercept.")
  }
  if (intercept && dimensions$k==0) warning("Intercept could not be added because k is set to 0. Set k>=1 if intercept=True.")

  # Construct families
  families <- generate_families(family, dimensions$p)

  # Generate data
  variables <- gen_Y(Z=Z, X=X, parameters = parameters, families = families)

  if (miss.prob != 0) {
    Y0 <- variables$Y
    variables$Y[runif(prod(dim(variables$Y))) < miss.prob] <- NA
    Miss <- is.na(variables$Y)
  } else {
    Miss <- NULL
  }

  gllvmprime <- new_gllvmprime(
    Y=variables$Y,
    Z=variables$Z,
    X=X,
    parameters=parameters,
    families=families,
    dimensions=dimensions,
    Miss = Miss
  )

  if(!is.null(Miss)) gllvmprime$Y0 <- Y0
  gllvmprime$intercept <- intercept
  gllvmprime$parameters_true <- gllvmprime$parameters

  validate_gllvmprime(gllvmprime)
  gllvmprime
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



#' Generate the families list as used in the other functions
#' @param family: either a string, or a vector of strings of length p specifying the name(s) of the families to be used: each element must be one of "gaussian", "binomial", or "poisson". For now, only the canonical link function is used so it needs not be specified.
generate_families <- function(family, p){
  stopifnot(is.character(family))
  if (!(length(family) %in% c(1,p))) stop("Length of the 'family' vector must be either 1 or p.")

  families <- list(
    objects = list(),
    id = list(),
    vec = vector(length=p)
  )

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

  for (i in seq_along(families$id)) {
    families$vec[families$id[[i]]] <- rep(families$objects[[i]]$family, length(families$id[[i]]))
  }

  families
}


# SOME TESTS
if(0) {
    devtools::load_all()
    set.seed(1234)
    poisson  <- 1000
    gaussian <- 1000
    binomial <- 1000
    nobs <- 100
    q <- 5
    p <- poisson + gaussian + binomial

    intercept <- T
    k <- 1
    if(k==0 & intercept) k <- 1
    family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial))
    set.seed(10030)
    fg <- gen_gllvmprime(nobs=nobs, p=p, q=q, k=k, family=family, intercept=intercept, phi=runif(p) + 0.5, miss.prob = 0, scale=1)

    plot(simulate(fg, conditional=T)$Y, fg$Y)
    # check initialization
    # set.seed(123)
    # fit <- gllvmprime(fg$Y, q = q, family=family,  intercept = T, hist=T, controls = list(maxit=100,alpha=.5, beta=0, eps=1e-10, learning_rate.args=list(end=0.01, method="constant")), median=F)
    # plot(fit)
    set.seed(1304)
    fit.simple <- gllvmprime(fg$Y, X= fg$X, q = q, family=family,  intercept = T, hist=T, controls = list(maxit=200, alpha=5, beta=0, eps=1e-10, learning_rate.args=list(end=0.01, method="spall", rate=2)), method="simple", median=.5, batch_size=12)
    # set.seed(1304)
    # fit.full <- gllvmprime(fg$Y, X= fg$X, q = q, family=family,  intercept = T, hist=T, controls = list(maxit=500, alpha=5, beta=0, eps=1e-10, learning_rate.args=list(end=0.01, method="spall", rate=10)), method="full", median=.2)

    plot(fit.simple)
    # plot(fit.full)
    plot(fit.simple$fit$controls$learning_rate(1:200))

    plot(fg$parameters$A, psych::Procrustes(fit.simple$parameters$A, fg$parameters$A)$loadings, col=3)
    points(fg$parameters$B, fit.simple$parameters$B, pch=2, col=3)
    # LTM TEST
    library(ltm)
    if(q==1) fit.ltm <- ltm(fg$Y ~ z1)
    if(q==2) fit.ltm <- ltm(fg$Y ~ z1 + z2)

    points(fg$parameters$A, psych::Procrustes(fit.ltm$coefficients[,2:(1+q)], fg$parameters$A)$loadings)#, xlim=c(-5,5), ylim=c(-5,5))
    points(fg$parameters$B, fit.ltm$coefficients[,1], pch=2)
    # points(fg$parameters$A, psych::Procrustes(fit.full$parameters$A, fg$parameters$A)$loadings, col=2)
    # points(fg$parameters$B, fit.full$parameters$B, pch=2, col=2)
    legend(x="bottomright",pch=1,  legend=c("ltm","full", "simple"), col=1:3)
    abline(0,1,col=2)
    abline(0,-1,col=2)
    compute_error(fit.ltm$coefficients[,2:(1+q), drop=F], fg$parameters$A, rotate=T)
    # compute_error(fit$parameters$A, fg$parameters$A, rotate=T)
    compute_error(fit.simple$parameters$A, fg$parameters$A, rotate=T)
    # compute_error(fit.full$parameters$A, fg$parameters$A, rotate=T)

    plot(fg$parameters$A, psych::Procrustes(fit.simple$parameters$A, fg$parameters$A)$loadings);abline(0,1,col=2)
    plot(fg$parameters$B, fit.simple$parameters$B);abline(0,1,col=2)


    plot(fg$Z, fit.simple$Z)

    # GAUSSIAN TEST
    fit.ffa <- ffa(fg$Y, q, iteratively_update_Psi = T)
    cat("\nError fg:", compute_error(fit.simple$parameters$A, fg$parameters$A, rotate=T))
    cat("\nError ffa:", compute_error(fit.ffa$A, fg$parameters$A, rotate=T))

    compute_cov <- function(A, Psi){
      cov <- A %*% t(A) + diag(Psi)
      cov <- diag(diag(cov)**-.5) %*% cov %*% diag(diag(cov)**-.5)
      cov
    }

    plot(fg$parameters$phi, fit.ffa$Psi);points(fg$parameters$phi, fit.simple$parameters$phi, col=2); abline(0,1,col=2)
    plot(cor(fg$Y), compute_cov(fit.ffa$A, fit.ffa$Psi))
    points(cor(fg$Y), compute_cov(fit.simple$parameters$A, fit.simple$parameters$phi), col=3); abline(0,1,col=2)

    fa <- factanal(fg$Y, q)
    fit <- fit.simple
    cov1 <- fit$parameters$A %*% t(fit$parameters$A) + diag(fit$parameters$phi)
    cov1 <- diag(diag(cov1)**-0.5)%*% cov1 %*%diag(diag(cov1)**-0.5)
    cov2 <- fa$loadings %*% t(fa$loadings) + diag(fa$uniquenesses)
    plot(cor(fg$Y), cov1)
    points(cor(fg$Y), cov2, col=2); abline(0,1)
    plot(fg$parameters$A, psych::Procrustes(fit$parameters$A, fg$parameters$A)$loadings);abline(0,1,col=2)
    points(fg$parameters$A, psych::Procrustes(fa$loadings, fg$parameters$A)$loadings, col=2);abline(0,1,col=2)




  # For binary
  library(mirtjml)
  fit.m <- mirtjml_expr(fg$Y, q, tol = 1e-2)

  fit.fg <- fit.simple

  compute_error(fit.m$A_hat, fg$parameters$A, rotate = T)
  compute_error(fit.fg$parameters$A, fg$parameters$A, rotate=T)


  plot(fg$Z, fit.fg$Z); abline(0,-1,col=2)
  plot(fg$par$A, psych::Procrustes(fit.fg$parameters$A, fg$parameters$A)$loadings); abline(0,1,col=2)
  points(fg$par$A, psych::Procrustes(fit.m$A_hat, fg$parameters$A)$loadings, col=2)
  plot(fg$par$B, fit.fg$parameters$B, ylim=range(fg$parameters$B*1.5)); abline(0,1,col=2)
  points(fg$par$B, fit.m$d_hat, col=2)

  ts.plot(fit.fg$fit$hist$A[,1:min(100, p*q)])
  ts.plot(fit.fg$fit$hist$B[,1:p])
  ts.plot(fit.fg$fit$hist$Z[,1:min(100, fg$dimensions$n*q)])

  compute_error(fit.fg$parameters$B, fg$parameters$B)
  compute_error(fit.m$d_hat, fg$parameters$B)

  # load a simulated dataset
  attach(data_sim)
  # run the exploratory analysis
  res <- mirtjml_expr(data_sim$response, data_sim$K, tol = 1e-1)
  fit.fg <- gllvmprime(data_sim$response, data_sim$K, family="binomial", controls=list(eps=1e-4))

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
  gllvmprime.fit <- function(Y, X, A.init=NULL, B.init=NULL, phi.init=NULL, Z.init=NULL, H=1, maxit=250 , tol=1e-5, learning_rate = NULL,  learning_rate.args = NULL, verbose = T ){
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
    class(gllvmprime) <- "gllvmprime"
    gllvmprime
  }
}

