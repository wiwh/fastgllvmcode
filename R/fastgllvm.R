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
new_fastgllvm <- function(Y, X, A, B, phi, families, fit=list(), linpar=list()) {
  stopifnot(is.matrix(Y))
  stopifnot(is.matrix(X))
  stopifnot(is.matrix(A))
  stopifnot(is.matrix(B))
  stopifnot(is.vector(phi))
  stopifnot(attr(family, "class") == "family")
  stopifnot(is.list(fit))

  n <- nrow(Y)
  p <- nrow(A)
  q <- ncol(A)
  k <- ncol(B)
  if (q >= p)
    stop("The number of latent variables (q) must be strictly smaller than the number of observed variables (p).")
  if (nrow(B) != p)
    stop("The number of rows of B must be equal to the number of observed variables (p).")
  if (!is.vector(phi) || length(phi) != p)
    stop("phi must be a vector of length p.")

  if (!is.matrix(X)) {
    if (X==0 & length(B) > 0)
      stop("No intercept given, yet fixed effect coefficients are given. Should there be an intercept?")
    if (X==1 & ncol(B) !=1)
      stop("The dimensions of the fixed effect coefficients matrix B is inconsistent with the value of X implying only an intercept is required.")
  } else {
    if (nrow(X) != n || ncol(X) != k)
      stop("X must be  either 0 or 1, or a matrix of dimensions p * k.")
  }

  fastgllvm <- structure(
    list(Y=Y,
         X=X,
         A=A,
         B=B,
         phi=phi,
         families=families,
         dims=list(n=n, p=p, q=q, k=k),
         fit=fit,
         linpar=linpar),
    class="fastgllvm")
  fastgllvm
}


# validate a fastgllvm object
# TODO: families$id must total n!!!
validate_fastgllvm <- function(fastgllvm){
  with(fastgllvm,{
    stopifnot(is.matrix(A))
    stopifnot(is.matrix(B))
    stopifnot(is.vector(phi))
  })

  # check that the families have been correctly specified
  with(fastgllvm, {
    if(!all(sort(as.vector(do.call(c, families$id))) == 1:dims$p)) {
      stop("Family incorrectly specified: check the indices.")
    }
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
                          k=0,
                          family="gaussian",
                          A=NULL,
                          B=NULL,
                          phi=NULL,
                          X=NULL,
                          Z=NULL,
                          intercept=F){

  # Generate unspecified parameters
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
  if (is.null(phi)) {
    if(is.null(p)) stop("Either phi, or p must be supplied")
    phi <- rep(1, p)
  } else {
    stopifnot(length(phi) == p)
  }

  # Generate unspecified variates
  if (is.null(Z)) Z <- gen_Z(nobs, q)
  if (is.null(X)) {
    if (k==0) {
      X <- matrix(NA, nobs, 0)
    } else {
      X <- gen_X(nobs, k, intercept)
    }
  } else {
    if (intercept && any(X[,1] != 1)) warning("When X is supplied, the intercept argument is ignored. Add a column of ones as the first column if you want to have an intercept.")
  }
  if (intercept && k==0) warning("Intercept could not be added because k is set to 0. Set k>=1 if intercept=True.")

  # Construct families
  families <- generate_families(family, p)

  # Generate data
  linpar <- compute_linpar(Z, A, X=X, B=B)
  Y <- generate_y(linpar, phi, families, nobs, p, q)

  fastgllvm <- new_fastgllvm(
    Y=Y,
    X=X,
    A=A,
    B=B,
    phi=phi,
    families=families,
    linpar=linpar
  )

  fastgllvm$Z <- Z

  validate_fastgllvm(fastgllvm)
  fastgllvm
}

#' Generate the families list as used in the other functions
generate_families <- function(family, p){
  # objects are required for later use
  families <- list(
    objects = list(),
    id = list(),
    vec = vector(length=p)
  )

  if (!is.list(family)) {
    if (is.character(family))
      family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
      family <- family()
    if (is.null(family$family)) {
      print(family)
      stop("'family' not recognized")
    }
    families$objects[[family$family]] <- family
    families$id[[family$family]] <- 1:p
  } else {
    # TODO: test that the user entered the family correctly
    families$id <- family
    families$objects <- lapply(names(families$id), function(family_name)
      get(family_name, mode = "function", envir = parent.frame())())
    names(families$objects) <- names(families$id)
  }
  for (i in seq_along(families$id)) {
    families$vec[families$id[[i]]] <- rep(families$objects[[i]]$family, length(families$id[[i]]))
  }
  families
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

