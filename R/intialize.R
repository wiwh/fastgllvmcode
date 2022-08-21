initialize_parameters <- function(parameters, dimensions, method) {
  if (method == "full") {
    parameters <- initialize_parameters_full(parameters, dimensions)
  }
  if (method == "simple") {
    parameters <- initialize_parameters_simple(parameters, dimensions)
  }
  parameters
}

initialize_controls <- function(controls) {
  stopifnot(is.list(controls))
  if (is.null(controls[["maxit"]])) controls$maxit <- 100
  if (is.null(controls[["alpha"]])) controls$alpha <- 1
  if (is.null(controls[["beta"]]))  controls$beta  <- 0
  if (is.null(controls[["ma"]]))    controls$ma    <- .9
  if (is.null(controls[["eps"]]))   controls$eps   <- 1e-3
  if (is.null(controls[["safe.start"]]))   controls$safe.start   <- F
  if (is.null(controls[["minit"]])) {
    controls$minit <- min(20, 2/(1 - controls$ma))
  } else {
    if (controls$minit < 2/(1 - controls$ma)) warnings(paste0("It is recommended to set a minit of ", min(20, 2/(1 - controls$ma)), "or larger"))
  }

  if(controls$minit > controls$maxit) stop("Set a maxit larger than minit.")

  controls
}

initialize_gradients <- function(parameters, method){
  if (method == "full") {
    gradients <- initialize_gradients_full(parameters)
  }
  if( method == "simple") {
    gradients <- initialize_gradients_simple(parameters)
  }
  gradients
}

#' Function factory to return a function that generates a learning rate as a function of i
#'
#' @param method: the name of the method to use. One of "constant", "lin", "exp", "linexp", "spall".
#' @param method.args: list of method arguments
initialize_learning_rate <- function(maxit, learning_rate.args=list(method="spall")){
  if(is.null(learning_rate.args))       learning_rate.args <- list()
  if(is.null(learning_rate.args$start)) learning_rate.args$start <- 1
  if(is.null(learning_rate.args$end))   learning_rate.args$end <- 0.005
  if(is.null(learning_rate.args$constant)) learning_rate.args$constant <- 1

  if(learning_rate.args$method=="constant"){
    learning_rate <- function(i){
      rep(learning_rate.args$constant, length(i))
    }
  }
  if(learning_rate.args$method=="lin"){
    lr <- with(learning_rate.args,
               seq(start, end, l=maxit))
    learning_rate <- function(i){
      lr[i]
    }
  }
  if(learning_rate.args$method=="exp"){
    lr <- with(learning_rate.args, {
      rho <- exp(log(end/start)/(maxit-1))
      rep(start, maxit) * rho**(0:(maxit-1))
    })
    learning_rate <- function(i){
      lr[i]
    }
  }
  if(learning_rate.args$method=="linexp"){
    if(is.null(learning_rate.args$weight.exp))learning_rate.args$weight.exp <- .5
    lr.lin <- with(learning_rate.args,
                   seq(start, end, l=maxit))
    lr.exp <- with(learning_rate.args, {
      rho <- exp(log(end/start)/(maxit-1))
      rep(start, maxit) * rho**(0:(maxit-1))
    })
    lr <- (1-learning_rate.args$weight.exp) * lr.lin + learning_rate.args$weight.exp * lr.exp
    learning_rate <- function(i){
      lr[i]
    }
  }
  if(learning_rate.args$method=="spall"){
    if(is.null(learning_rate.args$rate)) learning_rate.args$rate <- 1000
    lr <- with(learning_rate.args, {
      b <- exp(log(end / start) / rate)
      b <- (maxit * b - 1) / (1 - b)
      start*((1+b)/(1:maxit + b))**rate
    })
    learning_rate <- function(i){
      lr[i]
    }
  }
  learning_rate
}



#' Given a fastgllvm object, compute good intial values.
#' @param fastgllvm: an object of class ''fastgllvm''
#' @param target: if non NULL, must be a matrix of loadings: will perform a Procrustes rotation to target
#' @param rescale: if true, both Z and A will be rescaled so that Z%*%t(A) remains unchanged but Z has variance identity.
compute_parameters_initial_values <- function(fastgllvm, target=NULL, rescale=F) {
  pb <- txtProgressBar(style=3, width=40)
  cat(" Initializing: ")
  with(fastgllvm, {
    # TODO: try initializing using gllvm
    Y.transformed <- Y
    if(length(families$id$poisson) > 0) {
      Y.transformed[,families$id$poisson] <- log(Y.transformed[,families$id$poisson] + 0.1) * 1.5
    }

    if(length(families$id$binomial) > 0) {
      Y.transformed[,families$id$binomial] <- (Y.transformed[,families$id$binomial] - 0.5) *6
    }

    NAs <- is.na(Y)
    Y.transformed[NAs] <- 0

    setTxtProgressBar(pb, 0.2)
    cat(" Initializing B...")

    # Initialize XB
    if(dimensions$k >0) {
      B <- t(lm(Y.transformed ~ 0+X)$coef) * (1+colMeans(NAs)) # TODO: ffa should initialize with missing values as well
      XB <- X %*% t(B)
      Y.transformed <- Y.transformed - XB
    } else {
      B <- NULL
      XB <- NULL
    }

    setTxtProgressBar(pb, 0.4, title="Initializing A")
    cat(" Initializing A and phi...")
    # Initialize A and phi
    Y.transformed[NAs] <- NA
    fit.ffa <- ffa(Y.transformed, dimensions$q, maxiter=100, iteratively_update_Psi = T) # TODO: there must be some Z here as well...
    phi <- rep(1, dimensions$p)
    phi[families$id$gaussian] <- fit.ffa$Psi[families$id$gaussian]


    A <- fit.ffa$A
    # transform to target
    if (!is.null(target)) {
      A <- psych::Procrustes(A, target)$loadings
    }
    setTxtProgressBar(pb, 0.9, title="Initializing Z")
    cat(" Initializing Z...")
    # initialize Z
    Z <- fit.ffa$Z
    if(rescale) {
      rescaled <- rescale(Z, A)
      A=rescaled$A
      Z=rescaled$Z
    }
    setTxtProgressBar(pb, 1, title="Initialization Complete.")
    cat(" Initialization complete.")
    close(pb)
    # re-scale zhat and A
    list(A=A, B=B, phi=phi, Z=Z, covZ=cov(Z))
  })
}

# rescale Z to have unit diagonal variance, and A so that ZA remains the same value
# target.cov: the target covariance matrix the new Z must match. If NULL, the variance is the identity matrix.
rescale <- function(parameters, rescale.AB = F, target.cov=NULL, intercept=F) {
  Z <- scale(parameters$Z, scale=F)
  b <- matrix(attr(Z, "scaled:center"), nrow=1)

  if (is.null(target.cov)) {
    C <- chol((t(Z) %*% Z)/nrow(Z))
    Cneg <- solve(C)
    Z <- Z %*% Cneg

  } else {
    Z.c <- chol((t(Z) %*% Z)/nrow(Z))
    target.c <- chol(target.cov)
    C <- solve(target.c) %*% Z.c
    Cneg <- solve(Z.c) %*% target.c
    Z <- Z %*% Cneg

  }
  if (rescale.AB) {
      parameters$A <- parameters$A %*% t(C)
      if(!is.null(parameters$B) && intercept) {
        parameters$B[,1] <- parameters$B[,1]  + as.vector(parameters$A %*% t(b))
      }
  }
  # add the rescaled bias (p.137 in blue book)
  if(!rescale.AB) {
    Z <- t(t(Z) + as.vector(b %*% Cneg))
  }
  parameters$Z <- Z

  parameters
}




# init_Z <- function(Y, A, phi) {
#   Y %*% (A/phi) %*% solve(t(A) %*% (A/phi))
# }
# rescale Z to have unit diagonal variance, and A so that ZA remains the same value
# target: a data matrix whose variance the transformed Z must match. If NULL, the variance is the identity matrix.
rescale.old <- function(Z, A=NULL, target.data=NULL, target.cov=NULL) {
  Z <- scale(Z, scale=F)
  b <- matrix(attr(Z, "scaled:center"), nrow=1)
  if (is.null(target.data) && is.null(target.cov)) {
    #TODO: rescale B as well here? if Z is not centered, recenter it and change B so that XB + ZA remains the same value
    C <- chol((t(Z) %*% Z)/nrow(Z))
    Cneg <- solve(C)
    Z <- Z %*% Cneg
    if (!is.null(A)) A <- A %*% t(C)
  } else {
    Z.c <- chol((t(Z) %*% Z)/nrow(Z))
    if(is.null(target.cov)) {
      target <- scale(target.data, scale=F)
      target.c <- chol(t(target) %*% target/nrow(target))
    } else {
      target.c <- chol(target.cov)
    }

    C <- solve(target.c) %*% Z.c
    Cneg <- solve(Z.c) %*% target.c
    Z <- Z %*% Cneg
    if (!is.null(A)) A <- A %*% t(C)
  }

  # add the rescaled bias (p.137 in blue book)
  Z <- t(t(Z) + as.vector(b %*% Cneg))
  # TODO: TEST THAT: Zt(A) (after centing Z) remains the same
  list(Z=Z, A=A)
}

if(0) {
  devtools::load_all()
  set.seed(1231)
  fg <- gen_fastgllvm(nobs=1000, p=1000, q=20, family=c(rep("poisson", 500), rep("gaussian", 0), rep("binomial", 500)), k=1, intercept=1, miss.prob = 0.5)
  fg <- gen_fastgllvm(nobs=1000, p=100, q=10, family=c(rep("poisson", 0), rep("gaussian", 0), rep("binomial", 100)), k=1, intercept=F, miss.prob = 0)
  init <- compute_parameters_initial_values(fg, target=fg$parameters$A, rescale=F)

  plot(fg$parameters$A, init$A); abline(0,1,col=2)
  points(fg$parameters$Z, init$Z, col=2)
  points(fg$parameters$phi, init$phi, col=3)
  points(fg$parameters$B, init$B, col=4)


  # TEST RESCALE
  A <- scale(matrix(rnorm(100*4), 100, 4), scale=F)
  B <- scale(matrix(rnorm(100*4), 100, 4), scale=F)

  AA <- t(A) %*% A/nrow(A)
  BB <- t(B) %*% B/nrow(B)

  A.c <- chol(AA)
  B.c <- chol(BB)

  all.equal(t(A.c) %*% A.c, AA)
  all.equal(t(B.c) %*% B.c, BB)

  A.unit <- A %*% solve(A.c)
  all.equal(t(A.unit) %*% A.unit/(nrow(A.unit)), diag(4))

  A.B <- A %*% (solve(A.c) %*% B.c)
  all.equal(t(A.B) %*% A.B/nrow(A.B), BB)



  Z <- scale(Z, scale=F)
  Z <- fg$Z
  ZA.before <- Z %*% t(fg$parameters$A)

  ZA.after <- rescale(Z, fg$parameters$A)
  all.equal(ZA.before, ZA.after$Z %*% t(ZA.after$A))

  X <- Z/2
  ZA.after <- rescale(Z, fg$parameters$A, target=X)
  all.equal(ZA.before, ZA.after$Z %*% t(ZA.after$A))
  all.equal(cov(ZA.after$Z), cov(X))

  # learning rate tests
  devtools::load_all()
  igrid <- 1:100
  learning_rate <- initialize_learning_rate(method="spall", 100, learning_rate.args = list(rate=2, end=.01))
  plot(igrid, learning_rate(igrid), main="spall")
  learning_rate <- initialize_learning_rate(method="constant", 100)
  plot(igrid, learning_rate(igrid), main="constant")
  learning_rate <- initialize_learning_rate(method="lin", 100)
  plot(igrid, learning_rate(igrid), main="lin")
  learning_rate <- initialize_learning_rate(method="exp", 100)
  plot(igrid, learning_rate(igrid), main="exp")
  learning_rate <- initialize_learning_rate(method="linexp", 100)
  plot(igrid, learning_rate(igrid), main="linexp")
}
