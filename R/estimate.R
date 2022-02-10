get_K <- function(A, phi){
  p <- nrow(A)
  q <- ncol(A)
  Ap <- A/phi
  solve(t(A)%*% Ap + diag(1, q)) %*% t(Ap)
}

#' Returns the (Monte-Carlo) approximations of E[Y] and E[Y Y]K.
#' @param A: the matrix of loadings
#' @param B: the matrix of fixed effect coefficients
#' @param X: the matrix of covariates
#' @param Z: a list of H generated Z samples
#' @param K: the matrix so that Z = K Y
#' @param family: an object of class "family"
get_expectations <- function(A, B, X, Z, K, phi, family){
  EY <- lapply(Z, function(Zh)family$linkinv(compute_natpar(A, B, Zh, X)))
  EYYK <- lapply(EY, function(EYh){
    EYh.c <- scale(EYh, scale=F) # TODO HEREEEEEEEEEEEE
    # EYYh.diag <- colMeans(EYh.c^2)
    EYh.var <- colMeans(family$variance(EYh)) * phi # TODO: check that we need to multiply by phi... or how to model overdispersion here
    EYYKh <- t(EYh.c) %*% (EYh.c %*% t(K/nrow(EYh.c)))
    # The diagonal of EYZi is wrong. We need to replace its current value (EYYi.diag) by its true value (EYi.var).
    # However, we do not want to compute EYY, so we remove its effect after the multiplication by t(K).
    # t(K) is pre-multiplied because resp.var and EYYi.diag are vectors; it would be post-multiplied if they were
    # diagonal matrices.
    EYYKh <- EYYKh + t(K) * (EYh.var) # the second element is equal to diag(EYh.var) %*% t(K)
    EYYKh
  })
  list(
    EY = Reduce("+", EY)/length(EY),
    EYYK = Reduce("+", EYYK)/length(EYYK)
  )
}

#' Returns the psi function approximated by H monte carlo samples
#' @param Y: uncentered data
#' @param Y.c: Y that are CENTERED, they must have an attribute named `scaled:center` that represent the removed mean
#' @param A: the matrix of loadings
#' @param B: the matrix of fixed effect coefficients
#' @param phi: the vector of response-specific scale parameters
#' @param X: the matrix of covariates
#' @param family: an object of class "family"
#' @param generate_Z a function returned by generate_Z_functionfactory that returns a list
get_Psi <- function(Y, Y.c, A, B, phi, X, family, generate_Z){
  stopifnot(!is.null(attr(Y.c, "scaled:center")))
  p <- nrow(A)
  q <- ncol(A)
  K  <- get_K(A, phi)
  Z  <- generate_Z()
  Exp <- get_expectations(A, B, X, Z, K, phi, family)

  Psi <- list(
    A = t(Y.c) %*% (Y.c %*% t(K/nrow(Y.c))) - Exp$EYYK,
    B = t(Y - Exp$EY) %*% X/nrow(X),
    phi = rep(0, p)  # TODO: implement
  )

  Psi
}

#' Returns a fitted object of class "fastgllvm".
#' @param Y: a matrix of dimensions n x q
#' @param q: the number of latent variables
#' @param X: either 0 (no covariates, no intercept), 1 (no covariates but an intercept), or a matrix of covariates. If the latter and an intercept is desired, it must be included in X as a column of 1.
#' @param family: either one of ("gaussian", "poisson", "binomial"), the function name or the corresponding object of class "family" obtained by calling the function.
#' @param family: one of "SA" or "SP"
#' @param H: how many samples of Z to draw: if in (0,1), then a batch method is used and only that proportion of data is used to estimate the model.
#'
#' @export
fastgllvm_old <- function(Y, q=1, X=0, family=binomial(), method="SA", H=1, A.init=NULL, B.init=NULL, phi.init=NULL, maxit=250, tol=1e-5, learning_rate = ifelse(method=="SA", "exp", "constant"), learning_rate.args = NULL, verbose=T){
  stopifnot(is.matrix(Y))
  n <- nrow(Y)
  p <- ncol(Y)

  if(length(X) == 1){
    stopifnot(X == 0 | X == 1)
    X <- matrix(X, n, 1)
  }

  k <- ncol(X)

  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  if(is.character(learning_rate)){
    learning_rate <- ff_learning_rate(method=learning_rate, maxit=maxit, learning_rate.args = learning_rate.args)
  }
  learning_rate.seq <- learning_rate(1:maxit)

  if(!is.null(A.init)) A <- A.init else A <- diag(1, p, q)
  if(!is.null(B.init)) B <- B.init else B <- matrix(0, p, k)
  if(!is.null(phi.init)) phi <- phi.init else phi <- rep(1, p)

  print(A)
  A.hist <- matrix(0, maxit, p*q)
  A.hist[1,] <- A
  B.hist <- matrix(0, maxit, p*k)
  B.hist[1,] <- B
  phi.hist <- matrix(0, maxit, p)
  phi.hist[1,] <- phi
  crit.hist <- rep(0, maxit)


  Y.c <- scale(Y, scale=F)
  generate_Z <- generate_Z_functionfactory(n, q, method=method, H=H)

  i <- 1
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
    A.hist[i, ] <- as.vector(A)
    B.hist[i, ] <- as.vector(B)
    phi.hist[i, ] <- as.vector(phi)

    crit.hist[i] <- learning_rate.seq[i] * norm(Psi$A)/(p*q)
    if(crit.hist[i] < tol) converged <- TRUE

    if(verbose)cat("\ni: ", i, " - norm:", crit.hist[i], " - learning rate:", learning_rate.seq[i])
    # check if the criterion is small enough to jump to the next "repetition", where the learning rate increases again
    if(i < maxit && crit.hist[i] < tol){
      # fill in the histories
      A.hist <- A.hist[1:i,]
      B.hist <- B.hist[1:i,]
      phi.hist <- phi.hist[1:i,]
    }
  }
  structure(list(
    Y=Y,
    X=X,
    A=A,
    B=B,
    phi=phi,
    family=family,
    n=n,
    p=p,
    q=q,
    k=k,
    A.hist=A.hist,
    B.hist=B.hist,
    phi.hist=phi.hist,
    crit.hist=crit.hist,
    learning_rate.seq = learning_rate.seq),
  class="fastgllvm")
}


#' Returns a fitted object of class "fastgllvm".
#' @param Y: a matrix of dimensions n x q
#' @param q: the number of latent variables
#' @param X: either 0 (no covariates, no intercept), 1 (no covariates but an intercept), or a matrix of covariates. If the latter and an intercept is desired, it must be included in X as a column of 1.
#' @param family: either one of ("gaussian", "poisson", "binomial"), the function name or the corresponding object of class "family" obtained by calling the function.
#' @param family: one of "SA" or "SP"
#' @param H: how many samples of Z to draw: if in (0, 1), then a batch method is used and only that proportion of data is used to estimate the model.
#'
#' @export
fastgllvm <- function(Y, q=1, X=0, family=binomial(), method="SA", H=1, A.init=NULL, B.init=NULL, phi.init=NULL, maxit=250, tol=1e-5, learning_rate = NULL, learning_rate.args = NULL, verbose=T){
  stopifnot(is.matrix(Y))
  n <- nrow(Y)
  p <- ncol(Y)

  if(length(X) == 1){
    stopifnot(X == 0 | X == 1)
    X <- matrix(X, n, 1)
  }

  k <- ncol(X)

  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  if(!is.null(A.init)) A <- A.init else A <- diag(1, p, q)
  if(!is.null(B.init)) B <- B.init else B <- matrix(0, p, k)
  if(!is.null(phi.init)) phi <- phi.init else phi <- rep(1, p)

  hist <- list(
    A = matrix(as.vector(A), 1, p*q),
    B = matrix(as.vector(B), 1, p*k),
    phi = matrix(phi, 1, p),
    crit = rep(0, 1)
  )

  fg <- new_fastgllvm(Y, X, A, B, phi, family, hist)
  validate_fastgllvm(fg)

  fg <- fit_fastgllvm(fg, method,
               H=H,
               maxit=maxit,
               tol=tol,
               learning_rate = learning_rate,
               learning_rate.args = learning_rate.args,
               verbose=T)
  fg
}
