devtools::load_all()
library(fastfactoranalysis)

genMask <- function(Y, nmiss, type="MCAR"){
  n <- nrow(Y)
  p <- ncol(Y)
  if(type == "MCAR"){
    s <- sample(n, nmiss)
    M <- matrix(0, n, p)
    M[s, 1] <- 1
  } else if(type == "trunc"){
    M <- matrix(0,n, p)
    M[Y[,1] > quantile(Y[,1], (1-nmiss/n)),1] <- 1
  } else {
    stop("imputation type unknown")
  }
  M>0
}

impute <- function(Y, M){
  Ymeans <- colMeans(Y)
  Yimp <- sapply(1:ncol(Y), function(j){
    Yj <- Y[,j]
    Yj[M[,j]] <- Ymeans[j]
    Yj
  })
  Yimp
}


genY <- function(n, p, type="IID", ...){
  if(type=="IID"){
    args <- list(...)
    if(is.null(args$mu)) mu <- rep(0, p) else mu <- args$mu
    if(is.null(args$psi)) psi <- rep(1, p) else psi <- args$psi
    Y <- sapply(1:p, function(j){
      rnorm(n, mu[j], psi[j])
    })
    args <- list(mu=mu)
  }
  if(type=="FA"){
    args <- list(...)
    if(is.null(args$q))stop("Please specify q.")
    if(is.null(args$psi)) psi <- rep(1, p) else psi <- args$psi
    if(is.null(args$mu)) mu <- rep(0, p) else mu <- args$mu
    if(is.null(args$A)) A <- matrix(rnorm(p*q), p, q) else A <- args$A
    if(is.null(args$Z)) Z <- matrix(rnorm(n*q), n, q) else Z <- args$Z
    linpar <- t(A %*% t(Z) + mu)
    Y <- sapply(1:p, function(j){
      rnorm(n, linpar[,j], psi[j])
    })
    args <- list(mu=mu, A=A, linpar=linpar, psi=psi, Z=Z)
  }
  list(Y=Y, args=args)
}

pihat <- function(Y, M){
  colMeans(impute(Y, M))
}

n <- 1000
p <- 3
q <- 1

nmiss <- .1*n
type <- "MCAR"
type <- "trunc"
set.seed(241411)
mu0 <- rep(0, p)
Y0 <- genY(n, p, mu=mu0)$Y
M0 <- genMask(Y0, nmiss, type=type)
Yimp <- impute(Y0, M0)

sims <- t(sapply(1:100, function(seed){
  set.seed(12412+seed)
  Y <- genY(n, p)$Y
  M <- genMask(Y, nmiss, type=type)
  pihat(Y, M)
}))



ib <- function(Y, M, seed=3123, H=10){
  theta <- target <- pihat(Y,M)
  for(i in 1:3){
    set.seed(seed)
    exp <- rowMeans(sapply(1:H, function(i){
      Y <- genY(n, p, mu=theta)$Y
      M <- genMask(Y, nmiss, type=type)
      pihat(Y, M)
    }))
    theta <- theta + (target - exp)
  }
  theta
}


sims.ib <- t(sapply(1:100, function(seed){
  set.seed(12412+seed)
  Y <- genY(n, p)$Y
  M <- genMask(Y, nmiss, type=type)
  ib(Y, M, seed=seed+4442323, H=10)
}))

boxplot(cbind(sims, sims.ib)); abline(0,0)


theta <- target <- pihat(Y0,M0)

for(i in 1:10){
  set.seed(1230)
  exp <- rowMeans(sapply(1:1, function(i){
    Y <- genY(n, p, mu=theta)$Y
    M <- genMask(Y, nmiss, type=type)
    pihat(Y, M)
  }))
  theta <- theta + (target - exp)
  print(theta)
}

get_Z <- function(Y, A, psi){
  scale(Y, scale=F) %*% diag(psi^-1) %*% A %*% solve(t(A) %*% diag(psi^-1) %*% A + diag(q))
}
# NEW SETTING
# with correlation, learning low rank, known mask


theta_to_par <- function(theta, p, q){
  A <- matrix(theta[1:(p*q)], p, q)
  psi <- theta[(p*q+1):(p*q+p)]
  mu  <- theta[(p*q+p + 1): length(theta)]
  list(A=A, psi=psi, mu=mu)
}
par_to_theta <- function(A, psi, mu){
  c(as.vector(A), psi, mu)
}

generator <- function(theta, n, p, q, missType){
  par <- theta_to_par(theta, p, q)
  Y <- genY(n, p, q=q, type="FA", A=par$A, psi=par$psi, mu=par$mu)$Y
  M <- genMask(Y, nmiss=nmiss, type=missType)
  list(Y=Y, M=M)
}

estimator <- function(Y, q){
  Yimp <- impute(Y$Y, Y$M)
  fit <- ffa(Yimp, nfactors=q)
  par_to_theta(A=fit$loadings, psi=fit$communalites, mu=colMeans(Yimp))
}


# estimator must take a list named Y as first parameter, the rest in est.args.
# generator must take theta, a vector, as first parameter, the rest in gen.args
ib <- function(Y, estimator, generator, H=1, step=0.1, maxit=100, gen.args=list(), est.args=list(), seed=12451, verbose=F, save=F){
 target <- theta <- do.call(estimator, c(list(Y=Y), est.args))
 if(save) hist <- matrix(0, maxit, length(target))
 for(i in 1:maxit){
   exp <- lapply(1:H, function(h){
     set.seed(seed+h)
     Y <- do.call(generator, c(list(theta=theta), gen.args))
     est <- do.call(estimator, c(list(Y=Y), est.args))
     est
   })
   exp <- Reduce("+", exp)/length(exp)
   theta <- theta + step *(target - exp)
   theta[is.na(theta)] <- -Inf
   cat(i, "\n", theta[c(1,6)])
   if(save) hist[i,] <- theta
 }
 list(theta=theta, target=target, save= if(save) hist else NULL)
}



n <- 10000
p <- 10
q <- 2
nmiss <- .2*n
type <- "MCAR"
missType <- "trunc"
set.seed(1023)
dat <- genY(n, p, q=q, type="FA", mu=rnorm(p))
fit <- ffa(dat$Y, nfactors=q)
dat$M <- genMask(dat$Y, nmiss = nmiss, type=missType)


fit.ib <- ib(dat[c("Y", "M")], estimator, generator, est.args=list(q=q),
             gen.args=list(n=n, p=p, q=q, missType=missType), maxit=10, step=.1, seed=12312, H=10)

tini <- theta_to_par(fit.ib$target, p, q)
that <- theta_to_par(fit.ib$theta, p, q)

plot(dat$args$A, psych::Procrustes(tini$A, dat$args$A)$loadings)
points(dat$args$A, psych::Procrustes(that$A, dat$args$A)$loadings, col=2)
abline(0,1,col=2)
abline(0,-1,col=2)

# #########################
# Part 3: ADDING THE MASK AND LEARNING EVERYTHING
# LOW RANK APPROXIMATION WITH FACTOR ANALYSIS


# pmiss: number of variables with missing values
generator <- function(theta, n, p, q){
  par <- theta_to_par(theta, p, q)
  Z <- matrix(rnorm(n*q), n ,q)
  linpar <- t(par$A %*% t(Z) + par$mu)
  Y <- sapply(1:(p/2), function(j){
    rnorm(n, linpar[,j], par$psi[j])
  })
  M <- sapply((p/2+1):p, function(j){
    rbinom(n, 1, 1/(1+exp(-linpar[,j])))
  })>0
  list(Y=Y, M=M)
}

estimator <- function(Y, q){
  p <- ncol(Y$Y)*2
  no.missing <- c(rep(T, p/2), colMeans(Y$M)!=0)
  Yimp <- impute(Y$Y, Y$M)
  fit <- ffa(cbind(Yimp, Y$M*4)[, no.missing], nfactors=q)

  # Add rows of 0 for the missing variables
  A <- matrix(0, p, q)
  A[no.missing, ] <- fit$loadings
  psi <- rep(.1, p)
  psi[no.missing] <- fit$communalites
  psi[(p/2+1):p] <- .1
  mu <- rep(-Inf, p)
  mu[no.missing] <- fit$means
  par_to_theta(A=A, psi=psi, mu=mu)
}

# set.seed(1231)
# n <- 1000
# p <- 20
# q <- 2
# A <- matrix(rnorm(p*q*2), p*2, q)
# A[c(1, p/2+1)] <- 2
# mu <- rnorm(p*2)
# mu[(p+2):(2*p)] <- -Inf
# psi <- rep(1, p*2)
# par <- list(A=A, psi=psi, mu=mu)




# gen data from another mechanism
n <- 1000
p <- 20
q <- 1
nmiss <- .1*n
missType <- "trunc"
set.seed(10243)

A <- matrix(rnorm(p*q), p, q)*1.5
A[1:2] <- 2

dat0 <- genY(n, p, q=q, type="FA", mu=rnorm(p), A=A)
dat0$M <- genMask(dat0$Y, nmiss = nmiss, type=missType)

hist(dat0$Y[dat0$M])
hist(dat0$Y[,1])

# recover now everything
fit <- ib(dat0, estimator, generator, H=5, step=.5, maxit=500,
   gen.args=list(n=n, p=p*2, q=q),
   est.args=list(q=q), save=T)

ts.plot(fit$save)

par.ib <- theta_to_par(fit$theta, p*2, q)
par.ib.avg <- theta_to_par(colMeans(fit$save[50:100,]), p*2, q)
par.tar <- theta_to_par(fit$target, p*2, q)


# true
A <- dat0$args$A

A.ib <- psych::Procrustes(par.ib$A[1:p,], A)$load
A.tar <- psych::Procrustes(par.tar$A[1:p, ], A)$load

plot(A, A.ib)
points(A, A.tar, col=2)
abline(0,1,col=2)

plot(dat0$args$mu, par.ib$mu[1:p])
points(dat0$args$mu, par.tar$mu[1:p], col=2); abline(0,1,col=2)

