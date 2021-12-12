# Functions to compute the ARE
devtools::load_all()
# b and c functions
bfunc <- function(natpar){
  natpar - log(sigmoid(natpar))
}
cfunc <- function(Y, psi){
  Y * 0
}


# conditional density of Y (n * p matrix) given Z (n * q matrix), returns an n-vector

fYZ <- function(Y, Z, par){
  natpar <- Z %*% t(par$A)
  dens <- exp(t(t(Y * natpar - bfunc(natpar))/par$Psi) + cfunc(Y, par$Psi))
  apply(dens, 1, prod)
}

# same as fYZ but with Y a vector, also there is no cfunc... TODO

fYiZ <- function(Yi, Z, par){
  natpar <- Z %*% t(par$A)
  dens <- exp(t(t(t(Yi * t(natpar)) - bfunc(natpar))/par$Psi))
  apply(dens, 1, prod)
}


# derviative of the log conditional distribution w.r.t A.
# Returns an n * (p*q) matrix, where each row is a p*q vector,
# which corresponds to the rows of the corresponding p *q matrix of partial derivatives that have been concatenated into a vector. Since there is one such vector per observations, we obtain $n * (p*q)$ matrix in the end

dlogfYZ <- function(Y, Z, par){
  n <- nrow(Y)
  p <- nrow(par$A)
  q <- ncol(par$A)
  natpar <- Z %*% t(par$A)
  t(sapply(1:n, function(i){
    Zimat <- matrix(Z[i,], nrow=p, ncol=q, byrow=T)
    matderiv <- (Y[i,] - sigmoid(natpar[i,])) * Zimat
    vecderiv <- as.vector(t(matderiv))
  }))
}

# same as dlogfYZ but for a single observation
dlogfYiZi <- function(Yi, Zi, par){
  natpar <- as.vector(Zi %*% t(par$A))
  as.vector(t((Yi - sigmoid(natpar)) %*% t(Zi)))
}

# same as dlogfYZ but with vector Y and matrix Z
dlogfYiZ <- function(Yi, Z, par){
  p <- nrow(par$A)
  q <- ncol(par$A)
  n <- nrow(Z)
  natpar <- Z %*% t(par$A)
  t(sapply(1:n, function(i){
    Zimat <- matrix(Z[i,], nrow=p, ncol=q, byrow=T)
    matderiv <- (Yi - sigmoid(natpar[i,])) * Zimat
    vecderiv <- as.vector(t(matderiv))
  }))
}

# Returns D() for each of the observations; returns an n-vector
compute_D <- function(Y, par, nD, seed=NULL){
  n <- nrow(Y)
  q <- ncol(par$A)

  if(is.null(seed))set.seed(seed)
  if(n==1){
    D <- mean(sapply(1:nD, function(na){
      Z <- gen_Z(n = n, q=q)
      fYZ(Y, Z, par)
    }))
  } else {
    D <-
      rowMeans(sapply(1:nD, function(na){
        Z <- gen_Z(n = n, q=q)
        fYZ(Y, Z, par)
      }))
  }
  D
}

# Compute N() for each of the observations, returns an n-vector
compute_N <- function(Y, par, nN, seed=NULL){
  n <- nrow(Y)
  q <- ncol(par$A)

  if(is.null(seed))set.seed(seed)
  t(sapply(1:n, function(i){
    Z <- gen_Z(nN, q)
    # dlogFYZ returns an n vector
    colMeans(dlogfYiZ(Y[i,], Z, par) * fYiZ(Y[i,], Z, par))
  }))
}

# finally, compute the scores as a n * (p*q) matrix:
compute_scores <- function(Y, par, nD, nN, seed=NULL){
  N <- compute_N(Y, par, nN, seed)
  D <- compute_D(Y, par, nD, seed)
  N/D
}

# Compute Jacobian:

jacobian <- function(score, theta, eps=1e-3, ...){
  pq <- length(theta)
  jacobian <- matrix(0, pq, pq)
  for(i in 1:pq){
    delta <- rep(0, pq)
    delta[i] <- eps
    jacobian[,i] <- (score(theta+delta, ...) - score(theta - delta, ...))/(2*eps)
  }
  jacobian
}



# or just one score from one observation, as a (p*q) vector:

compute_score <- function(Yi, par, nD, nN, seed=NULL){
  Y <- t(as.vector(Yi))
  N <- compute_N(Y, par, nN, seed=NULL)
  D <- compute_D(Y, par, nD, seed=NULL)
  N/D
}

get_I <- function(scores){
  # center scores
  scores <- scale(scores, scale=F)
  I <- t(scores_Y) %*% scores_Y
  I <- I/nrow(scores_Y)
  I
}

compute_I <- function(par, n, nD, nN, seed=NULL){
  p <- nrow(par$A)
  q <- ncol(par$A)
  gllvm_dat <- gen_gllvm(n, p, q, k=0, family="bernoulli", par=par)
  Y <- gllvm_dat$Y
  scores <- compute_scores(Y, par, nD, nN, seed=NULL)
  I <- get_I(scores)
  Ineg <- solve(I)
  liste(I=I, Ineg =Ineg)
}

# transforms theta from vector to matrix and vice versa
theta_to_vec <- function(theta) as.vector(theta)
theta_to_mat <- function(theta, p, q) matrix(theta, p, q)

# compute all possible outcomes of Y1, ..., Yp
compute_outcomes <- function(p){
  t(sapply(0:(2**p-1),function(x){ as.integer(intToBits(x))})[p:1, ])
}

# Compute the scores at Y
compute_scores_Y <- function(Y, par, nD, nN, seed=NULL, outcomes_scores=NULL){
  n <- nrow(Y)
  p <- nrow(par$A)
  q <- ncol(par$A)

  num_outcomes <- 2**p

  # only compute the scores for the different outcomes
  if(is.null(outcomes_scores)){
    outcomes <- compute_outcomes(p)
    outcomes_scores <- compute_scores(outcomes, par, nD, nN, seed=seed)
  }

  bin <- 2^((p-1):0)
  Y.outcomes <- as.vector(Y %*% bin)
  Y.outcomes.freq  <- table(Y.outcomes)

  scores_Y <- matrix(0, n, p*q)

  for(i in 1:num_outcomes){
    scores_Y[which(Y.outcomes == (i-1)), ] <- matrix(outcomes_scores[i, ],
                                                     nrow=Y.outcomes.freq[i],
                                                     ncol=p*q,
                                                     byrow = T)
  }
  scores_Y
}


# Functions for our method:


compute_Z <- function(Y, par){
  K <- compute_K(par)
  Z <- Y %*% t(K)
}

compute_profile_scores <- function(Y, par, nE=nrow(Y), seed=NULL){
  p <- nrow(par$A)
  q <- ncol(par$A)

  Zstar <- compute_Z(Y, par)
  deriv <- dlogfYZ(Y, Zstar, par)
  # compute expectation under par_true
  if(!is.null(seed)) set.seed(seed)
  Y0 <- gen_gllvm(nE, p, q, family="bernoulli", par=par)$Y
  Zstar0 <- compute_Z(Y0, par)
  deriv0 <- dlogfYZ(Y0, Zstar0, par)
  deriv0.mean <- colMeans(deriv0)

  t(t(deriv) - deriv0.mean)
}

variance_scores <- function(Y, par, nE){
  scores <- compute_profile_scores(Y, par, nE)
  n <- nrow(scores)
  scores <- scale(scores, scale=F)
  t(scores) %*% scores/n
}

score_to_derive <- function(Y, theta, par, seed=NULL){
  p <- nrow(par$A)
  q <- ncol(par$A)
  par2 <- par
  par2$A <- theta_to_mat(theta, p, q)
  colMeans(compute_profile_scores(Y, par2, nE, seed=seed))
}


simulate_ARE <- function(par1=0){
  # Simulate:
  set.seed(3523)
  p <-  5
  q <- 1
  n <- 1e4
  nN <- nD <- 10000
  par <- gen_par(p, q, family="bernoulli")
  par$A <- matrix(c(par1, -1, -1, 1, 1, .2, 1, .4)[1:p], p, 1)

  gllvm <- gen_gllvm(n, p, q, family = "bernoulli", par=par)
  Y <- gllvm$Y
  par <- gllvm$par
  scores_Y <- compute_scores_Y(Y, par, nD, nN, seed=4345)

  boxplot(scores_Y[sample(n, 100, repl=T),], main="Scores at the true value of the parameter.")
  points(1:(p*q), colMeans(scores_Y), col=2, pch=1, pty=3)
  abline(h=0, col=2)

  I <- get_I(scores_Y)

  Ineg <- solve(I)
  Ineg


  nE <- n
  theta <- theta_to_vec(par$A)


  # now we need the derivative

  A <- - jacobian(score_to_derive, theta, eps=.5, Y=Y, par=par, seed=532)
  # A2 <- - numDeriv::jacobian(score_to_derive, theta, method.args=list(eps=1e-1))
  B <- variance_scores(Y, par, nE)
  # now compute the avar
  Aneg <- solve(A)
  Avar <- Aneg %*% B %*% t(Aneg)

  ARE = sum(diag(Ineg))/sum(diag(Avar))
  ARE
}

par1 <- seq(-.5,.5,l=11)

Ares <- sapply(par1, function(par)simulate_ARE(par))
plot(Ares)
