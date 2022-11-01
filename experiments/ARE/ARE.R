# Functions to compute the ARE
devtools::load_all()
# b and c functions
bfunc <- function(natpar){
  natpar - log(sigmoid(natpar))
}
cfunc <- function(Y, phi){
  Y * 0
}

sigmoid <- function(x) 1/(1+exp(-x))


# conditional density of Y (n * p matrix) given Z (n * q matrix), returns an n-vector

fYZ <- function(Y, Z, par){
  natpar <- Z %*% t(par$A)
  dens <- exp(t(t(Y * natpar - bfunc(natpar))/par$phi) + cfunc(Y, par$phi))
  apply(dens, 1, prod)
}

# same as fYZ but with Y a vector, also there is no cfunc... TODO

fYiZ <- function(Yi, Z, par){
  natpar <- Z %*% t(par$A)
  dens <- exp(t(t(t(Yi * t(natpar)) - bfunc(natpar))/par$phi))
  apply(dens, 1, prod)
}


# deriviative of the log conditional distribution w.r.t A.
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


# or just one score from one observation, as a (p*q) vector:
compute_score <- function(Yi, par, nD, nN, seedD=NULL, seedN=NULL){
  Y <- t(as.vector(Yi))
  N <- compute_N(Y, par, nN, seed=seedN)
  D <- compute_D(Y, par, nD, seed=seedD)
  N/D
}

# cl is the cluster
compute_scores_parallel <- function(cl, Y, par, nD, nN, seed=NULL){
  parallel::clusterExport(cl, c("Y", "par", "nD", "nN", "seed", "gen_Z"), envir = environment())
  t(parallel::parSapply(cl, 1:nrow(Y), function(i) {
    compute_score(Y[i,], par, nD, nN, seedD=seed+i, seedN=nrow(Y)+seed+i)
  }))
}


gen_Z <- function(nobs, q){
  matrix(stats::rnorm(nobs*q), nobs, q)
}


# Compute Jacobian: (not used anymore)

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




get_I <- function(scores){
  # center scores NOOO?! dont!
  scores <- scale(scores, scale=F)
  I <- t(scores) %*% scores
  I <- I/nrow(scores)
  solve(I)
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
  Z <- Y %*% K
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


compute_K <- function(par){
  A <- par$A
  phi <- par$phi
  p <- nrow(A)
  q <- ncol(A)
  t(solve(t(A) %*% diag(phi**-1) %*% A + diag(q)) %*% t(A) %*% diag(phi**-1))
}



compute_psi_AB <- function(Y, ZX, phi, linpar_bprime, Miss){
  if (!is.null(Miss)) {
    Y[Miss] <- 0 # this is correct... check the NA rmd file for examples
    linpar_bprime[Miss] <- 0
  }
  t(sapply(1:nrow(Y), function(i){
    (Y[i,] - linpar_bprime[i,])/phi * ZX[i,]
  }))
}

compute_psi_AB_old <- function(Y, ZX, phi, linpar_bprime, Miss){
  if (!is.null(Miss)) {
    Y[Miss] <- 0 # this is correct... check the NA rmd file for examples
    linpar_bprime[Miss] <- 0
  }
  (t(Y - linpar_bprime)/phi) %*% ZX
}



compute_I_PROMES <- function(scores_Y, gllvm) {

  PSI <- compute_scores_PROMES(gllvm)
  PSI_C <- scale(PSI, scale=F) # for nobs large enough this is justified because the score needs to be centered by its expectation, which can be obtained by taking the empirical mean...
  # for the "A", centering is also justified because the expecation of the MLE scores is 0.
  scores_Y_C <- scale(scores_Y, scale=F)

  A <- t(PSI_C) %*% scores_Y_C
  B <- t(PSI_C) %*% PSI_C

  Aneg <- solve(A)
  list(Avar = Aneg %*% B %*% t(Aneg) * nrow(gllvm$Y), scores=PSI)
}

compute_I_SPROMES <- function(scores_Y, gllvm) {

  PSI <- scale(compute_scores_SPROMES(gllvm), scale=F)

  PSI_C <- scale(PSI, scale=F) # this is justified because the score needs to be centered by its expectation, which can be obtained by taking the empirical mean...
  # for the "A", centering is also justified because the expecation of the MLE scores is 0.
  scores_Y_C <- scale(scores_Y, scale=F)


  # # compute its expectation NO NEED, just center ()
  # PSI_E <- lapply(1:1000, function(seedi){
  #   if(!is.null(seed)) set.seed(seed + seedi)
  #   gllvmsim <- gen_fastgllvm(n, p, q, family = "binomial", A= gllvm$parameters$A, intercept = F, k=0)
  #   compute_scores_SPROMES(gllvmsim)
  # })
  # PSI_E <- Reduce("+", PSI_E)/length(PSI_E)
  #
  # PSI_C <- PSI - PSI_E
  #
  # A <- t(scale(PSI, scale=F)) %*% scale(scores_Y, scale=F)
  # B <- t(scale(PSI_C, scale=F)) %*% scale(PSI_C, scale=F)
  # B2 <- t(scale(PSI, scale=F)) %*% scale(PSI, scale=F)

  A <- t(PSI_C) %*% scores_Y_C
  B <- t(PSI_C) %*% PSI_C

  Aneg <- solve(A)
  list(Avar=Aneg %*% B %*% t(Aneg) * nrow(gllvm$Y), scores = PSI)
}


compute_scores_PROMES <- function(gllvm) {
  # for A:
  zstar <- compute_zstar(gllvm$Y, A=gllvm$parameters$A, phi=gllvm$parameters$phi, B=gllvm$parameters$B, X=gllvm$X, families = gllvm$families, start = gllvm$parameters$Z)
  zstar <- zstar$Zstar
  # zstar <- gllvm$parameters$Z
  # zstar <- compute_Z(Y = gllvm$Y, par = gllvm$parameters)

  ZX <- cbind(zstar, gllvm$X)
  linpar <- compute_linpar(Z = zstar, A = gllvm$parameters$A, X = gllvm$parameters$X, B = gllvm$parameters$B)$linpar
  linpar_bprime <- compute_linpar_bprime(linpar, gllvm$families)
  linpar_bprimeprime <- compute_linpar_bprimeprime(linpar, gllvm$families)
  t(sapply(1:nrow(gllvm$Y), function(i) ((gllvm$Y[i,]-linpar_bprime[i,]) * ZX[i,])/gllvm$parameters$phi))
}

compute_scores_SPROMES <- function(gllvm) {
  zstar <- compute_zstar(gllvm$Y, A=gllvm$parameters$A, phi=gllvm$parameters$phi, B=gllvm$parameters$B, X=gllvm$X, families = gllvm$families, start = gllvm$parameters$Z)
  zstar <- zstar$Zstar

  ZX <- cbind(zstar, gllvm$X)

  t(sapply(1:nrow(gllvm$Y), function(i) (gllvm$Y[i,] * ZX[i,])/gllvm$parameters$phi))
}

sim_ARE <- function(cl, p, seed) {
  set.seed(seed)
  par <- generate_parameters(A = matrix(c(-2, -1, rep(0, p-4), 1,2), p, 1), B=NULL, phi=NULL, p=p, q=q, k=0)
  set.seed(seed)
  par <- generate_parameters(A = matrix(c(-2, -1, runif(p-4, -.5, .5), 1,2), p, 1), B=NULL, phi=NULL, p=p, q=q, k=0)
  set.seed(seed)
  par <- generate_parameters(A = matrix(runif(p, -1, 1)), B=NULL, phi=NULL, p=p, q=q, k=0)
  set.seed(seed)
  gllvm <- gen_fastgllvm(n, p, q, family = "binomial", A= par$A, intercept = F, k=0)
  zstar <- compute_zstar(gllvm$Y, A=gllvm$parameters$A, phi=gllvm$parameters$phi, B=gllvm$parameters$B, X=gllvm$X, families = gllvm$families, start = gllvm$parameters$Z)$Zstar

  scores_Y <- compute_scores_parallel(cl, gllvm$Y, gllvm$parameters, nD, nN, seed=seed+124312)

  I_MLE <- get_I(scores_Y)
  I_PROMES <- compute_I_PROMES(scores_Y, gllvm)
  I_SPROMES <- compute_I_SPROMES(scores_Y, gllvm)

  list(
    PROMES1 = mean(diag(I_MLE %*% solve(I_PROMES$Avar))),
    PROMES2 = mean(diag(I_MLE))/mean(diag(I_PROMES$Avar)),
    SPROMES1 = mean(diag(I_MLE %*% solve(I_SPROMES$Avar))),
    SPROMES2 = mean(diag(I_MLE))/mean(diag(I_SPROMES$Avar)),
    I_MLE = I_MLE,
    I_PROMES = I_PROMES,
    I_SPROMES = I_SPROMES
  )
}

library(parallel)
# Detect the number of available cores and create cluster
cl <- parallel::makeCluster(detectCores()-2)
# Run parallel computation
parallel::clusterExport(cl, ls())
parallel::clusterExport(cl, c("gen_fastgllvm", "compute_zstar", "Z_update"))
# Simulate:

q <- 1
n <- 1e4
nN <- nD <- 1e3

ARE <- lapply(c(5, 10, 20, 30, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500), function(p){
  cat("\nComputing for p=", p)
  sim_ARE(cl, p, seed=1231+p)}
  )

saveRDS(ARE, file="./Rmd/ARE/ARE_n_1e4_seed1231.rds")
