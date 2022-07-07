# Test that we can estimate zstar

# TODO: test that this converges
set.seed(2131)
n <- 1000
p <- 100
q <- 5
A <- matrix(rnorm(p*q), p, q)
B <- matrix(rnorm(p*q), p, q)
family <- "poisson"
family <- "binomial"

family <- list(poisson=c(1:50), binomial=c(51:100))

# f <- gen_fastgllvm(n=1000, p = 8, q=2, k=4, intercept=T, phi=c(rep(1, p/2), (1:(p/2))/1), family=list(binomial=1:(p/4), poisson=(p/4+1):(p/2), gaussian=(p/2+1):p))
set.seed(123)
f <- gen_fastgllvm(n=1000, p =p, q = q, k=1, A=A, phi=rep(1, p), intercept=F, family=family)
Zstar <- compute_zstar(f$Y, f$A, f$phi, f$linpar$XB, f$families, f$dims, save=F, thresh=1e-3, start=f$Z)
# Zstar <- compute_zstar(f$Y, f$A, f$phi, f$linpar$XB, f$families, f$dims, save=F, thresh=1e-3)

plot(f$Z, Zstar$Zstar, main = paste("Converged in ", Zstar$niter, " iterations."))

# test that the psi are computed correctly
list2env(f, .GlobalEnv)
# Test that the

# newton for A
A <- t(sapply(1:p, function(j){
  lm(log(Y[,j]+0.01)~0+Z, offset=X%*%t(B)[,j])$coef
}))*1.5
# A <- matrix(0, p, q)
hist <- list()
for(i in 1:20){
  A <- A - compute_psi_AB(Y, Z, X, B, A, phi, families, lambda = 2)$psi_A
  hist[[i]] <- as.vector(A)
}
hist <- do.call(rbind, hist)
ts.plot(hist)
abline(h=f$A, col=2)

plot(f$A, A)

A <- NULL
B <- NULL
phi <- NULL
set.seed(52345)
Z <- generate_Z()[[1]]

pi_0 <- compute_pi(Y=Y, Z=Z_0_hat, X=X, A=A, B=B, phi=phi, families=families, maxit=100)
A <- A + 1/i*(pi_0$A - pi_gen$A)
B <- B + 1/i*(pi_0$B - pi_gen$B)
phi <- phi - 1/i*(pi_0$phi - pi_gen$phi)

# A <- matrix(rnorm(p*q), p, q)*.1
# B <- matrix(rnorm(p*ncol(X)), p, ncol(X)) *.1
# phi <- rep(1, p)
sol <- compute_pi(Y=f$Y, Z=Z, X=f$X, A=f$A, B=f$B, phi=f$phi, families=f$families, maxit=100)
A <- sol$A
B <- sol$B
phi <- sol$phi


for(i in 1:10){

  Z_0_hat <- compute_zstar(Y, A, phi, X %*% t(B), families, start=Z, dims, scale=F, maxit=100)$Zstar

  Z_gen <- generate_Z()[[1]]
  linpar <- compute_linpar(Z_gen, A, X=X, B=B)
  Y_gen <- generate_y(linpar, phi, families, n, p, q)

  Z_gen_hat <- compute_zstar(Y_gen, A, phi, X %*% t(B), families, start=Z_gen, dims, scale=F, maxit=100, thresh=1e-20)$Zstar

  pi_gen <- compute_pi(Y=Y_gen, Z=Z_gen_hat, X=X, A=A, B=B, phi=phi, families=families, maxit=100)



  A.proc <- psych::Procrustes(f$A, A)$loadings
  plot(f$A, A.proc)
  abline(0,1,col=2)

  i <- i+1
}
plot(f$Z, Z)
