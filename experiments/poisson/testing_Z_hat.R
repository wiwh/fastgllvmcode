rexpfam <- function(n, mu, family, scale=1){
  switch(family$family,
         "poisson"=rpois(n, mu),
         "binomial"=rbinom(n,1, mu),
         "gaussian"=rnorm(n, mu, scale))
}

genY <- function(A){
  browser()
  z <- matrix(rnorm(n*q), n, q)
  lp <- z %*% t(A)
  y <- apply(lp, 1:2, function(ij) rexpfam(1,fam$linkinv(ij), fam))
  list(Z=z, Y=y)
}
get_Zhat <- function(Y, A){
   Zhat <- Y %*%  (A %*% solve( t(A) %*% A))
   matrix(Zq[(order(order(Zhat)))], nrow(Y), ncol(A))
}

get_Zhat_glm <- function(Y, A, family, maxit=3){
  t(apply(Y, 1, function(Yi){
    glm(Yi~0 + A, family=family, maxit=maxit)$coef
  }))
}

get_Ahat_glm <- function(Y, Z, family, maxit=3){
  t(apply(Y, 2, function(Yi) glm(Yi~0 + Z, family=family, control=list(maxit=maxit))$coef))
}


fam <- gaussian()
fam <- binomial()
fam <- poisson()

n <- 1000
p <- 200
q <- 5
A0 <- matrix(sample(c(-2,-1,0,1,2), size=p*q, repl=T), p, q)
A0 <- matrix(rnorm(p*q), p, q)
set.seed(414123)
dat <-genY(A0)
Y0 <- dat$Y
Z0 <- dat$Z
Zq <- qnorm((1:(n*q))/(1+n*q))

Zh <- Zq[order(order(Y0 %*% solve(A0 %*% t(A0) + diag(p))%*% A0))]

plot(Z0, Y0 %*% solve(fam$linkinv(A0 %*% t(A0))+diag(p))%*%A0);abline(0,1,col=2)

plot(Z0, get_Zhat_glm(Y0, A0, family = fam))



Ahat <- matrix(rnorm(p*q), p, q)*.5

par(mfrow=c(2,1))
for(i in 1:10){
Zhat <- get_Zhat_glm(Y0, Ahat, fam, maxit=3)
Ahat <- get_Ahat_glm(Y0, Zhat, fam, maxit=3)
Ahat <- psych::Procrustes(Ahat, A0)$loadings
plot(A0, Ahat);abline(0,1,col=2)
plot(Z0, Zhat);abline(0,1,col=2)
}
par(mfrow=c(1,1))


get_Z2 <- function(Y, A)log(Y+1) %*% A/norm(A)
get_A2 <- function(Y, Z, family){
  t(log(Y+1)) %*% Z /norm(Z)
}

mu <- fam$linkinv(Z0 %*% t(A0))

# Q <- solve(t(mu) %*% mu) %*% t(mu) %*% Z0
# plot(Z0, Y0 %*% Q)

Zhat <- scale(Y0, scale=F) %*% solve(exp(A0) %*% t(exp(A0)) + diag(p)) %*% exp(A0)
Zhat <- scale(Y0, scale=F) %*% solve(exp(A0 %*% t(A0)) + diag(p)) %*% A0
plot(Z0, Zhat)




A <- matrix(rnorm(p*q), p, q)*.1

A <- get_A2(Y0, get_Z2(Y0,A), family=family)
A <- psych::Procrustes(A, A0)$loadings
plot(A0, A); abline(0,1,col=2)


onestep <- function(A){
  # from observed sample
  Z.obs <- get_Zhat(Y0, A)
  A.obs <- get_Ahat(Y0, Z.obs, family = fam)
  # from generated sample
  Y.gen <- genY(A)$Y
  Z.gen <- get_Zhat(Y.gen, A)
  A.gen <- get_Ahat(Y.gen, Z.gen, family=fam)

  A.obs - A.gen
}

A <- matrix(rnorm(p*q), p, q)

Z <- get_Zhat(Y0, A)
A <- get_Ahat(Y0, Z, family)

A <- psych::Procrustes(A, A0)$loadings
plot(A0, A); abline(0,1,col=2)


for(i in 1:10){
  set.seed(123)
  A <- A + .5* onestep(A)
  A <- psych::Procrustes(A, A0)$loadings
  plot(A, A0)
}
plot(A, A0); abline(0,1,col=2)
