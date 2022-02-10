devtools::load_all()
family <- binomial()
family <- poisson()
library(gmf)
library(glmnet)

n <- 10000

q <- 1
p <- 10
k <- 0
set.seed(1241)
fg <- gen_fastgllvm(n = n, p=p, q=q, k=k,family =family, intercept = T, A=matrix(rnorm(p*q), p, q)) #, B=matrix(sample(-1:5, p, repl=T)))


A <- matrix(rnorm(p*q), p, q)

natpar <- compute_natpar(A, fg$B, fg$Z, fg$X)
mus <- family$linkinv(natpar)
EE <- cov(mus) + diag(colMeans(family$variance(mus)))
d <- sqrt(diag(EE))
c1 <- t(t(EE / d)/d)
YY <- cov(fg$Y)
d <- sqrt(diag(YY))
c2 <- t(t(YY / d) / d)
A <- A + .1*(c1 - c2)%*% A
plot(fg$A, A)


# here we test the Newton Method
np <- family$linkinv(fg$natpar)**2
w <- family$mu.eta(fg$natpar)*family$linkinv(fg$natpar) + family$linkinv(fg$natpar)
# boxplot((t(t(np)/colMeans(w))), outline=F)

plot(fg$Z, t(t(scale(fg$Y, scale=F))/colMeans(w))%*% fg$A/p)
plot(fg$Z, t(t(scale(fg$Y, scale=F))/colMeans(w))%*% fg$A/p)


(t(np)/colMeans(w)) %*% fg$Z / n

fit <- fit_fastgllvm(fg, learning_rate.args = list(start=1e-2, end=1e-4), maxit = 50)
plot_fastgllvm(fit)

# true prediction and wrong predictions
if(0){
  Z <- predict.fastgllvm(fg)[,1]
  Z1 <- (fg$Y - fg$X %*% t(fg$B)) %*% fg$A[,1,drop=F]
  Zq <- qnorm((1:n)/(n+1))
  Zhat <- Zq[order(order(Z1))]
  plot(fg$Z[,1], Zhat)
  points(fg$Z[,1], Z, col=2)
}


# whole fit using the juju's trick
Z <- matrix(rnorm(n*q), n, q)
for(i in 1:10){
  if(ncol(Z)==1){
    A <-  matrix(sapply(1:p, function(j)glm(fg$Y[,j] ~0 + Z, family=family)$coef))
  } else {
    A <-  t(sapply(1:p, function(j)glm(fg$Y[,j] ~0 + Z, family=family)$coef))
  }
  Z <- sapply(1:q, function(j){
    Z <- fg$Y %*% A[,j, drop=F]
    Zq[order(order(Z))]
  })
  plot(fg$A, psych::Procrustes(A, fg$A)$load, main=i)
  # plot(fg$Z, Z, main=i)
}
A <- psych::Procrustes(A, fg$A)$load
plot(fg$A, A)
#re-compute after procrustes!
plot(fg$Z, Z)

fit <- fastgllvm(fg$Y, q=q, fg$X, method="SP", H=10, family="poisson",
                  learning_rate = "exp", maxit=100,
                 learning_rate.args = list(start=1e-1, end=1e-1), tol=0)
plot_fastgllvm(fit)


ffa <- fastfactoranalysis::ffa(scale(fg$Y),nfactors = 1)

f1 <- psych::fa(fg$Y, nfactors = 1)

f2 <- gmf(fg$Y, X=fg$X, p=1, family= poisson())


fit <- fit_fastgllvm(fit, learning_rate.args = list(start=1e-1, end=1e-5))
Z <- predict.fastgllvm(fit)
Zq <- qnorm((1:n)/(n+1))
Zhat <- Zq[order(Z)]



plot(fg$A, psych::Procrustes(fit$A, fg$A)$loadings)
abline(0,1,col=2)
plot(fg$Z, Z)
# Zhat <- predict(fg, method="ridge", lambda=.001)
# plot(-fg$Z, Zhat); abline(0,1,col=2)

n <- fg$n
H <- 1
A <- fg$A
B <- fg$B
phi <- fg$phi
q <- fg$q
p <- fg$p
method <- "SA"

H <- 10
generate_Z <- generate_Z_functionfactory(n, q, H, method)
Y <- fg$Y
Y.c <- scale(Y, scale=F)
X <- fg$X
family <- fg$family

# sims <- lapply(1:1000, function(na) get_Psi(Y, Y.c, A, B, phi, X, family, generate_Z))
# A.sim <- t(sapply(sims, function(simi)simi$A))
# boxplot(A.sim)



K <- get_K(A, phi)
E <- get_expectations(A, B, X, generate_Z(), K, phi, family)
plot(t(Y.c) %*% (Y.c %*% t(K))/nrow(Y.c), E$EYYK); abline(0,1,col=2)


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
