devtools::load_all()
n <- 10000
p <- 100
q <- 1

family <- binomial()

Z <- matrix(rnorm(n*q), n, q)
A <- matrix(rnorm(p*q), p, q)
B <- matrix(1, p, 1)
X <- matrix(1, n, q)
phi <- rep(1,p)
linpar <- compute_natpar(A, B, Z, X)
K <- get_K(A, phi)


mu <- family$linkinv(linpar)
Y <- sapply(1:p, function(j){
  if(family$family=="binomial")Yj <- rbinom(n, 1, family$linkinv(mu[,j]))
  if(family$family=="poisson")Yj <- rpois(n, family$linkinv(mu[,j]))
  Yj
})

Y.c <- scale(Y, scale=F)
mu.c <- scale(mu, scale=F)
# all.equal(family$variance(mu) , sapply(1:p, function(j) family$variance(mu[,j])))
var <- colMeans(family$variance(mu))

YY <- t(Y.c) %*% Y.c/nrow(Y.c)
mumu <- t(mu.c) %*% mu.c/nrow(mu.c) + diag(var)
plot(YY, mumu); abline(0,1,col=2)



#slow
mumuK <- (t(mu.c) %*% (mu.c)/nrow(mu.c) + diag(var)) %*% t(K)
YYK <- (t(Y.c) %*% Y.c /nrow(Y.c)) %*% t(K)
# fast
YYK2 <- t(Y.c) %*% (Y.c %*% t(K)) / nrow(Y.c)
mumuK2 <- t(mu.c) %*% (mu.c %*% t(K))/nrow(mu.c) + t(K) * var

E <- get_expectations(A, B, X, list(Z), K, phi, family=family)

mumuK3 <- E$EYYK

all.equal(mumuK, mumuK2, mumuK3)

plot(YYK, mumuK); abline(0,1,col=2)
