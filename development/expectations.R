devtools::load_all()
n <- 1000
p <- 250
q <- 2

family <- poisson()

Z <- matrix(rnorm(n*q), n, q)
A <- matrix(rnorm(p*q), p, q) * ifelse(family$family == "poisson", 2, 1)
B <- matrix(runif(p, -1, 1), p, 1)
X <- matrix(1, n, 1)
phi <- rep(1,p)
linpar <- compute_natpar(A, B, Z, X)
K <- get_K(A, phi)


mu <- family$linkinv(linpar)
Y <- sapply(1:p, function(j){
  if(family$family=="binomial")Yj <- rbinom(n, 1,mu[,j])
  if(family$family=="poisson")Yj <- rpois(n, mu[,j])
  Yj
})

Y.c <- scale(Y, scale=F)
mu.c <- scale(mu, scale=F)
# all.equal(family$variance(mu) , sapply(1:p, function(j) family$variance(mu[,j])))
var <- colMeans(family$variance(mu))



#slow
YY <- t(Y.c) %*% Y.c/nrow(Y.c)
mumu <- t(mu.c) %*% mu.c/nrow(mu.c) + diag(var)
YYK <- YY %*% t(K)
mumuK <- mumu %*% t(K)
# fast
YYK2 <- t(Y.c) %*% (Y.c %*% t(K)) / nrow(Y.c)
mumuK2 <- t(mu.c) %*% (mu.c %*% t(K))/nrow(mu.c) + t(K) * var

get_mumuK_fast <- function(){
  mu <- family$linkinv(compute_natpar(A, B, Z, X))
  K  <- get_K(A, phi)
  mu.c <- scale(mu, scale=F)
  var <- colMeans(family$variance(mu))
  t(mu.c) %*% (mu.c %*% t(K))/nrow(mu.c) + t(K) * var
}
mumuK3 <- get_mumuK_fast()  # this is same as 2 but without pre-computing quantities

E <- get_expectations(A, B, X, list(Z), K, phi, family=family)
mumuK4 <- E$EYYK

all.equal(mumuK, mumuK2) && all.equal(mumuK2, mumuK3) && all.equal(mumuK3, mumuK4)

microbenchmark::microbenchmark(get_mumuK_fast(), get_expectations(A, B, X, list(Z), K, phi, family=family), times=10)

par(mfrow=c(2,1))
plot(YY, mumu, main=all.equal(YY, mumu)); abline(0,1,col=2)
plot(YYK, mumuK, main=all.equal(YY, mumu)); abline(0,1,col=2)
par(mfrow=c(1,1))
