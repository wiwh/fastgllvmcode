devtools::load_all()
# missing data completely at random


n <- 10000
p <- 20
beta <- matrix(rnorm(p))
beta.miss <- matrix(rnorm(p*p), p, p)

fam <- gaussian()

X <- matrix(rnorm(n*p), n, p)
X[,1] <- 1
linpar <- X %*% beta
Y <- rnorm(n, linpar, 1)
M <- matrix(runif(n*p) < prob.m, n, p)*1
eta.miss <- X %*% beta.miss
cov <- rWishart(1, p, diag(p))[,,1]/p/5
eta.miss <- eta.miss %*% t(chol(cov))

M <- apply(eta.miss, 2, function(eta_j)rbinom(n, 1, sigmoid(eta_j)))



fit <- fastgllvm(M, q=3, X=X, learning_rate.args = list(start=10, end=10))
plot_fastgllvm(fit)
plot(cov(M), cov(simulate_fastgllvm(fit, n=1000)$Y))

plot(fit$B, t(beta.miss))


