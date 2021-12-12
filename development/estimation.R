devtools::load_all()

n <- 1000
p <- 100
q <- 1
dat <- gen_gllvm(n,family="bernoulli", p=p, q=q, k=10, scale=1)
hist(sigmoid(dat$natpar))

image(cov(dat$Y))

fit <- bernoulli.ffa.stochastic(dat$Y, q)
plot(fit$A, dat$par$A)


fit <- bernoulli.ffa.covariate(dat$Y, q, dat$X)

plot(fit$A, dat$par$A)
