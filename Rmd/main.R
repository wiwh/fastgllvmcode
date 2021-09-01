# example
n <- 1000
p <- 10
q <- 3

par <- gen_par(p, q, family="bernoulli")

dat <- gen_gllvm(n, p, q, family="bernoulli", par=par)

# conditional distribution
