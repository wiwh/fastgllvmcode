sim<- readRDS("experiments/GLLVM_binary/large_dimensions/simres/n500_p700_seed6_settingA.rds")

MPE(sim$prime$A, sim$model$parameters$A)
MPE(sim$sprime$A, sim$model$parameters$A)
MPE(sim$gmf$A, sim$model$parameters$A)
MPE(sim$mirtjml$A, sim$model$parameters$A)

compute_mu <- function(Z, A, X, B){
  1/(1+exp(-(Z %*% t(A) + X %*% t(B))))
}

mu <- compute_mu(sim$prime$Z, sim$prime$A, sim$model$X, sim$prime$B)
MD(sim$model$Y, mu, family="binomial")
mu <- compute_mu(sim$sprime$Z, sim$sprime$A, sim$model$X, sim$sprime$B)
MD(sim$model$Y, mu, family="binomial")
mu <- compute_mu(sim$mirtjml$Z, sim$mirtjml$A, sim$model$X, sim$mirtjml$B)
MD(sim$model$Y, mu, family="binomial")
mu <- compute_mu(sim$gmf$Z, sim$gmf$A, sim$model$X, sim$gmf$B)
MD(sim$model$Y, mu, family="binomial")


plot(sim$model$parameters$A, psych::Procrustes(sim$prime$A, sim$model$parameters$A)$loadings)
points(sim$model$parameters$A, psych::Procrustes(sim$sprime$A, sim$model$parameters$A)$loadings, col=2)
points(sim$model$parameters$A, psych::Procrustes(sim$gmf$A, sim$model$parameters$A)$loadings, col=3)
points(sim$model$parameters$A, psych::Procrustes(sim$mirtjml$A, sim$model$parameters$A)$loadings, col=4)
abline(0,1,col=2)

sim$mirtjml$time
sim$sprime$time
