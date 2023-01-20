devtools::load_all()
set.seed(1234)
poisson  <- 0
gaussian <- 10
binomial <- 0
q <- 1
intercept <- F
k <- ifelse(intercept, 1, 0)

p <- poisson + gaussian + binomial
family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial))
set.seed(102)
fg <- gen_gllvmprime(nobs=1000, p=p, q=q, family=family, phi=3*(1:p)/p, k=k, A=matrix((seq(-2, 2, l=p*q)),p, q), intercept=intercept, miss.prob = 0, scale=1)

# Some testing when the parameters are not correct (to simulate what happens in training)
#----------------------
weight_A <- 1.2
weight_B <- 1.
bias_A   <- -0.4
bias_B   <- 0.


# set.seed(1231)
model_true <- simulate(fg, return_object = T)
model_true <- compute_mean(model_true, return_object = T)
# pretend the parameters are wrong
model_wrong <- model_true
model_wrong$parameters$A <- model_wrong$parameters$A * weight_A + bias_A
if(!is.null(model_wrong$parameters$B)) {
  model_wrong$parameters$B <- model_wrong$parameters$B * weight_B + bias_B
}

# Compare the mean computed under the true and wrong model
model_true <- compute_Z(model_true, return_object = T)
model_true <- compute_mean(model_true, return_object = T)

model_wrong <- compute_Z(model_wrong, return_object = T)
model_wrong <- compute_mean(model_wrong, return_object = T)


# The two Z are, of course, not the same...
plot(model_true$Z, model_wrong$Z); abline(0,1,col=2)

# ... but the two means are almost the same!
plot(model_true$mean, model_wrong$mean); abline(0,1,col=2)

# What about the dot product? It is *not* the same! And herein lies our information!
plot(t(model_true$mean) %*% model_true$Z, t(model_wrong$mean) %*% model_wrong$Z); abline(0,1,col=2)

# Information in  the centered means (biased gradient of PRIME)
plot(t(model_true$Y - model_true$mean) %*% model_true$Z, t(model_wrong$Y-model_wrong$mean) %*% model_wrong$Z); abline(0,1,col=2)
points(t(model_true$Y - model_true$mean) %*% model_true$X, t(model_wrong$Y-model_wrong$mean) %*% model_wrong$X, col=5)

# Information in the simplified scores (biased gradient of S-PRIME)
plot(t(model_true$Y) %*% model_true$Z, t(model_wrong$Y) %*% model_wrong$Z); abline(0,1,col=2)
points(t(model_true$Y) %*% model_true$X, t(model_wrong$Y) %*% model_wrong$X, col=5)

# Centering
# ---------
# To center, we simulate at the wrong model
# set.seed(1231123)
model_wrong_sim <- simulate(model_wrong, return_object = T)
model_wrong_sim <- compute_Z(model_wrong_sim, return_object = T)
model_wrong_sim <- compute_mean(model_wrong_sim, return_object = T)


par(mfrow=c(1,2))
# The information is not clear in the centered gradient of PRIME
plot(cbind(model_true$parameters$A, model_true$parameters$B), gradient_full(model_wrong)$AB -  gradient_full(model_wrong_sim)$AB, col=5); abline(h=0,col=2)
points(cbind(model_true$parameters$A), (gradient_full(model_wrong)$AB- gradient_full(model_wrong_sim)$AB)[,1:q], col=1)
# What about the S-PRIME? It is extremely clear in the centered gradient of S-PRIME!
plot(cbind(model_true$parameters$A, model_true$parameters$B), gradient_simple(model_wrong)$AB -  gradient_simple(model_wrong_sim)$AB, col=5); abline(h=0,col=2)
points(cbind(model_true$parameters$A), (gradient_simple(model_wrong)$AB- gradient_simple(model_wrong_sim)$AB)[,1:q], col=1)
par(mfrow=c(1,1))


# Now simulate to get the expected gradient at the given values of parameters!
# ----------------------------------------------------------------------------

simulate_centered_gradients <- function(gradient, weight_A=1, weight_B=1, bias_A=0, bias_B=0, rescale=F) {
  sim <- sapply(1:100, function(i){
    set.seed(i+3)
    if(i%%10==0 || i==1) cat("\n", i)
    model_true <- simulate(fg, return_object = T)
    model_true <- compute_Z(model_true, return_object = T)
    model_true <- compute_mean(model_true, return_object = T)

    # Pretend the parameters are not correct (to simulate what happens in training)
    #----------------------
    model_wrong <- model_true
    model_wrong$parameters$A <- model_wrong$parameters$A * weight_A + bias_A
    if (!is.null(model_wrong$parameters$B)) {
      model_wrong$parameters$B <- model_wrong$parameters$B * weight_B + bias_B
    }

    # Compare the mean computed under the true and wrong model
    model_true <- compute_Z(model_true, start = model_true$Z, return_object = T)
    model_true <- compute_mean(model_true, return_object = T)

    model_wrong <- compute_Z(model_wrong, start = model_true$Z, return_object = T)
    model_wrong <- compute_mean(model_wrong, return_object = T)

    model_wrong_sim <- simulate(model_wrong, return_object = T)
    model_wrong_sim <- compute_Z(model_wrong_sim, start = model_wrong_sim$Z, return_object = T)
    model_wrong_sim <- compute_mean(model_wrong_sim, return_object = T)

    gradient_wrong <- gradient(model_wrong)
    gradient_wrong_sim <- gradient(model_wrong_sim)


    gradient_diff <- sapply(seq_along(gradient_wrong), function(i) {
      gradient_wrong[[i]] - gradient_wrong_sim[[i]]
    }, simplify=F)
    names(gradient_diff) <- names(gradient_wrong)

    if (rescale) {
      hessian <- simulate_hessian_AB(model_wrong, H=10)
      gradient_diff$AB <- mult_invHessian_dAB(gradient_diff$AB, hessian)
    }

    names(gradient_diff) <- names(gradient_wrong)
    gradient_diff
  }, simplify=F)


  sim <- sapply(names(sim[[1]]),
                function(parname) do.call(rbind, sapply(sim, function(simi) as.vector(simi[[parname]]),
                                                        simplify=F)), simplify=F)
  sim
}


plot_sim <- function(sim) {
  if(is.null(fg$parameters$B)) {
    boxplot(sim$AB[,(order(fg$parameters$A))], outline=F); abline(h=0, col=2)
    points(colMeans(sim$AB[,(order(fg$parameters$A))]), col=2)
  } else {
    par(mfrow=c(2,1))
    boxplot(sim$AB[,(order(fg$parameters$A))], outline=F); abline(h=0, col=2)
    points(colMeans(sim$AB[,(order(fg$parameters$A))]), col=2)
    boxplot(sim$AB[,-(order(fg$parameters$A)), drop=F][,order(fg$parameters$B)], outline=F); abline(h=0, col=2)
    points(colMeans(sim$AB[,-(order(fg$parameters$A)), drop=F][,order(fg$parameters$B)]), col=2)
    par(mfrow=c(1,1))
  }
}

# Simulation setting
devtools::load_all()
set.seed(1234)
poisson  <- 10
gaussian <- 0
binomial <- 10
q <- 1
intercept <- F
k <- ifelse(intercept, 1, 0)

p <- poisson + gaussian + binomial
family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial))
set.seed(102)
fg <- gen_gllvmprime(nobs=1000, p=p, q=q, family=family, phi=rep(1,p), k=k, A=matrix((seq(-2, 2, l=p*q)),p, q), intercept=intercept, miss.prob = 0, scale=1)

# Some testing when the parameters are not correct (to simulate what happens in training)
#----------------------

weight_A <- 1
weight_B <- 1
bias_A   <- .5
bias_B   <- 0

sim_full <- simulate_centered_gradients(
  gradient_full,
  weight_A = weight_A,
  weight_B = weight_B,
  bias_A = bias_A,
  bias_B = bias_B,
  rescale= F
)
sim_full_rescaled <- simulate_centered_gradients(
  gradient_full,
  weight_A = weight_A,
  weight_B = weight_B,
  bias_A = bias_A,
  bias_B = bias_B,
  rescale= T
)
sim_simple <- simulate_centered_gradients(
  gradient_simple,
  weight_A = weight_A,
  weight_B = weight_B,
  bias_A = bias_A,
  bias_B = bias_B,
  rescale = F
)
sim_simple_rescaled <- simulate_centered_gradients(
  gradient_simple,
  weight_A = weight_A,
  weight_B = weight_B,
  bias_A = bias_A,
  bias_B = bias_B,
  rescale = T
)

plot_sim(sim_full)
plot_sim(sim_full_rescaled)
plot_sim(sim_simple)
plot_sim(sim_simple_rescaled)

# saving
sims_full <- list(
  sim_full = sim_full,
  sim_full_rescaled = sim_full_rescaled
)

sims_simple <- list(
  sim_simple = sim_simple,
  sim_simple_rescaled = sim_simple_rescaled
)

saveRDS(sims_full, file="./experiments/gradients/sims_gradient_full.rds")
saveRDS(sims_simple, file="./experiments/gradients/sims_gradient_simple.rds")


# Draw 2 graphs, 1 for
