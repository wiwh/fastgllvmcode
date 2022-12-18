est.gmf <- function(dat){
  tm <- proc.time()
  model.newton = gmf(Y = dat$Y, X = dat$X, p=dat$dimensions$q, maxIter = 200,
                     family=poisson(), intercept=F, method="quasi", tol=1e-5)
  tmdiff <- proc.time()-tm

  A <- model.newton$v
  B <- t(model.newton$beta)
  Z <- model.newton$u

  list(A=A, B=B, Z=Z, time=tmdiff)
}

est.gllvm <- function(dat){
  tm <- proc.time()
  fit.gllvm <- gllvm(y= dat$Y, X= dat$X, formula = ~0 +., num.lv=dat$dimensions$q, family=poisson(), method = "EVA", sd.errors=F)
  # fit.gllvm <- gllvm(y = dat$Y, family="poisson", method="VA")
  tmdiff <- proc.time()-tm

  coefs <- coefficients(fit.gllvm)
  A <- coefs$theta
  B <- coefs$Xcoef
  Z <- as.matrix(fit.gllvm$lvs)
  rownames(Z) <- NULL
  colnames(Z) <- NULL

  list(A=A, B=B, Z=Z, time=tmdiff)
}

est.sprime <- function(dat){
  tm <- proc.time()
  fit.sa <- fastgllvm(dat$Y, dat$dimensions$q, family = "poisson", intercept = T, method="simple", alpha=.15, maxit=100)
  tmdiff <- proc.time()-tm
  fit.sa <- update(fit.sa, alpha=.05, maxit=50)
  A <- fit.sa$parameters$A
  B <- fit.sa$parameters$B
  Z <- fit.sa$Z

  list(A=A, B=B, Z=Z, time=tmdiff)
}


est.prime <- function(dat){
  tm <- proc.time()
  fit.sa <- fastgllvm(dat$Y, dat$dimensions$q, family = "poisson", intercept = T, method="full", alpha=.15, maxit=100)
  tmdiff <- proc.time()-tm
  fit.sa <- update(fit.sa, alpha=.05)
  A <- fit.sa$parameters$A
  B <- fit.sa$parameters$B
  Z <- fit.sa$Z

  list(A=A, B=B, Z=Z, time=tmdiff)
}

est.all <- function(dat){
  list(prime=est.prime(dat),
       sprime=est.sprime(dat),
       # gllvm = est.gllvm(dat),
       gmf=est.gmf(dat)
  )
}


# gen A according to the two sparsity settings (see file "loading_structure.R" in the same directory as this current one)

gen_A <- function(p, q, setting="A", prop=.4) {
  # setting "A" has only the 100 first loadings that are non-zero
  A <- matrix(runif(p*q, -2, 2), p, q)
  if (setting == "A") {
    if(p>100) {
      A[101:p, ] <- 0
    }
  } else if (setting =="B") {
    for (k in 1:q) {
      nonzeros <- (1:p) %in% ((1:(p*prop)) + round((k-1) * p * (1-prop)/4 ))
      A[!nonzeros,k] <- 0
    }
  }
  A
}

onesim <- function(seed, n, p, A, B, setting){
  library(gllvm)
  library(gmf)
  devtools::load_all()
  set.seed(seed)
  dat <-  gen_fastgllvm(nobs = n, p=p, q=q, k=1, family="poisson", intercept = T, A = A, B=B)
  simres <- tryCatch(c(model=list(dat), est.all(dat), seed=seed, n=n, p=p, setting=setting), error=function(e) paste0("Error with seed=", seed, ", n=", n, ", p=",p, ", q=", 2 ))
  saveRDS(simres, file=paste0("experiments/GLLVM_poisson/simres/n", n, "_p", p, "_seed", seed, "_setting", setting,".rds"))
}


devtools::load_all()
library(gmf)
library(gllvm)

family <- binomial()

# Simulation settings:
q <- 5
p.list <- c(100, 200, 200, 300, 400, 500)
n.list <- c(100, 500)
setting <- c("A", "B")


# We compare 5 methods:

settings <- expand.grid(n.list, p.list, setting)
colnames(settings) <- c("n", "p", "setting")
settings$setting <- as.character(settings$setting)

library(parallel)
# Detect the number of available cores and create cluster
cl <- parallel::makeCluster(detectCores()-2)
# Run parallel computation
parallel::clusterExport(cl, ls())

for(i in 1:nrow(settings)){
  n <- as.numeric(settings[i,1])
  p <- as.numeric(settings[i,2])
  setting <- settings[i,3]

  set.seed(123132 + i)
  B <- matrix(runif(p, -2,2), p,1)
  A <- gen_A(p, 5, setting=setting)

  parallel::parLapply(cl, 1:50, onesim, n=n, p=p, A=A, B=B, setting)
  # Close cluster
}
parallel::stopCluster(cl)


# example
if(0){
  devtools::load_all()
  poisson  <- 100
  gaussian <- 0
  binomial <- 100
  nobs <- 500
  q <- 5
  p <- poisson + gaussian + binomial

  intercept <- T
  k <- 1
  if(k==0 & intercept) k <- 1


  family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial))
  set.seed(14240)
  A <- gen_A(p, 5, "A", prop=.4)
  B <- matrix(runif(p*k, -2, 2), p, k)
  fg <- dat <- gen_fastgllvm(nobs=nobs, p=p, q=q, k=k, family=family, intercept=intercept, A=A, B=B)

  fit.fg <- fastgllvm(fg$Y, q=q, family="poisson", intercept=T, alpha=.2, maxit = 50)
  plot(fit.fg)

  library(mirtjml)
  fit.m <- mirtjml_expr(fg$Y, K=q, tol = 2)

  library(gmf)
  fit.gmf <- gmf(fg$Y, p=q, intercept=T, family=binomial(), penaltyU = .1)

  MPE(fit.fg$parameters$A, fg$parameters$A)
  MPE(fit.m$A_hat, fg$parameters$A)
  MPE(fit.gmf$v, fg$parameters$A)

  plot(fg$parameters$A, psych::Procrustes(fit.fg$parameters$A, fg$parameters$A)$loadings)
  points(fg$parameters$A, psych::Procrustes(fit.m$A_hat, fg$parameters$A)$loadings, col=2)
  points(fg$parameters$A, psych::Procrustes(fit.gmf$v, fg$parameters$A)$loadings, col=3)
  abline(0,1,col=2)
}
