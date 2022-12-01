

est.gmf <- function(dat){
  tm <- proc.time()
  model.newton = gmf(Y = dat$Y, X = dat$X, p=dat$dimensions$q, maxIter = 100,
                     family=binomial(), tol = 1e-5, intercept = F, method="quasi")
  tmdiff <- proc.time()-tm

  A <- model.newton$v
  B <- t(model.newton$beta)

  list(A=A, B=B, time=tmdiff)
}

est.gllvm <- function(dat){
  tm <- proc.time()
  fit.gllvm <- gllvm(y= dat$Y, X= dat$X, formula = ~0 +., num.lv=dat$dimensions$q, family=binomial(link="logit"))
  tmdiff <- proc.time()-tm

  coefs <- coefficients(fit.gllvm)
  A <- coefs$theta
  B <- coefs$Xcoef

  list(A=A, B=B, time=tmdiff)
}

est.sprime <- function(dat){
  tm <- proc.time()
  fit.sa <- fastgllvm(dat$Y, dat$dimensions$q, family = "binomial", intercept = T, method="simple", alpha=.2, maxit=50, hist=100, batch_size = 500)
  tmdiff <- proc.time()-tm
  fit.sa <- update(fit.sa, alpha=.05)
  A <- fit.sa$parameters$A
  B <- fit.sa$parameters$B

  list(A=A, B=B, time=tmdiff)
}


est.ltm <- function(dat) {
  library(ltm)
  tm <- proc.time()
  if(q==1) fit.ltm <- ltm(dat$Y ~ z1)
  if(q==2) fit.ltm <- ltm(dat$Y ~ z1 + z2)
  tmdiff <- proc.time()-tm

  A <- as.matrix(fit.ltm$coefficients[,2:(1 + dat$dimensions$q), drop=F])
  B <- as.matrix(fit.ltm$coefficients[,1, drop=F])
  list(A=A, B=B, time=tmdiff)
}

est.prime <- function(dat){
  tm <- proc.time()
  fit.sa <- fastgllvm(dat$Y, dat$dimensions$q, family = "binomial", intercept = T, method="full", alpha=.2, maxit=50, hist=100, batch_size = 500)
  tmdiff <- proc.time()-tm
  fit.sa <- update(fit.sa, alpha=.05)
  A <- fit.sa$parameters$A
  B <- fit.sa$parameters$B

  list(A=A, B = B, time=tmdiff)
}


est.mirtjml <- function(dat) {
  library(mirtjml)
  tm <- proc.time()
  fit <- mirtjml::mirtjml_expr(dat$Y, dat$dimensions$q, tol=1e-3)
  tmdiff <- proc.time()-tm
  A <- fit$A_hat
  B <- fit$d_hat

  list(A=A, B = B, time=tmdiff)
}

est.all <- function(dat){
  list(prime=est.prime(dat),
       sprime=est.sprime(dat),
       gmf=est.gmf(dat),
       gllvm=est.gllvm(dat),
       ltm=est.ltm(dat),
       mirtjml = est.mirtjml(dat)
       )
}



gen_dat <- function(n, p, q, A, B){
  dat <- gen_fastgllvm(nobs = n, p=p, q=q, k=1, family="binomial", intercept = T, A = A, B=B)
  dat
}


onesim <- function(seed, n, p, A, B){
  library(gmf)
  library(gllvm)
  library(ltm)
  library(mirtjml)
  devtools::load_all()
  set.seed(seed)
  dat <- gen_dat(n, p, q=2, A=A, B=B)
  simres <- tryCatch(c(model=list(dat), est.all(dat), seed=seed, n=n, p=p), error=function(e) paste0("Error with seed=", seed, ", n=", n, ", p=",p, ", q=", 2 ))
  saveRDS(simres, file=paste0("experiments/GLLVM_binary/low_dimensions/simres/n", n, "_p", p, "_seed", seed,".rds"))
}



devtools::load_all()
library(ltm)
library(mirtjml)
library(gmf)
library(gllvm)

family <- binomial()

# Simulation settings:
rep <- 50
q <- 2
p.list <- c(20, 40)
n.list <- c(1000)


# We compare 5 methods:

settings <- expand.grid(n.list, p.list)
colnames(settings) <- c("n", "p")

library(parallel)
# Detect the number of available cores and create cluster
cl <- parallel::makeCluster(detectCores()-2)
# Run parallel computation
parallel::clusterExport(cl, ls())

for(i in 1:nrow(settings)){
  n <- as.numeric(settings[i,1])
  p <- as.numeric(settings[i,2])

  A <- matrix(runif(p*q,-2,2), p, q)
  B <- matrix(runif(p, -1,1), p,1)

  parallel::parLapply(cl, 1:50, onesim, n=n, p=p, A=A, B=B)
  # Close cluster
}
parallel::stopCluster(cl)
