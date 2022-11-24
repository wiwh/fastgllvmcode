devtools::load_all()
library(ltm)
library(mirtjml)
library(gmf)
library(gllvm)

family <- binomial()

# Simulation settings:
rep <- 1000
q <- 2
p.list <- c(20, 40)
n.list <- c(50, 100, 500)

# We compare 4 methods:

est.gmf <- function(dat){
  tm <- proc.time()
  model.newton = gmf(Y = dat$Y, X = dat$X, d = dat$dimentions$k, p=dat$dimensions$q, maxIter = 100,
                     family=family, tol = 1e-5, intercept = F)
  tmdiff <- proc.time()-tm

  A.est <- model.newton$v

  proc <- MPE(A.est, dat$parameters$A)
  list(A=proc$Lrot, error=proc$Error, time=tmdiff)
}

est.gllvm <- function(dat){
  tm <- proc.time()
  fit.gllvm <- gllvm(y= dat$Y, X= dat$X, formula = ~., num.lv=dat$dimensions$q, family=binomial(link="logit"))
  tmdiff <- proc.time()-tm

  coefs <- coefficients(fit.gllvm)
  A.est <- coefs$theta

  proc <- procrustes(A.est, dat$par$A)
  list(A=proc$Lrot, error=proc$Error, time=tmdiff)
}

est.sprime <- function(dat){
  tm <- proc.time()
  fit.sa <- fastgllvm(dat)
  tmdiff <- proc.time()-tm

  A.est <- fit.sa$A

  proc <- procrustes(A.est, dat$par$A)
  list(A=proc$Lrot, error=proc$Error, time=tmdiff)
}

est.prime <- function(dat){
  tm <- proc.time()
  fit.sa <- fastgllvm(dat$Y, dat$dimensions$q, family = "binomial", intercept = T, method="full", alpha=1, maxit=100, hist=100)
  tmdiff <- proc.time()-tm
  for(i in 1:4) fit.sa <- update(fit.sa, fit.sa$controls$a*10, H=10, maxit=10)

  plot(fit.sa)

  A.est <- fit.sa$parameters$A

  MPE(A.est, dat$parameters$A)

  proc <- psych::Procrustes(A.est, dat$parameters$A)
  list(A=proc$Lrot, error=proc$Error, time=tmdiff)
}

est.all <- function(dat){
  list(sp=est.sp(dat),
       sa=est.sa(dat),
       gmf=est.gmf(dat),
       gllvm=est.gllvm(dat))
}


gen_dat <- function(n, p, q){
  A <- matrix(0, p, q)
  A[,1] <- seq(2,-2, l=p)
  if(q==2) {
    A[,2] <- seq(-1,1, l=p)
  }
  dat <- gen_fastgllvm(nobs = n, p=p, q=q, k=1, family="binomial", intercept = T, A = A)
  dat
}

onesim <- function(seed, n, p){
  set.seed(seed)
  library(gmf)
  library(gllvm)
  library(ltm)
  library(mirtjml)
  devtools::load_all()
  dat <- gen_dat(n, p)
  c(dat=dat, est.all(dat), seed=seed, n=n, p=p)
}

settings <- expand.grid(n.list, p.list)
colnames(settings) <- c("n", "p")

library(parallel)
# Detect the number of available cores and create cluster
cl <- parallel::makeCluster(detectCores()-2)
# Run parallel computation
parallel::clusterExport(cl, ls())

bigsim <- vector(mode="list", l=nrow(settings))

for(i in seq_along(bigsim)){
  n <- as.numeric(settings[i,1])
  p <- as.numeric(settings[i,2])
  bigsim[[i]] <- parallel::parLapply(cl,
                                     1:500,
                                     onesim, n=n, p=p)
  # Close cluster
}
save(bigsim, file="./experiments/low_dimensions/low_dimensions_results.Rdata")
parallel::stopCluster(cl)
