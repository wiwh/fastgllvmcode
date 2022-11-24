# This reproduces the experimental setup of 3.1 dai.2020.matrix-free-method-high-dimensional

# For each setting, we compute the mean procrustes error.

compute_FFA_error <- function(Y, A, Psi) {
  K <- ffa_comp_K(A, Psi)
  SK <- t(Y) %*% (Y %*% K)/nrow(Y)
  EK <- A %*% (t(A) %*% K) + K * Psi
  norm(EK-SK, type = "F")/norm(EK, type="F")
}
#times

# settings : (n, p, q, r): we test from 1 to 10 p
settings <- list(
  c(100, 2000, 5),
  c(200, 4000, 5),
  c(400, 8000, 5)
)

settings <- lapply(settings, function(setting) {
  t(sapply(1:10, function(r) c(setting, r)))
})

settings <- do.call(rbind, settings)
colnames(settings) <- c("n", "p", "q", "r")

# we run 50 of each

sims <- lapply(1:nrow(settings), function(j) {
  setting <- settings[j,]

  n <- setting[1]
  p <- setting[2]
  q <- setting[3]
  r <- setting[4]

  lapply(1:20, function(i) {
    cat("\nSetting", j, ", (n, p, q, r) = ", paste(settings[j,], collapse=","), " iter ", i, "/20.")
    set.seed(1234 + 20*j + i)
    fg <- gen_fastgllvm(nobs=n, p=p, q=q, k=1, family="gaussian", intercept=T, phi=runif(p, .2, .8))
    fg$Y <- scale(fg$Y, scale=F)

    fg$Y <- scale(fg$Y, scale=F)

    fit.fad <- fad::fad(fg$Y, factors = q)
    fit.ffa <- ffa(fg$Y, q=q, maxiter = 100, eps=1e-10, savepath = T, rotate_updates = F)

    sds <- sqrt(apply(fg$Y, 2, var) * ((n-1)/n))

    fad.time <- system.time({fit.fad <- fad::fad(fg$Y, factors = r)})
    ffa.time <- system.time({fit.ffa <- ffa(fg$Y, q=r, maxiter = 100, eps=1e-10, savepath = F, rotate_updates = F)})

    fad.error.equality <- compute_FFA_error(fg$Y, fit.fad$loadings * sds, fit.fad$uniquenesses * sds**2)
    ffa.error.equality <- compute_FFA_error(fg$Y, fit.ffa$A, fit.ffa$Psi)

    list(n=n, p=p, q=q, r=r, seed=1234 + 20*j + i, fad.time=fad.time, ffa.time=ffa.time, fad.error.equality=fad.error.equality, ffa.error.equality = ffa.error.equality)
  })
})
saveRDS(sims, file="./experiments/FFA/sims_speed_fad_5.rds")


if(0) {


K <- ffa_comp_K(fg$parameters$A, fg$parameters$phi)
Sig <- comp_Sig(fg$parameters$A, fg$parameters$phi)

all.equal(fg$parameters$A, Sig %*% K)

comp_Sig <- function(A, phi) {
  A %*% t(A) + diag(phi)
}

n <- 400
p <- 8000
q <- 5


fg <- gen_fastgllvm(nobs=n, p=p, q=q, k=1, family="gaussian", intercept=T, phi=rep(1, p))

# fit50 <- ffa(fg$Y, q, savepath = T, verbose=T, maxiter=50, eps=-1)
# fit500 <- ffa(fg$Y, q, savepath = T, verbose=T, maxiter=500, eps=-1)

# all.equal(fit50$A, psych::Procrustes(fit500$A, fit50$A)$loadings)

# plot(fit50$A, fit500$A)

fg$Y <- scale(fg$Y, scale=F)

q <- 5

fit.fad <- fad::fad(fg$Y, factors = q)
fit.ffa <- ffa(fg$Y, q=5, maxiter = 100, eps=1e-10, savepath = F, rotate_updates = F)
print(fit.ffa$niter)


ts.plot(fit.ffa$path$A)
ts.plot(fit.ffa$path$Psi)




mgm <- microbenchmark::microbenchmark(fad::fad(fg$Y, factors = q),
                               ffa(fg$Y, q=q, maxiter = 200, eps=1e-10, savepath = T, rotate_updates = F), times=20)
plot(mgm)

sds <- sqrt(apply(fg$Y, 2, var) * (n-1)/n)

compute_FFA_error(fg$Y, fit.fad$loadings*sds, fit.fad$uniquenesses*(sds**2))
compute_FFA_error(fg$Y, fit.ffa$A, fit.ffa$Psi)


plot(cov(fg$Y), comp_Sig(fit.fad$loadings*sds, fit.fad$uniquenesses*(sds**2)))
points(cov(fg$Y), comp_Sig(fit.ffa$A, fit.ffa$Psi), col=2)



a <- system.time({fgtest <- ffa(fg$Y, 10, maxiter = 100, eps=1e-10, savepath = T)})
b <- system.time({fgtest <- ffa(fg$Y, 10, maxiter = 100, eps=1e-10, savepath = T)})

fg$Y <- scale(fg$Y, scale=F)

sds <- apply(fg$Y, 2, sd)
tic()
fit2 <- ffa(fg$Y, 10, maxiter = 100, eps=1e-10, savepath = T)
b <- toc()
ts.plot(fit2$path$Psi)


}
