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
  c(100, 2000, 3),
  c(200, 4000, 3),
  c(400, 8000, 3),
  c(100, 2000, 5),
  c(200, 4000, 5),
  c(400, 8000, 5)
)

settings <- lapply(settings, function(setting) {
  t(sapply(1:10, function(r) c(setting, r)))
})

settings <- do.call(rbind, s)
colnames(settings) <- c("n", "p", "q", "r")

# we run 50 of each

lapply(1:nrow(settings), function(j) {
  cat("\nSetting", j, ", (n, p, q, r) = ", paste(settings[j,], collapse=","))
  setting <- settings[j,]

  n <- setting[1]
  p <- setting[2]
  q <- setting[3]
  r <- setting[4]

  lapply(1:50, function(i) {
    set.seed(1234 + 50*j + i)
    fg <- gen_fastgllvm(nobs=nobs, p=p, q=q, k=k, family=family, intercept=intercept, miss.prob = 0, scale=1, phi=rep(1, p))
    fg$Y <- scale(fg$Y, scale=F)
    sds <- apply(fg$Y, 2, sd)

    fad.time <- system.time({fit.fad <- fad::fad(fg$Y, factors = r)})
    ffa.time <- system.time({fit.ffa <- ffa(fg$Y, q=r, maxiter = 20, eps=1e-5, savepath = T)})

    fad.error <- compute_error(fit1$loadings  * sds, fg$parameters$A)
    ffa.error <- compute_error(fit2$A, fg$parameters$A)
  })
})



fg <- gen_fastgllvm(nobs=nobs, p=p, q=q, k=k, family=family, intercept=intercept, miss.prob = 0, scale=1, phi=rep(1, p))
a <- system.time({fgtest <- ffa(fg$Y, 10, maxiter = 100, eps=1e-10, savepath = T)})
b <- system.time({fgtest <- ffa(fg$Y, 10, maxiter = 100, eps=1e-10, savepath = T)})

fg$Y <- scale(fg$Y, scale=F)

sds <- apply(fg$Y, 2, sd)
tic()
fit2 <- ffa(fg$Y, 10, maxiter = 100, eps=1e-10, savepath = T)
b <- toc()
ts.plot(fit2$path$Psi)




