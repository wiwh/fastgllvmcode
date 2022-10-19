devtools::load_all()

nobs <- 100
p <- 1000
q <- 10
family <- "binomial"
set.seed(123132)
a <- gen_fastgllvm(nobs=nobs, p=p, q=q, family=family)

sims <- lapply(1:100, function(i) {
  new_gllvm <- with(a, gen_fastgllvm(nobs = nobs, p=p, q=q, k=1, family = family, A=parameters$A,B=parameters$B, X=X, Z=parameters$Z, phi=parameters$phi))
  Zstar <- with(new_gllvm, compute_zstar(Y, parameters$A, parameters$phi, X, parameters$B, families)$Zstar)
  # Zstar.2 <- scale(Zstar, center=F)
  # Zstar <- (Zstar + Zstar.2 )/2
  k2 <- as.vector(t(new_gllvm$Y) %*% (Zstar)/nrow(new_gllvm$Y))
  k2
})

sims <- do.call(rbind, sims)

sims <- sims[, order(as.vector(a$parameters$A))]

lp <- with(a, 1/(1+exp(-compute_linpar(parameters$Z, parameters$A, X, parameters$B)$linpar)))
# lp <- with(a, compute_linpar(parameters$Z, parameters$A, X, parameters$B)$linpar)

Zstar <- with(a, compute_zstar(Y, parameters$A, parameters$phi, X, parameters$B, families)$Zstar)

k1 <- as.vector((t(lp) %*% a$parameters$Z)/nrow(a$Y))
k1 <- k1[order(as.vector(a$parameters$A))]

s <- sort(sample(p*q, 100, replace=T))
boxplot(sims[,s])
points(1:length(s), k1[s], col=2)


qqplot(colMeans(sims), k1); abline(0,1,col=2)

# scale_quantiles <- function(z){
#   apply(Z,2,  function(zi){
#     n <- length(zi)
#     qn <- qnorm((1:n)/(n+1))
#     zi <- qn[rank(z)]
#
#   })
# }
#
# scale_quantiles(a$Z)
