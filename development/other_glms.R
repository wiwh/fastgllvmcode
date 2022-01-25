devtools::load_all()
family <- binomial()

n <- 100

q <- 1
p <- 2000
k <- 20
set.seed(1241)
fg <- gen_fastgllvm(n = n, p=p, q=q, k=k,family = family, intercept = T)

fit <- fastgllvm(fg$Y, q=q, fg$X, method="SA", H=1, family=family,
                 A.init <- diag(1, p,q),
                   learning_rate = "spall", maxit=100,
                 learning_rate.args = list(start=4e1, end=1e-2), tol=1e-6)

par(mfrow=c(2,1))
ts.plot(fit$A.hist)
ts.plot(fit$B.hist)
plot(fg$A, psych::Procrustes(fit$A, fg$A)$loadings)
abline(0,1,col=2)
plot(fg$B, fit$B)
abline(0,1,col=2)
par(mfrow=c(1,1))


# Zhat <- predict(fg, method="ridge", lambda=.001)
# plot(-fg$Z, Zhat); abline(0,1,col=2)

n <- fg$n
H <- 1
A <- fg$A
B <- fg$B
phi <- fg$phi
q <- fg$q
p <- fg$p
method <- "SA"

H <- 10
generate_Z <- generate_Z_functionfactory(n, q, H, method)
Y <- fg$Y
Y.c <- scale(Y, scale=F)
X <- fg$X
family <- fg$family

# sims <- lapply(1:1000, function(na) get_Psi(Y, Y.c, A, B, phi, X, family, generate_Z))
# A.sim <- t(sapply(sims, function(simi)simi$A))
# boxplot(A.sim)


K <- get_K(A, B, phi)
E <- get_expectations(A, B, X, generate_Z(), K, phi, family)
plot(t(Y.c) %*% (Y.c %*% t(K))/nrow(Y.c), E$EYYK); abline(0,1,col=2)
