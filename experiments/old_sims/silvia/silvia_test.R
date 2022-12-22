data <- read.table("experiments/silvia/dataset.dat", header=T)

F1<-factanal(x = ~x1 + x2 + x3 + x9, factors = 1, data = data[, c(1:3, 9)])

Y<-data[,c(1:3,9)]

Y<-as.matrix(Y)
Y <- scale(Y)


set.seed(311)
fit <- fastgllvm(Y, q = 1, family="gaussian",  intercept = T, method="simple")
plot(fit)
fit <- update(fit, H=10, maxit=10)

fit2 <- ffa(Y, 1, iteratively_update_Psi = T)

plot(fit2$A, psych::Procrustes(fit$parameters$A, fit2$A)$loadings); abline(0,1,col=2)
points(fit2$Psi, fit$parameters$phi, col=2); abline(0,1,col=2)
points(fit2$A, psych::Procrustes(F1$loadings, fit2$A)$loadings, col=3)


fit$parameters$A
fit2$A
F1$loadings

fit.psych <- psych::fa(Y, 1)
