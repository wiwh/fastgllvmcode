devtools::load_all()


fg <- gen_fastgllvm(family = "binomial", q=2, p=20)
Y <- simulate_fastgllvm(fg)$Y

fit <- fg

set.seed(12431)
fit <- fit_fastgllvm(fit, method="SP", verbose=T)
plot_fastgllvm(fit)


A <- NULL
fit <- fg

for(i in 1:10){
  fit <- fit_fastgllvm(fit, H=20, method="SP", verbose=T)
  A <- rbind(A, as.vector(fit$A))
}

for(i in 1:10){
  fit <- fit_fastgllvm(fit, method="SA", verbose=T)
  A <- rbind(A, as.vector(fit$A))
}


