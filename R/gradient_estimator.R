compute_gradients_glm <- function(Y, X, parameters, families, Miss, ...) {
  if(!is.null(parameters$B) && all(X[,1]==1)) rescale.B=1 else rescale.B=FALSE
  # recenter Z and B
  parameters_sim <- parameters_sam <- recenter(parameters, 1)

  # Generate sim
  Y_sim <- generate_y(
    linpar = NULL,
    phi = parameters$phi,
    families = families,
    A = parameters$A,
    B = parameters$B,
    X = X,
    Z = NULL,
    nobs = nrow(Y)
  )

  # Compute psi for sam
  parameters_sam$Z <- compute_zstar(Y, parameters_sam$A, parameters_sam$phi, X, parameters_sam$B, families, start=parameters_sam$Z, Miss=Miss)$Zstar
  # parameters_sam <- recenter(parameters_sam, intercept.id=1)

  ZX_sam <- ZX_join(parameters_sam$Z, X)
  est_sam <- sapply(1:ncol(Y), function(j){
    glm(Y[,j] ~ 0 + ZX_sam, family=families$vec[j])$coef
  }, simplify=F)
  est_sam <- do.call(rbind, est_sam)
  AB <- AB_separate(est_sam, ncol(parameters$A))
  list(A=AB$A, B=AB$B)
}



if(0) {
  devtools::load_all()
  set.seed(1234)
  poisson  <- 10
  gaussian <- 5
  binomial <- 25
  q <- 4
  p <- poisson + gaussian + binomial
  family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial))
  set.seed(120303)
  fg <- gen_fastgllvm(nobs=1000, p=p, q=q, family=family, phi=3*(1:p)/p, k=1, intercept=T, miss.prob = 0, scale=1)
  psi <- compute_psi_simple(fg$Y, fg$X, fg$parameters$Z, fg$parameters, fg$families, fg$Miss)
  gradient <- compute_gradients_glm(fg$Y, fg$X, fg$parameters, fg$families, fg$Miss)
  plot(fg$parameter$A, psych::Procrustes(gradient$A,fg$parameters$A)$loadings)

  # this must have expectaiton 0 under the true model. check this
  sim <- sapply(1:1000, function(i){
    set.seed(i)
    if(i%%10==0 || i==1) cat("\n", i)
    fg <- gen_fastgllvm(nobs=1000, p=p, q=q, A= fg$parameters$A, B = fg$parameters$B, phi=fg$parameters$phi, family=family, k=1, intercept=T, miss.prob = 0, scale=1)
    psi <- compute_gradients_simple_rescale(fg$Y, fg$X, fg$parameters, fg$families, fg$Miss)
    psi
  }, simplify=F)

  sim <- sapply(names(sim[[1]]),
                function(parname) do.call(rbind, sapply(sim, function(simi) as.vector(simi[[parname]]),
                                                        simplify=F)), simplify=F)

  boxplot(sim$A, outline=F); abline(h=0, col=2)
  boxplot(sim$B,outline=F); abline(h=0, col=2)
  boxplot(sim$phi);abline(h=0, col=2)
}
