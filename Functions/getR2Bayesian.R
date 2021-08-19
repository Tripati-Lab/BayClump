

getR2Bayesian <- function(model = NULL, calibrationData=NULL, hasMaterial=F) {
  if(hasMaterial==F){
  mcmc <- model$BUGSoutput$sims.matrix
  Xmat = model.matrix( ~ Temperature, calibrationData)
  coefs = mcmc[, c("alpha", "beta")]
  fit = coefs %*% t(Xmat)
  resid = sweep(fit, 2, calibrationData$D47, "-")
  var_f = apply(fit, 1, var)
  var_e = apply(resid, 1, var)
  R2 = var_f / (var_f + ifelse(is.na(var_e), 0,var_e))
  cbind.data.frame(
    median = median(R2),
    lwr = quantile(R2, 0.025),
    upr = quantile(R2, 0.975)
  )
  }else{
    mcmc <- model$BUGSoutput$sims.matrix
    Xmat = model.matrix( ~ Temperature, calibrationData)
    coefs = as.data.frame(mcmc[,-ncol(mcmc) ])
    new <- rowSums(coefs[,grep('alpha', names(coefs)), drop=FALSE])
    new2 <- rowSums(coefs[,grep('beta', names(coefs)), drop=FALSE])
    coefs<-cbind(new,new2)
    fit = coefs %*% t(Xmat)
    resid = sweep(fit, 2, calibrationData$D47, "-")
    var_f = apply(fit, 1, var)
    var_e = apply(resid, 1, var)
    R2 = var_f / (var_f + var_e)
    cbind.data.frame(
      median = median(R2),
      lwr = quantile(R2, 0.025),
      upr = quantile(R2, 0.975)
    )
    
  }
}
