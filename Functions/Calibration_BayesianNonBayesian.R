#' This function generate temperature predictions (in 10^6/T2) based on a 
#' calibration dataset and target D47. Note that this alternative function
#' propagates uncertainty around the target D47.
#' 
#' @param calibrationData The calibration dataset
#' @param n.iter number of MCMC iterations
#' @param priors Informative, difusse, or NonInformative on the beta and alpha
#' @param D47error the column in calibrationData containing the uncertainty in D47
#' @param init.values Use initial values for runs in JAGS?


fitClumpedRegressions <<- function(calibrationData, 
                                n.iter = 5000, 
                                priors = "Informative",
                                D47error = "D47error",
                                samples=NULL,
                                init.values = FALSE){
  
  if(! priors %in% c("Informative", "Difusse", "NonInformative") ){ 
    stop("Priors must be in `Informative`, `Difusse` or `NonInformative`")
  }
  
  if(is.null(samples)){
    warning("Using the full dataset in the calibration step.")
    samples <- nrow(calibrationData)
  }else{
    warning("Sampling ", samples, " from the dataset.")
    calibrationData <- calibrationData[sample(1:nrow(calibrationData), samples), ]
  }
     
  
  if(priors == "Informative"){
    alphaBLM1 = "dnorm(0.231,0.065)" 
    betaBLM1 = "dnorm(0.039,0.004)"}
  
  if(priors == "Difusse"){
    alphaBLM1 = "dnorm(0, 0.01)" 
    betaBLM1 = "dnorm(0, 0.01)"
  }
  
  if(priors == "NonInformative"){
    alphaBLM1 = "dnorm(0.231, 0.195)" 
    betaBLM1 = "dnorm(0.231, 0.012)"
  }
  
  
  ##Models
  BLM1 <- paste(" model{
    # Diffuse normal priors for predictors
    alpha ~ ", alphaBLM1," \n ",
              "beta ~ ", betaBLM1," \n ", 
              "
    tau <- pow(sigma, -2) 
    sigma ~ dunif(0, 100)                             
    
    for (i in 1:N){
        x[i] ~ dnorm(11,0.01)
    }
    # Likelihood
    for (i in 1:N){
        obsy[i] ~ dnorm(y[i],pow(erry[i],-2))
        y[i] ~ dnorm(mu[i],tau)
        obsx[i] ~ dnorm(x[i],pow(errx[i],-2))
        mu[i] <- alpha + beta*x[i]
    }
    
    ##Log-likelihood
  for(i in 1:N){ 
   regression_residual[i] <- y[i] - mu[i]
   zloglik[i] <- logdensity.norm(y[i], mu[i], tau)
  }
}")
  
  
  BLM1_NoErrors <- paste("model{
                # Diffuse normal priors for predictors
                alpha ~ ", alphaBLM1," \n ",
                       "beta ~ ", betaBLM1," \n ",
                       "
    tau <- pow(sigma, -2) 
    sigma ~ dunif(0, 100)                             
    
  # calibration
  for(i in 1:N){   
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + beta * x[i]
  }
  
  ##Log-likelihood
  for(i in 1:N){ 
   regression_residual[i] <- y[i] - mu[i]
   zloglik[i] <- logdensity.norm(y[i], mu[i], tau)
  }
}")
  
  ##Mixed Model (interaction effects; multiple betas and alphas)
  BLM3 <- paste(" model{
  
    # Diffuse normal priors for predictors
        for (i in 1:K) {
            beta[i] ~  ", betaBLM1 ," \n ",
              
              " alpha[i] ~  ",alphaBLM1 ," \n ",
              " }
              
    # Prior for standard deviation
    tau <- pow(sigma, -2) 
    sigma ~ dunif(0, 100)                             
    
    # Diffuse normal priors for true x
    for (i in 1:N){
        x1[i] ~ dnorm(11,0.01)
    }

    # Likelihood function
    for (i in 1:N){
        obsy[i] ~ dnorm(y[i],pow(erry[i],-2))
        y[i] ~ dnorm(mu[i],tau)
        obsx1[i] ~ dnorm(x1[i],pow(errx1[i],-2))

        mu[i] <- alpha[type[i]] + beta[type[i]] * x1[i]
    }
    
    ##R2s (mod from)
    #http://samcarcagno.altervista.org/stat_notes/r2_lmm_jags/r_squared_lmm.html
  for (s in 2:N){
    subjEff[s] ~ dnorm(0, 1/zSubjEffSigma^2)
  }
  subjEff[1] <- -sum(subjEff[2:N]) #sum-to-zero constraint
  zSubjEffSigma ~ dunif(0.00001, 10)
  subjEffSigma <- zSubjEffSigma*sd(y)

  for (j in 1:N){
     yPredFixed[j] <-  sum(alpha[type[j]] + beta[type[j]] * x1[j])
  }

  varFixed <- (sd(yPredFixed))^2
  varResidual <- sigma^2 # get the variance of residuals
  varRandom <- subjEffSigma^2  # get the variance of random plot effect
  # calculate marginal R^2
  marginalR2 <- varFixed / (varFixed + varRandom + varResidual)
  # calculate conditional R^2
  conditionalR2 <- (varRandom + varFixed) / (varFixed + varRandom + varResidual) 

##Log-likelihood
  for(i in 1:N){ 
   regression_residual[i] <- y[i] - mu[i]
   zloglik[i] <- logdensity.norm(y[i], mu[i], tau)
  }

}")
  
  #Data

    ANCOVA2_Data <- list(obsx1 = calibrationData$Temperature, 
                         obsy = calibrationData$D47 , 
                         errx1 = abs(calibrationData$TempError), 
                         erry = calibrationData[,D47error], 
                         K = length(unique(calibrationData$Material)),
                         N = nrow(calibrationData),
                         type = as.numeric(calibrationData$Material))
    
    LM_Data <- list(obsx = calibrationData$Temperature, 
                    obsy = calibrationData$D47 , 
                    errx = abs(calibrationData$TempError), 
                    erry = calibrationData[,D47error], 
                    N = nrow(calibrationData))
    
    LM_No_error_Data <- list(x = calibrationData$Temperature, 
                             y = calibrationData$D47,
                             N = nrow(calibrationData))
    
    
    #Inits
    initsMixed <- function () {
      list(alpha = rnorm(ANCOVA2_Data$K,0.231,0.065),
           beta = rnorm(ANCOVA2_Data$K,0.039,0.004))
      
    }
    
    initsSimple <- function () {
      list(alpha = rnorm(1,0.231,0.065),
           beta = rnorm(1,0.039,0.004))
      
    }
    
    init.values
    
    #Fit models
    BLM3_fit <- jags(data = ANCOVA2_Data, 
                     inits = if(init.values){initsMixed}else{NULL},
                     parameters = c("alpha","beta","conditionalR2", "marginalR2", "tau"), 
                     model = textConnection(BLM3), 
                     n.chains = 3,
                     n.iter = n.iter)
    BLM3_fit <- autojags(BLM3_fit)
    
    BLM1_fit <- jags(data = LM_Data, 
                     inits = if(init.values){initsSimple}else{NULL},
                     parameters = c("alpha","beta", "tau"),
                     model = textConnection(BLM1), 
                     n.chains = 3, 
                     n.iter = n.iter)
    BLM1_fit <- autojags(BLM1_fit)
    
    BLM1_fit_NoErrors <- jags(data = LM_No_error_Data, 
                              inits = if(init.values){initsSimple}else{NULL},
                              parameters = c("alpha","beta", "tau"),
                              model = textConnection(BLM1_NoErrors), 
                              n.chains = 3,
                              n.iter = n.iter)
    BLM1_fit_NoErrors <- autojags(BLM1_fit_NoErrors)
    
    #Extract relevant descriptors
    R2sComplete<-rbind.data.frame(getR2Bayesian(BLM1_fit, calibrationData=calibrationData),
                                  getR2Bayesian(BLM1_fit_NoErrors, calibrationData=calibrationData),
                                  getR2Bayesian(BLM3_fit, calibrationData=calibrationData,hasMaterial=T)
                                  )
    
    R2sComplete <- cbind.data.frame(model=c("BLM1_fit", "BLM1_fit_NoErrors", "BLM3_fit"), R2sComplete)
    
    aMErrors <- sum(dic.samples(BLM1_fit$model, n.iter=n.iter, thin = 1, "pD")$deviance) 
    aMNoErrors <- sum(dic.samples(BLM1_fit_NoErrors$model, n.iter=n.iter, thin = 1, "pD")$deviance)
    aM <- sum(dic.samples(BLM3_fit$model, n.iter=n.iter, thin = 1, "pD")$deviance) 
    DICs<-c(aMErrors, aMNoErrors, aM)
    names(DICs)<-c("BLM1_fit", "BLM1_fit_NoErrors", "BLM3_fit")
    
    CompleteModelFit<-list("BLM1_fit"=BLM1_fit,"BLM1_fit_NoErrors"=BLM1_fit_NoErrors, "BLM3_fit"=BLM3_fit)

    attr(CompleteModelFit, "data") <- calibrationData 
    attr(CompleteModelFit, "R2s") <- R2sComplete 
    attr(CompleteModelFit, "DICs") <- DICs 
  
  return(CompleteModelFit)
}


#' Bootstrap York regression models from a calibration dataset
#' 
#' @param data The calibration dataset
#' @param replicates Number of bootstrap replicates
#' @param samples Number of samples per replicate
#' @param D47error The column in data containing the errors in D47

simulateYork_measured <<- function(data, 
                                   replicates, 
                                   samples = NULL, 
                                   D47error = "D47error"){
  do.call(rbind,lapply(1:replicates, function(x){
    dataSub <- data[sample(seq_along(data[,1]), if(is.null(samples)){nrow(data)}else{samples}, replace = T),]
    dataSub$y_SE <- dataSub[,D47error]
    dataSub$x_SE <- abs(dataSub$TempError)
    Reg <- york(cbind.data.frame(dataSub$Temperature, dataSub$x_SE, dataSub$D47, dataSub$y_SE))
    cbind.data.frame("alpha"=Reg$a[1],"beta"=Reg$b[1])
  }))
}


#' Bootstrap an OLS regression models from a calibration dataset
#' 
#' @param data The calibration dataset
#' @param replicates Number of bootstrap replicates
#' @param samples Number of samples per replicate
#' @param D47error The column in data containing the errors in D47

simulateLM_measured <<- function(data, 
                               replicates, 
                               samples = NULL, 
                               D47error="D47error"){
  
  a<-lapply(1:replicates, function(x){
    dataSub<-data[sample(seq_along(data[,1]), if(is.null(samples)){nrow(data)}else{samples}, replace = T),]
    dataSub$y_SE<-dataSub[,D47error]
    dataSub$x_SE<-dataSub$TempError
    Reg<-summary(lm(D47 ~ Temperature,  dataSub))
    res<-cbind.data.frame("alpha"=Reg$coefficients[1,1],"beta"=Reg$coefficients[2,1])
    attr(res, "R2") <- Reg$r.squared
    res
  })
  
  R2s<-unlist(lapply(a, function(x) attributes(x)$R2))
  R2s<-data.frame(median=median(R2s), lwr=quantile(R2s, 0.025), upr=quantile(R2s, 0.975))
  a<-do.call(rbind,a)
  attr(a, "R2") <- R2s
  return(a)
}


#' Bootstrap a weighted OLS regression models from a calibration dataset
#' 
#' @param data The calibration dataset
#' @param replicates Number of bootstrap replicates
#' @param samples Number of samples per replicate
#' @param D47error The column in data containing the errors in D47


simulateLM_inverseweights <<- function(data, 
                                     replicates, 
                                     samples = NULL, 
                                     D47error="D47error"){
  a<-lapply(1:replicates, function(x){
    dataSub<-data[sample(seq_along(data[,1]), if(is.null(samples)){nrow(data)}else{samples}, replace = T),]
    dataSub$y_SE<-dataSub[,D47error]
    dataSub$x_SE<-abs(dataSub$TempError)
    Reg0<-lm(D47 ~ Temperature,  dataSub)
    wt <- 1 / lm(abs(Reg0$residuals) ~ Reg0$fitted.values)$fitted.values^2
    Reg<-summary(lm(D47 ~ Temperature,  dataSub, weights = wt))
    res<-cbind.data.frame("alpha"=Reg$coefficients[1,1],"beta"=Reg$coefficients[2,1])
    attr(res, "R2") <- Reg$r.squared
    res
  })
  
  R2s<-unlist(lapply(a, function(x) attributes(x)$R2))
  R2s<-data.frame(median=median(R2s), lwr=quantile(R2s, 0.025), upr=quantile(R2s, 0.975))
  a<-do.call(rbind,a)
  attr(a, "R2") <- R2s
  return(a)
}


#' Bootstrap Deming regression models from a calibration dataset
#' 
#' @param data The calibration dataset
#' @param replicates Number of bootstrap replicates
#' @param samples Number of samples per replicate
#' @param D47error The column in data containing the errors in D47


simulateDeming <<- function(data, 
                          replicates, 
                          samples = NULL, 
                          D47error="D47error", 
                          multicore=FALSE){
  
  if(multicore){
  
  ncores = parallel::detectCores()
  do.call(rbind,mclapply(1:replicates, mc.cores = ncores, function(x){
    dataSub<-data[sample(seq_along(data[,1]), if(is.null(samples)){nrow(data)}else{samples}, replace = T),]
    dataSub$y_SE<-abs(dataSub[,D47error])/sqrt(nrow(dataSub))
    dataSub$x_SE<-abs(dataSub$TempError)/sqrt(nrow(dataSub))
    Reg<-deming(D47 ~ Temperature, dataSub, xstd= 1/x_SE^2, ystd= 1/y_SE^2)
    cbind.data.frame("alpha"=Reg$coefficients[1],"beta"=Reg$coefficients[2])
  }))
  
  }else{
    do.call(rbind,lapply(1:replicates, function(x){
      dataSub<-data[sample(seq_along(data[,1]), if(is.null(samples)){nrow(data)}else{samples}, replace = T),]
      dataSub$y_SE<-abs(dataSub[,D47error])/sqrt(nrow(dataSub))
      dataSub$x_SE<-abs(dataSub$TempError)/sqrt(nrow(dataSub))
      Reg<-deming(D47 ~ Temperature, dataSub, xstd= 1/x_SE^2, ystd= 1/y_SE^2)
      cbind.data.frame("alpha"=Reg$coefficients[1],"beta"=Reg$coefficients[2])
    }))
  }
  
}

#' Estimate the R2 for Bayesian models
#' 
#' @param model The model (from R2jags)
#' @param calibrationData the calibration dataset used to fit the model
#' @param hasMaterial Is it the mixed model?

getR2Bayesian <<- function(model, 
                           calibrationData, 
                           hasMaterial=FALSE) {
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
      mean = mean(R2),
      lwr = quantile(R2, 0.025),
      upr = quantile(R2, 0.975)
    )
  }else{
    mcmc <- model$BUGSoutput$sims.matrix
    Xmat = model.matrix( ~ Temperature, calibrationData)
    coefs = as.data.frame(mcmc[,-ncol(mcmc) ])
    new <- rowSums(coefs[,grep("alpha", names(coefs)), drop=FALSE])
    new2 <- rowSums(coefs[,grep("beta", names(coefs)), drop=FALSE])
    coefs<-cbind(new,new2)
    fit = coefs %*% t(Xmat)
    resid = sweep(fit, 2, calibrationData$D47, "-")
    var_f = apply(fit, 1, var)
    var_e = apply(resid, 1, var)
    R2 = var_f / (var_f + var_e)
    cbind.data.frame(
      mean = mean(R2),
      lwr = quantile(R2, 0.025),
      upr = quantile(R2, 0.975)
    )
    
  }
}


#' This function is used to generate CI estimates at given intervals. It is currently
#' used for plotting in BayClump.
#' 
#' @param data A data.frame with two columns labeled beta and alpha. 
#' This should be the result of bootstrapping a given calibration set.
#' @param from the lower limit in x
#' @param to the upper limit in x
#' @param length.out the number of breaks


RegressionSingleCI<-function(data, from, to, length.out=100){
  
  sampleDataReplicates<- as.data.frame(data)
  
  Theta_mat <- sampleDataReplicates
  
  # Points where to evaluate the model
  x_eval <- seq(from, to, length.out = length.out)
  
  # Matrix with the predictions
  fun <- function(x, theta) as.numeric(theta["alpha"]) + (x * as.numeric(theta["beta"]))
  Pred_mat <- apply(Theta_mat, 1, function(theta) fun(x_eval, theta))
  
  
  # Pack the estimates for plotting
  uncertaintyModels <- cbind(
    x = x_eval, 
    as.data.frame(t(apply(Pred_mat, 1, function(y_est) c(
      median_est = median(y_est), 
      ci_lower_est = quantile(y_est, probs = 0.025, names = FALSE), 
      ci_upper_est = quantile(y_est, probs = 0.975, names = FALSE)
    ))))
  )
  
  return(list(uncertaintyModels))
  
}

#' Generate a dataset reflecting the priors used to run the analyses
#' 
#' @param prior Informative or not
#' @param n number of observations to simulate

generatePriorDistCalibration <- function(prior, n=1000){
  if(prior == "Informative"){
    params <- cbind.data.frame(parameter=c("alpha", "beta"),
                               mean=c(0.231,0.039), 
                               sd=c(0.065,0.004))
    params
    } else {
      params <- cbind.data.frame(parameter=c("alpha", "beta"),
                                 mean=c(0,0.01), 
                                 sd=c(0,0.01))
      params
    }
  
  data <- cbind.data.frame(alpha=rnorm(n, params[1,2], params[1,3]), 
                   beta=rnorm(n, params[2,2], params[2,3]))
  attr(data, "priors") <-  prior
  attr(data, "params") <-  params
  data
}



