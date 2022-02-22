fitClumpedRegressionsPredictions_D47errors<-function(calibrationData, hasMaterial=F, 
                                           n.iter= 20000, burninFrac=0.5,
                                           useInits=T, 
                                           D47Pred,
                                           D47Prederror,
                                           materialsPred,
                                           priors='informative',
                                           degC=T){
  
  if(priors == 'informative'){
    alphaBLM1='dnorm(0.231,0.065)' 
    betaBLM1= 'dnorm(0.039,0.004)'}else{
      alphaBLM1='dnorm(0, 0.01)' 
      betaBLM1= 'dnorm(0, 0.01)'
    }
  
  
  ##Models
  BLM1<-paste(" model{
    # Diffuse normal priors for predictors
    alpha ~ ", alphaBLM1," \n ",
              "beta ~ ", betaBLM1," \n ", 
              "
    sigma <- 1/sqrt(tau)                              
    tau ~ dgamma(0.01, 0.01)
    
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
    
      # prediction
  mux2 ~ dnorm(7, 0.0001)
  taux2 ~ dgamma(0.01, 0.01)

  for(i in 1:N2){ 
    y2[i] ~ dnorm(mu2[i], tau)
    mu2[i] <- alpha  + beta * x2[i]
    x2[i] ~ dnorm(mux2, taux2)
  }
    
}")
 

BLM1_NoErrors<-paste("model{
                # Diffuse normal priors for predictors
                alpha ~ ", alphaBLM1," \n ",
                       "beta ~ ", betaBLM1," \n ",
                       "
  sigma2 <- 1 / tau
  tau ~ dgamma(0.01, 0.01)
  
  # calibration
  for(i in 1:N){   
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + beta * x[i]
  }
  
  # prediction
  mux2 ~ dnorm(7, 0.0001)
  taux2 ~ dgamma(0.01, 0.01)

  for(i in 1:N2){ 
    y2[i] ~ dnorm(mu2[i], tau)
    mu2[i] <- alpha  + beta * x2[i]
    x2[i] ~ dnorm(mux2, taux2)
  }
  
}")



  ##Mixed Model (interaction effects; multiple slopes and intercepts)
BLM3<-paste(" model{
  
    # Diffuse normal priors for predictors
        for (i in 1:K) {
            beta[i] ~  ", betaBLM1 ," \n ",
            
            " alpha[i] ~  ",alphaBLM1 ," \n ",
            " }
              
    # Gamma prior for standard deviation
    tau ~ dgamma(0.1, 0.1) # precision
    sigma <- 1 / sqrt(tau) # standard deviation

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

    # prediction
  mux2 ~ dnorm(7, 0.0001)
  taux2 ~ dgamma(0.01, 0.01)

  for(i in 1:N2){ 
    y2[i] ~ dnorm(mu2[i], tau)
    mu2[i] <- alpha[typePred[i]]  + beta[typePred[i]] * x2[i]
    x2[i] ~ dnorm(mux2, taux2)
  }

}")


  ##Fit linear models
  if(hasMaterial == T){
    ANCOVA2_Data <- list(obsx1 = calibrationData$Temperature , obsy = calibrationData$D47 , 
                         errx1 = calibrationData$TempError, erry = calibrationData$D47error, 
                         K=length(unique(calibrationData$Material)),
                         N=nrow(calibrationData),
                         type= as.numeric(factor(calibrationData$Material)), 
                         typePred= as.numeric(factor(materialsPred)), 
                         N2=length(D47Pred),
                         y2=D47Pred,
                         x2=rep(NA,length(D47Pred)))
    
    BLM3_fit <- jags(data = ANCOVA2_Data,# inits = inits,
                     parameters = c("x2"), 
                     model = textConnection(BLM3), n.chains = 3,
                     n.iter = n.iter,  n.burnin = n.iter*burninFrac)
    
    BLM1_fit <- BLM3_fit
      
    BLM1_fit_NoErrors <- BLM3_fit
    
    CompleteModelFit<-list("BLM1_fit"=BLM1_fit,'BLM1_fit_NoErrors'=BLM1_fit_NoErrors, "BLM3_fit"=BLM3_fit)
  }else{
    
    LM_No_error_Data <- list(x = calibrationData$Temperature , y = calibrationData$D47,
                             N=nrow(calibrationData), 
                             N2=length(D47Pred),
                             y2=D47Pred,
                             x2=rep(NA,length(D47Pred)))
    
    BLM1_fit_NoErrors <- jags(data = LM_No_error_Data,#inits = inits,
                              parameters = c("x2"),
                              model = textConnection(BLM1_NoErrors), n.chains = 3,
                              n.iter = n.iter,  n.burnin = n.iter*burninFrac)
    
    
    LM_Data <- list(obsx = calibrationData$Temperature , obsy = calibrationData$D47 , 
                    errx = abs(calibrationData$TempError), erry = calibrationData$D47error, 
                    N=nrow(calibrationData), 
                    N2=length(D47Pred),
                    y2=D47Pred,
                    x2=rep(NA,length(D47Pred))
    )
    
    BLM1_fit <- jags(data = LM_Data,#inits = inits,
                     parameters = c("x2"),
                     model = textConnection(BLM1), n.chains = 3, 
                     n.iter = n.iter,  n.burnin = n.iter*burninFrac)
    
    
    CompleteModelFit<-list('BLM1_fit'=BLM1_fit,'BLM1_fit_NoErrors'=BLM1_fit_NoErrors)
  }
  attr(CompleteModelFit, 'data') <- calibrationData 
  return(CompleteModelFit)
}


fitClumpedRegressionsPredictions_D47<-function(calibrationData, hasMaterial=F, 
                                                     n.iter= 20000, burninFrac=0.5,
                                                     useInits=T, 
                                                     D47Pred,
                                                     D47Prederror,
                                                     materialsPred,
                                                     priors='informative',
                                                     degC=T){
  
  if(priors == 'informative'){
    alphaBLM1='dnorm(0.231,0.065)' 
    betaBLM1= 'dnorm(0.039,0.004)'}else{
      alphaBLM1='dnorm(0, 0.01)' 
      betaBLM1= 'dnorm(0, 0.01)'
    }
  
  
  ##Models
  BLM1<-paste(" model{
    # Diffuse normal priors for predictors
    alpha ~ ", alphaBLM1," \n ",
              "beta ~ ", betaBLM1," \n ", 
              "
    sigma <- 1/sqrt(tau)                              
    tau ~ dgamma(0.01, 0.01)
    
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
    
      # prediction
  mux2 ~ dnorm(7, 0.0001)
  taux2 ~ dgamma(0.01, 0.01)

  for(i in 1:N2){ 
    y2Obs[i] ~ dnorm(y2[i],pow(y2Obserr[i],-2) )
    y2[i] ~ dnorm(mu2[i], tau)
    mu2[i] <- alpha  + beta * x2[i]
    x2[i] ~ dnorm(mux2, taux2)
  }
    
}")
  
  
  BLM1_NoErrors<-paste("model{
                # Diffuse normal priors for predictors
                alpha ~ ", alphaBLM1," \n ",
                       "beta ~ ", betaBLM1," \n ",
                       "
  sigma2 <- 1 / tau
  tau ~ dgamma(0.01, 0.01)
  
  # calibration
  for(i in 1:N){   
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + beta * x[i]
  }
  
  # prediction
  mux2 ~ dnorm(7, 0.0001)
  taux2 ~ dgamma(0.01, 0.01)

  for(i in 1:N2){ 
    y2Obs[i] ~ dnorm(y2[i],pow(y2Obserr[i],-2) )
    y2[i] ~ dnorm(mu2[i], tau)
    mu2[i] <- alpha  + beta * x2[i]
    x2[i] ~ dnorm(mux2, taux2)
  }
  
}")
  
  
  
  ##Mixed Model (interaction effects; multiple slopes and intercepts)
  BLM3<-paste(" model{
  
    # Diffuse normal priors for predictors
        for (i in 1:K) {
            beta[i] ~  ", betaBLM1 ," \n ",
              
              " alpha[i] ~  ",alphaBLM1 ," \n ",
              " }
              
    # Gamma prior for standard deviation
    tau ~ dgamma(0.1, 0.1) # precision
    sigma <- 1 / sqrt(tau) # standard deviation

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

    # prediction
  mux2 ~ dnorm(7, 0.0001)
  taux2 ~ dgamma(0.01, 0.01)

  for(i in 1:N2){ 
    y2Obs[i] ~ dnorm(y2[i],pow(y2Obserr[i],-2) )
    y2[i] ~ dnorm(mu2[i], tau)
    mu2[i] <- alpha[typePred[i]]  + beta[typePred[i]] * x2[i]
    x2[i] ~ dnorm(mux2, taux2)
  }

}")
  
  
  ##Fit linear models
  if(hasMaterial == T){
    ANCOVA2_Data <- list(obsx1 = calibrationData$Temperature , obsy = calibrationData$D47 , 
                         errx1 = calibrationData$TempError, erry = calibrationData$D47error, 
                         K=length(unique(calibrationData$Material)),
                         N=nrow(calibrationData),
                         type= as.numeric(factor(calibrationData$Material)), 
                         typePred= as.numeric(factor(materialsPred)), 
                         N2=length(D47Pred),
                         y2Obs=D47Pred,
                         y2Obserr=D47Prederror,
                         x2=rep(NA,length(D47Pred)))
    
    BLM3_fit <- jags(data = ANCOVA2_Data,# inits = inits,
                     parameters = c("x2"), 
                     model = textConnection(BLM3), n.chains = 3,
                     n.iter = n.iter,  n.burnin = n.iter*burninFrac)
    
    BLM1_fit <- BLM3_fit
    
    BLM1_fit_NoErrors <- BLM3_fit
    
    CompleteModelFit<-list("BLM1_fit"=BLM1_fit,'BLM1_fit_NoErrors'=BLM1_fit_NoErrors, "BLM3_fit"=BLM3_fit)
  }else{
    
    LM_No_error_Data <- list(x = calibrationData$Temperature , y = calibrationData$D47,
                             N=nrow(calibrationData), 
                             N2=length(D47Pred),
                             y2Obs=D47Pred,
                             y2Obserr=D47Prederror,
                             x2=rep(NA,length(D47Pred)))
    
    BLM1_fit_NoErrors <- jags(data = LM_No_error_Data,#inits = inits,
                              parameters = c("x2"),
                              model = textConnection(BLM1_NoErrors), n.chains = 3,
                              n.iter = n.iter,  n.burnin = n.iter*burninFrac)
    
    
    LM_Data <- list(obsx = calibrationData$Temperature , obsy = calibrationData$D47 , 
                    errx = abs(calibrationData$TempError), erry = calibrationData$D47error, 
                    N=nrow(calibrationData), 
                    N2=length(D47Pred),
                    y2Obs=D47Pred,
                    y2Obserr=D47Prederror,
                    x2=rep(NA,length(D47Pred))
    )
    
    BLM1_fit <- jags(data = LM_Data,#inits = inits,
                     parameters = c("x2"),
                     model = textConnection(BLM1), n.chains = 3, 
                     n.iter = n.iter,  n.burnin = n.iter*burninFrac)
    
    
    CompleteModelFit<-list('BLM1_fit'=BLM1_fit,'BLM1_fit_NoErrors'=BLM1_fit_NoErrors)
  }
  attr(CompleteModelFit, 'data') <- calibrationData 
  return(CompleteModelFit)
}


