fitClumpedRegressionsPredictions<-function(calibrationData, hasMaterial=F, 
                                           returnModels=T, n.iter= 50000, burninFrac=0.1,
                                           alphaBLM1='dnorm(0.231,0.065)', betaBLM1= "dnorm(0.039,0.004)",
                                           useInits=T, 
                                           D47Prederror,
                                           D47Pred){
  
  ##Models
  BLM1<-paste(" model{
    # Diffuse normal priors for predictors
    alpha ~ ", alphaBLM1," \n ",
              "beta ~ ", betaBLM1," \n ", 
              "# Uniform prior for standard deviation
    tauy <- pow(sigma, -2)                               # precision
    sigma ~ dunif(0, 100)                                # diffuse prior for standard deviation
    
    # Diffuse normal priors for true x
    for (i in 1:N){
        x[i] ~ dnorm(0,1e-3)
    }
    
    # Likelihood
    for (i in 1:N){
        obsy[i] ~ dnorm(y[i],pow(erry[i],-2))
        y[i] ~ dnorm(mu[i],tauy)
        obsx[i] ~ dnorm(x[i],pow(errx[i],-2))
        mu[i] <- alpha+beta*x[i]
    }

    # Diffuse normal priors for true D47
    for (i in 1:NPred){
        D47[i] ~ dunif(0.5,0.8)
    }

    ##Predictions 
     for ( i in 1:NPred) {
     sigmaC[i] ~ dunif(0, 1) # standard deviation
	   tauC[i] <- 1 / (sigmaC[i] * sigmaC[i]) # convert to precision
     }
	   
	   for ( i in 1:NPred) {
	   D47Pred[i]~dnorm(D47[i],pow(D47Prederror[i],-2))
	   tw[i]<- (D47[i]-alpha)/beta
		 Tcpropagated[i] ~ dnorm(tw[i], tauC[i])
	   }
	   
}")
  

  BLM1_NoErrors<-paste("model{
                # Diffuse normal priors for predictors
                alpha ~ ", alphaBLM1," \n ",
                       "beta ~ ", betaBLM1," \n ",
                       "# Gamma prior for scatter
                 tau ~ dgamma(1e-3, 1e-3)
                epsilon <- pow(tau, -2)
                for (i in 1:N){
                y[i] ~ dnorm(mu[i],tau)
                mu[i]<- eta[i]
                eta[i] <- alpha + inprod(x[i],beta)
                }


    # Diffuse normal priors for true D47
    for (i in 1:NPred){
        D47[i] ~ dunif(0.5,0.8)
    }

    ##Predictions 
     for ( i in 1:NPred) {
     sigmaC[i] ~ dunif(0, 1) # standard deviation
	   tauC[i] <- 1 / (sigmaC[i] * sigmaC[i]) # convert to precision
     }
	   
	   for ( i in 1:NPred) {
	   D47Pred[i]~dnorm(D47[i],pow(D47Prederror[i],-2))
	   tw[i]<- (D47[i]-alpha)/beta
		 Tcpropagated[i] ~ dnorm(tw[i], tauC[i])
	   }
  
}")
  
  
  ##Mixed Model (interaction effects; multiple slopes and intercepts)
  BLM3<-paste(" model{
  
    # Diffuse normal priors for predictors
        for (i in 1:K) {
            beta[i] ~  ", alphaBLM1," \n ",
              
              " alpha[i] ~  ", betaBLM1," \n ",
              " }
        

    # Gamma prior for standard deviation
    tau ~ dgamma(1e-3, 1e-3) # precision
    sigma <- 1 / sqrt(tau) # standard deviation

    # Diffuse normal priors for true x
    for (i in 1:N){
        x1[i] ~ dnorm(0,1e-3)
    }

    # Likelihood function
    for (i in 1:N){
        obsy[i] ~ dnorm(y[i],pow(erry[i],-2))
        y[i] ~ dnorm(mu[i],tau)
        obsx1[i] ~ dnorm(x1[i],pow(errx1[i],-2))

        mu[i] <- alpha[type[i]] + beta[type[i]] * x1[i]
    }
    
    # Diffuse normal priors for true D47
    for (i in 1:NPred){
        D47[i] ~ dunif(0.5,0.8)
    }

    ##Predictions 
     for ( i in 1:NPred) {
     sigmaC[i] ~ dunif(0, 1) # standard deviation
	   tauC[i] <- 1 / (sigmaC[i] * sigmaC[i]) # convert to precision
     }
	   
	   for ( i in 1:NPred) {
	   D47Pred[i]~dnorm(D47[i],pow(D47Prederror[i],-2))
	   tw[i]<- (D47[i]-alpha)/beta
		 Tcpropagated[i] ~ dnorm(tw[i], tauC[i])
	   }
	
}")
  
  
  LM_Data <- list(obsx = calibrationData$T2 , obsy = calibrationData$D47 , 
                  errx = calibrationData$TempError, erry = calibrationData$D47error, 
                  N=nrow(calibrationData), 
                  NPred=length(D47Prederror), D47Prederror=D47Prederror,
                  D47Pred=D47Pred)
  
  LM_No_error_Data <- list(x = calibrationData$T2 , y = calibrationData$D47,
                           N=nrow(calibrationData), 
                           NPred=length(D47Prederror), D47Prederror=D47Prederror,
                           D47Pred=D47Pred)
  
  ##Fit linear models
  if(hasMaterial == T){
    
    
    ##Create the calibrationDatasets for Bayesian Models
    LM_Data <- list(obsx = calibrationData$T2 , obsy = calibrationData$D47 , 
                    errx = calibrationData$TempError, erry = calibrationData$D47error, 
                    N=nrow(calibrationData), 
                    NPred=length(D47Prederror), D47Prederror=D47Prederror,
                    D47Pred=D47Pred)
    
    ANCOVA2_Data <- list(obsx1 = calibrationData$T2 , obsy = calibrationData$D47 , 
                         errx1 = calibrationData$TempError, erry = calibrationData$D47error, 
                         K=length(unique(calibrationData$Material)),
                         N=nrow(calibrationData),
                         type= as.numeric(calibrationData$Material), 
                         NPred=length(D47Prederror), D47Prederror=D47Prederror,
                         D47Pred=D47Pred)
    
    ##Fit the models
    inits <- if(useInits==T){ function () {
      #list(alpha = rnorm(1,0,.01),
      #    beta = rnorm(1,0,.01))
      
      list(alpha = rnorm(1,0.231,0.065),
           beta = rnorm(1,0.039,0.004),
           D47=rnorm(1,0.5,0.8))
      
    }}else{NULL}
    
    BLM1_fit <- jags(data = LM_Data,inits = inits,
                     parameters = c("Tcpropagated"),
                     model = textConnection(BLM1), n.chains = 3, 
                     n.iter = n.iter, n.burnin = n.iter*burninFrac, n.thin = 10)
    #BLM1_fit <- update(BLM1_fit,n.iter = n.iter,  n.burnin = n.iter*burninFrac)
    BLM1_fit <- autojags(BLM1_fit,n.iter = n.iter,  n.burnin = n.iter*burninFrac, Rhat = 1.1, n.update = 2, n.thin = 4) 
    
    BLM1_fit_NoErrors <- jags(data = LM_No_error_Data,inits = inits,
                              parameters = c("Tcpropagated"),
                              model = textConnection(BLM1_NoErrors), n.chains = 3,
                              n.iter = n.iter,  n.burnin = n.iter*burninFrac, n.thin = 10)
    
    BLM1_fit_NoErrors <- autojags(BLM1_fit_NoErrors,n.iter = n.iter,  n.burnin = n.iter*burninFrac, Rhat = 1.1, n.update = 2, n.thin = 4)
    
    
    ##ANCOVA 2
    inits <- if(useInits==T){ function () {
      #list(alpha = rnorm(ANCOVA2_Data$K,0,0.01),
      #     beta = rnorm(ANCOVA2_Data$K, 0, 0.01))
      
      list(alpha = rnorm(ANCOVA2_Data$K,0.231,0.065),
           beta = rnorm(ANCOVA2_Data$K,0.039,0.004))
      
    }}else{NULL}
    
    BLM3_fit <- jags(data = ANCOVA2_Data,inits = inits,
                     parameters = c("Tcpropagated"), 
                     model = textConnection(BLM3), n.chains = 3,
                     n.iter = n.iter,  n.burnin = n.iter*burninFrac, n.thin = 10)
    
    #BLM3_fit <- update(BLM3_fit,n.iter = n.iter,  n.burnin = n.iter*burninFrac)
    BLM3_fit <- autojags(BLM3_fit,n.iter = n.iter,  n.burnin = n.iter*burninFrac, Rhat = 1.1, n.update = 2, n.thin = 4) 
    
    
    CompleteModelFit<-list("BLM1_fit"=BLM1_fit,'BLM1_fit_NoErrors'=BLM1_fit_NoErrors, "BLM3_fit"=BLM3_fit)
  }else{
    
    ##Fit the models
    inits <- if(useInits==T){ function () {
      #list(alpha = rnorm(1,0,.01),
      #     beta = rnorm(1,0,.01))
      list(alpha = rnorm(1,0.231,0.065),
           beta = rnorm(1,0.039,0.004))
    }}else{NULL}
    
    BLM1_fit <- jags(data = LM_Data,inits = inits,
                     parameters = c("Tcpropagated"),
                     model = textConnection(BLM1), n.chains = 3, 
                     n.iter = n.iter,  n.burnin = n.iter*burninFrac, n.thin = 10)
    #BLM1_fit <- update(BLM1_fit,n.iter = n.iter,  n.burnin = n.iter*burninFrac)
    BLM1_fit <- autojags(BLM1_fit,n.iter = n.iter,  n.burnin = n.iter*burninFrac, Rhat = 1.1, n.update = 2, n.thin = 4) 
    
    BLM1_fit_NoErrors <- jags(data = LM_No_error_Data,inits = inits,
                              parameters = c("Tcpropagated"),
                              model = textConnection(BLM1_NoErrors), n.chains = 3,
                              n.iter = n.iter,  n.burnin = n.iter*burninFrac, n.thin = 10)
    
    BLM1_fit_NoErrors <- autojags(BLM1_fit_NoErrors,n.iter = n.iter,  n.burnin = n.iter*burninFrac, Rhat = 1.1, n.update = 2, n.thin = 4)
    
    CompleteModelFit<-list('BLM1_fit'=BLM1_fit,'BLM1_fit_NoErrors'=BLM1_fit_NoErrors)
  }
  attr(CompleteModelFit, 'data') <- calibrationData 
  return(CompleteModelFit)
}

