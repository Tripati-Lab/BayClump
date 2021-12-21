fitClumpedRegressionsPredictions<-function(calibrationData, hasMaterial=F, 
                                           n.iter= 20000, burninFrac=0.5,
                                           alphaBLM1='dnorm(0.231,0.065)', 
                                           betaBLM1= 'dnorm(0.039,0.004)',
                                           useInits=T, 
                                           D47Pred,
                                           D47Prederror,
                                           materialsPred){
  
  ##Models
  BLM1<-paste(" model{
    # Diffuse normal priors for predictors
    alpha ~ ", alphaBLM1," \n ",
              "beta ~ ", betaBLM1," \n ", 
              "
    sigma <- 1/sqrt(tauy)                              
    tauy ~ dgamma(0.1, 0.1)   
    
    # Diffuse normal priors for true x
    for (i in 1:N){
        x[i] ~ dnorm(11,0.01)
    }
    
    # Likelihood
    for (i in 1:N){
        obsy[i] ~ dnorm(y[i],pow(erry[i],-2))
        y[i] ~ dnorm(mu[i],tauy)
        obsx[i] ~ dnorm(x[i],pow(errx[i],-2))
        mu[i] <- alpha+beta*x[i]
    }
      
    #T_C = sqrt((beta * 10^6) / (y - alpha)) - 273.15

    # Diffuse normal priors for true D47 pred
    for (i in 1:NPred){
        truepredD47[i] ~ dnorm(0.6,0.001)
    }

	   for ( i in 1:NPred) {
	   D47Pred[i] ~ dnorm(truepredD47[i], pow(D47Prederror[i],-2))
	   tw[i] <-  (truepredD47[i] - alpha)/beta
		 Tcpropagated[i] ~ dnorm(tw[i], tauy)
	   }
	   
}")
  
  
  BLM1_NoErrors<-paste("model{
                # Diffuse normal priors for predictors
                alpha ~ ", alphaBLM1," \n ",
                       "beta ~ ", betaBLM1," \n ",
                       "
    sigma <- 1/sqrt(tau)                              
    tau ~ dgamma(0.1, 0.1)  
                for (i in 1:N){
                y[i] ~ dnorm(mu[i],tau)
                mu[i]<- eta[i]
                eta[i] <- alpha + inprod(x[i],beta)
                }
 
     # Diffuse normal priors for true D47 pred
    for (i in 1:NPred){
        truepredD47[i] ~ dnorm(0.6,0.001)
    }

	   for ( i in 1:NPred) {
	   D47Pred[i] ~ dnorm(truepredD47[i], pow(D47Prederror[i],-2))
	   tw[i] <- (truepredD47[i] - alpha)/beta
		 Tcpropagated[i] ~ dnorm(tw[i], tau)
	   }
  
}")
  
  
  ##Mixed Model (interaction effects; multiple slopes and intercepts)
  BLM3<-paste(" model{
  
    # Diffuse normal priors for predictors
        for (i in 1:K) {
            beta[i] ~  ",betaBLM1 ," \n ",
              
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
        y[i] ~ dnorm(mu[i], tau)
        obsx1[i] ~ dnorm(x1[i],pow(errx1[i],-2))
        mu[i] <- alpha[type[i]] + beta[type[i]] * x1[i]
    }
    
	   # Diffuse normal priors for true D47 pred
    for (i in 1:NPred){
        truepredD47[i] ~ dnorm(0.6,0.001)
    }

	   for ( i in 1:NPred) {
	   D47Pred[i] ~ dnorm(truepredD47[i], pow(D47Prederror[i],-2))
	   tw[i] <- (truepredD47[i] - alpha[typePred[i]])/beta[typePred[i]]
		 Tcpropagated[i] ~ dnorm(tw[i], tau)
	   }
}")
  
  
  
  LM_Data <- list(obsx = calibrationData$Temperature , obsy = calibrationData$D47 , 
                  errx = calibrationData$TempError, erry = calibrationData$D47error, 
                  N=nrow(calibrationData), 
                  NPred=length(D47Pred), 
                  D47Pred=D47Pred,
                  D47Prederror= D47Prederror
                  )
  
  LM_No_error_Data <- list(x = calibrationData$Temperature , y = calibrationData$D47,
                           N=nrow(calibrationData), 
                           NPred=length(D47Pred),
                           D47Pred=D47Pred,
                           D47Prederror= D47Prederror)
  
  ##Fit linear models
  if(hasMaterial == T){
    
    
    ##Create the calibrationDatasets for Bayesian Models
    LM_Data <- list(obsx = calibrationData$Temperature , obsy = calibrationData$D47 , 
                    errx = calibrationData$TempError, erry = calibrationData$D47error, 
                    N=nrow(calibrationData), 
                    NPred=length(D47Pred),
                    D47Prederror= D47Prederror,
                    D47Pred=D47Pred)
    
    ANCOVA2_Data <- list(obsx1 = calibrationData$Temperature , obsy = calibrationData$D47 , 
                         errx1 = calibrationData$TempError, erry = calibrationData$D47error, 
                         K=length(unique(calibrationData$Material)),
                         N=nrow(calibrationData),
                         type= as.numeric(factor(calibrationData$Material)), 
                         NPred=length(D47Pred),
                         D47Pred=D47Pred,
                         D47Prederror= D47Prederror,
                         typePred= as.numeric(factor(materialsPred)))
    
    ##Fit the models
    inits <- if(useInits==T){ function () {
      list(alpha = rnorm(1,0.231,0.065),
           beta = rnorm(1,0.039,0.004),
           truepredD47 = rnorm(LM_Data$NPred,0.6,0.01) )
      
    }}else{NULL}
    
    BLM1_fit <- jags(data = LM_Data,inits = inits,
                     parameters = c("Tcpropagated"),
                     model = textConnection(BLM1), n.chains = 3, 
                     n.iter = n.iter, n.burnin = n.iter*burninFrac)
    
    BLM1_fit_NoErrors <- jags(data = LM_No_error_Data,inits = inits,
                              parameters = c("Tcpropagated"),
                              model = textConnection(BLM1_NoErrors), n.chains = 3,
                              n.iter = n.iter,  n.burnin = n.iter*burninFrac)
    
    
    ##ANCOVA 2
    inits <- if(useInits==T){ function () {
      list(alpha = rnorm(ANCOVA2_Data$K,0.231,0.065),
           beta = rnorm(ANCOVA2_Data$K,0.039,0.004),
           truepredD47 = rnorm(LM_Data$NPred,0.6,0.01))
    }}else{NULL}
    
    BLM3_fit <- jags(data = ANCOVA2_Data, inits = inits,
                     parameters = c("Tcpropagated"), 
                     model = textConnection(BLM3), n.chains = 3,
                     n.iter = n.iter,  n.burnin = n.iter*burninFrac)
    
    CompleteModelFit<-list("BLM1_fit"=BLM1_fit,'BLM1_fit_NoErrors'=BLM1_fit_NoErrors, "BLM3_fit"=BLM3_fit)
  }else{
    
    ##Fit the models
    inits <- if(useInits==T){ function () {
      list(alpha = rnorm(1,0.231,0.065),
           beta = rnorm(1,0.039,0.004))
    }}else{NULL}
    
    BLM1_fit <- jags(data = LM_Data,inits = inits,
                     parameters = c("Tcpropagated"),
                     model = textConnection(BLM1), n.chains = 3, 
                     n.iter = n.iter,  n.burnin = n.iter*burninFrac)
    
    BLM1_fit_NoErrors <- jags(data = LM_No_error_Data,inits = inits,
                              parameters = c("Tcpropagated"),
                              model = textConnection(BLM1_NoErrors), n.chains = 3,
                              n.iter = n.iter,  n.burnin = n.iter*burninFrac)
    
    CompleteModelFit<-list('BLM1_fit'=BLM1_fit,'BLM1_fit_NoErrors'=BLM1_fit_NoErrors)
  }
  attr(CompleteModelFit, 'data') <- calibrationData 
  return(CompleteModelFit)
}

