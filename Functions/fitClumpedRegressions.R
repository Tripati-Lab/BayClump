fitClumpedRegressions<-function(calibrationData, predictionData=NULL,hasMaterial=F, 
                                returnModels=T, n.iter= 50000, burninFrac=0.1,
                                alphaBLM1='dnorm(0.231,0.065)', betaBLM1= "dnorm(0.039,0.004)",
                                useInits=T, D47error="D47error"){
  
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
}")
  
  
  LM_No_error_Data <- list(x = calibrationData$Temperature , y = calibrationData$D47,
                           N=nrow(calibrationData))
  
  ##Fit linear models
  if(hasMaterial == T){
    
    
    Y= NULL#IsoplotR::york(cbind(calibrationData[,c('Temperature','TempError','D47', 'D47error')]))
    M0=NULL#lm(D47 ~ Temperature, calibrationData)
    M1=NULL#lm(D47 ~ Temperature+Material, calibrationData)
    M2=NULL#lm(D47 ~ Temperature*Material, calibrationData)
    
    
    ##Create the calibrationDatasets for Bayesian Models
    LM_Data <- list(obsx = calibrationData$Temperature , obsy = calibrationData$D47 , 
                    errx = calibrationData$TempError, erry = calibrationData[,D47error], 
                    N=nrow(calibrationData))
    
    ANCOVA2_Data <- list(obsx1 = calibrationData$Temperature , obsy = calibrationData$D47 , 
                         errx1 = calibrationData$TempError, erry = calibrationData[,D47error], 
                         K=length(unique(calibrationData$Material)),
                         N=nrow(calibrationData),
                         type= as.numeric(calibrationData$Material))
    
    ##Fit the models
    inits <- if(useInits==T){ function () {
      #list(alpha = rnorm(1,0,.01),
      #    beta = rnorm(1,0,.01))
      
      list(alpha = rnorm(1,0.231,0.065),
           beta = rnorm(1,0.039,0.004))
      
    }}else{NULL}
    
    BLM1_fit <- jags(data = LM_Data,inits = inits,
                     parameters = c("alpha","beta", "tauy"),
                     model = textConnection(BLM1), n.chains = 3, 
                     n.iter = n.iter, n.burnin = n.iter*burninFrac, n.thin = 10)
    #BLM1_fit <- update(BLM1_fit,n.iter = n.iter,  n.burnin = n.iter*burninFrac)
    BLM1_fit <- autojags(BLM1_fit,n.iter = n.iter,  n.burnin = n.iter*burninFrac, Rhat = 1.1, n.update = 2, n.thin = 4) 
    
    BLM1_fit_NoErrors <- jags(data = LM_No_error_Data,inits = inits,
                              parameters = c("alpha","beta", "tau"),
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
                     parameters = c("alpha","beta"), 
                     model = textConnection(BLM3), n.chains = 3,
                     n.iter = n.iter,  n.burnin = n.iter*burninFrac, n.thin = 10)
    #BLM3_fit <- update(BLM3_fit,n.iter = n.iter,  n.burnin = n.iter*burninFrac)
    BLM3_fit <- autojags(BLM3_fit,n.iter = n.iter,  n.burnin = n.iter*burninFrac, Rhat = 1.1, n.update = 2, n.thin = 4) 
    
    
    CompleteModelFit<-list('Y'=Y,"M0"=M0,"M1"=M1,"M2"=M2,"BLM1_fit"=BLM1_fit,'BLM1_fit_NoErrors'=BLM1_fit_NoErrors, "BLM3_fit"=BLM3_fit)
  }else{
    Y=NULL#IsoplotR::york(calibrationData[,c('Temperature','TempError','D47','D47error')])
    M0=NULL#lm(D47 ~ Temperature, calibrationData)
    LM_Data <- list(obsx = calibrationData$Temperature , obsy = calibrationData$D47 , 
                    errx = calibrationData$TempError, erry = calibrationData[,D47error], 
                    N=nrow(calibrationData))
    ##Fit the models
    inits <- if(useInits==T){ function () {
      #list(alpha = rnorm(1,0,.01),
      #     beta = rnorm(1,0,.01))
      list(alpha = rnorm(1,0.231,0.065),
           beta = rnorm(1,0.039,0.004))
    }}else{NULL}
    
    BLM1_fit <- jags(data = LM_Data,inits = inits,
                     parameters = c("alpha","beta", "tauy"),
                     model = textConnection(BLM1), n.chains = 3, 
                     n.iter = n.iter,  n.burnin = n.iter*burninFrac, n.thin = 10)
    #BLM1_fit <- update(BLM1_fit,n.iter = n.iter,  n.burnin = n.iter*burninFrac)
    BLM1_fit <- autojags(BLM1_fit,n.iter = n.iter,  n.burnin = n.iter*burninFrac, Rhat = 1.1, n.update = 2, n.thin = 4) 
    
    BLM1_fit_NoErrors <- jags(data = LM_No_error_Data,inits = inits,
                              parameters = c("alpha","beta", "tau"),
                              model = textConnection(BLM1_NoErrors), n.chains = 3,
                              n.iter = n.iter,  n.burnin = n.iter*burninFrac, n.thin = 10)
    
    BLM1_fit_NoErrors <- autojags(BLM1_fit_NoErrors,n.iter = n.iter,  n.burnin = n.iter*burninFrac, Rhat = 1.1, n.update = 2, n.thin = 4)
    
    CompleteModelFit<-list('Y'=Y,'M0'=M0,'BLM1_fit'=BLM1_fit,'BLM1_fit_NoErrors'=BLM1_fit_NoErrors)
  }
  attr(CompleteModelFit, 'data') <- calibrationData 
  return(CompleteModelFit)
}

  
