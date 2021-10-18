fitClumpedRegressions<-function(calibrationData, predictionData=NULL, hasMaterial=F, 
                                n.iter= 5000, burninFrac=0.5,
                                alphaBLM1='dnorm(0.231,0.065)', betaBLM1= 'dnorm(0.039,0.004)',
                                useInits=T, D47error='D47error'){
  
  ##Models
  BLM1<-paste(" model{
    # Diffuse normal priors for predictors
    alpha ~ ", alphaBLM1," \n ",
              "beta ~ ", betaBLM1," \n ", 
    "
    sigma <- 1/sqrt(tauy)                              
    tauy ~ dgamma(0.1, 0.1)                                
    
    for (i in 1:N){
        x[i] ~ dnorm(11,0.01)
    }
    # Likelihood
    for (i in 1:N){
        obsy[i] ~ dnorm(y[i],pow(erry[i],-2))
        y[i] ~ dnorm(mu[i],tauy)
        obsx[i] ~ dnorm(x[i],pow(errx[i],-2))
        mu[i] <- alpha + beta*x[i]
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
      list(alpha = rnorm(1,0.231,0.065),
           beta = rnorm(1,0.039,0.004))
      
    }}else{NULL}
    
    BLM1_fit <- jags(data = LM_Data, inits = inits,
                     parameters = c("alpha","beta", "tauy"),
                     model = textConnection(BLM1), n.chains = 3, 
                     n.iter = n.iter, n.burnin = n.iter*burninFrac)

    BLM1_fit_NoErrors <- jags(data = LM_No_error_Data,inits = inits,
                              parameters = c("alpha","beta", "tau"),
                              model = textConnection(BLM1_NoErrors), n.chains = 3,
                              n.iter = n.iter,  n.burnin = n.iter*burninFrac)
    
    
    ##ANCOVA 2
    inits <- if(useInits==T){ function () {
      list(alpha = rnorm(ANCOVA2_Data$K,0.231,0.065),
           beta = rnorm(ANCOVA2_Data$K,0.039,0.004))
      
    }}else{NULL}
    
    BLM3_fit <- jags(data = ANCOVA2_Data,inits = inits,
                     parameters = c("alpha","beta","conditionalR2", "marginalR2"), 
                     model = textConnection(BLM3), n.chains = 3,
                     n.iter = n.iter,  n.burnin = n.iter*burninFrac)

    R2sComplete<-rbind.data.frame(getR2Bayesian(BLM1_fit, calibrationData=calibrationData, hasMaterial = F),
    getR2Bayesian(BLM1_fit_NoErrors, calibrationData=calibrationData, hasMaterial = F))
    
    BLM3_fitR2<-t(as.data.frame(BLM3_fit$BUGSoutput$summary[grep('conditional',row.names(BLM3_fit$BUGSoutput$summary)),c(5,3,7)]))
    colnames(BLM3_fitR2)<-names(R2sComplete)
    R2sComplete<-rbind(R2sComplete,BLM3_fitR2)
    row.names(R2sComplete)<-NULL
    BLM3_fitR2M<-t(as.data.frame(BLM3_fit$BUGSoutput$summary[grep('marginal',row.names(BLM3_fit$BUGSoutput$summary)),c(5,3,7)]))
    colnames(BLM3_fitR2M)<-names(R2sComplete)
    row.names(BLM3_fitR2M)<-NULL
    R2sComplete<-rbind(R2sComplete,BLM3_fitR2M)
    R2sComplete$model<-c("BLM1_fit", "BLM1_fit_NoErrors", "BLM3_fit","BLM3_fit")
    R2sComplete$class<-c("Conditional", "Conditional", "Conditional","Marginal")
    
    
    DICs<-c(BLM1_fit$BUGSoutput$DIC, BLM1_fit_NoErrors$BUGSoutput$DIC,  BLM3_fit$BUGSoutput$DIC)
    names(DICs)<-c("BLM1_fit", "BLM1_fit_NoErrors", "BLM3_fit")
    
    CompleteModelFit<-list('Y'=Y,"M0"=M0,"M1"=M1,"M2"=M2,"BLM1_fit"=BLM1_fit,'BLM1_fit_NoErrors'=BLM1_fit_NoErrors, "BLM3_fit"=BLM3_fit)
  }else{
    Y=NULL#IsoplotR::york(calibrationData[,c('Temperature','TempError','D47','D47error')])
    M0=NULL#lm(D47 ~ Temperature, calibrationData)
    LM_Data <- list(obsx = calibrationData$Temperature , obsy = calibrationData$D47 , 
                    errx = calibrationData$TempError, erry = calibrationData[,D47error], 
                    N=nrow(calibrationData))
    ##Fit the models
    inits <- if(useInits==T){ function () {
      list(alpha = rnorm(1,0.231,0.065),
           beta = rnorm(1,0.039,0.004))
    }}else{NULL}
    
    BLM1_fit <- jags(data = LM_Data,inits = inits,
                     parameters = c("alpha","beta", "tauy"),
                     model = textConnection(BLM1), n.chains = 3, 
                     n.iter = n.iter,  n.burnin = n.iter*burninFrac)

    BLM1_fit_NoErrors <- jags(data = LM_No_error_Data,inits = inits,
                              parameters = c("alpha","beta", "tau"),
                              model = textConnection(BLM1_NoErrors), n.chains = 3,
                              n.iter = n.iter,  n.burnin = n.iter*burninFrac)
    
    R2sComplete<-rbind.data.frame(getR2Bayesian(BLM1_fit, calibrationData=calibrationData),
                       getR2Bayesian(BLM1_fit_NoErrors, calibrationData=calibrationData))
    R2sComplete$model<-c("BLM1_fit", "BLM1_fit_NoErrors")
    DICs<-c(BLM1_fit$BUGSoutput$DIC, BLM1_fit_NoErrors$BUGSoutput$DIC)
    names(DICs)<-c("BLM1_fit", "BLM1_fit_NoErrors")
    
    
    CompleteModelFit<-list('Y'=Y,'M0'=M0,'BLM1_fit'=BLM1_fit,'BLM1_fit_NoErrors'=BLM1_fit_NoErrors)
  }
  attr(CompleteModelFit, 'data') <- calibrationData 
  attr(CompleteModelFit, 'R2s') <- R2sComplete 
  attr(CompleteModelFit, 'DICs') <- DICs 
  
  return(CompleteModelFit)
}

  
