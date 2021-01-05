fitClumpedRegressions<-function(calibrationData, predictionData=NULL,hasMaterial=F, 
                                returnModels=T, n.iter= 50000, burninFrac=0.1,
                                alphaBLM1='dnorm(0,1e-3)', betaBLM1= "dnorm(0,1e-3)",
                                alphaBLM2='dnorm(0,1e-3)', mu0BLM2='dnorm(0,1)', tau0BLM='dunif(1e-1,5)',
                                mu0BLM3='dnorm(0,1)', tau0BLM3='dunif(1e-1,5)', mu1BLM3='dnorm(0,1)', tau1BLM3='dunif(1e-1,5)',
                                useInits=T){
  
  ##Models
  BLM1<-paste("model{
  # Diffuse normal priors for predictors 
  alpha ~ ", alphaBLM1," \n ",
              "beta ~ ", betaBLM1," \n ", 
              "# Gamma prior for scatter 
  tau ~ dunif(0, 100) ##NOTE TO SELF: IT DOESN'T DO THE ADAPTIVE PHASE WITH GAMMA!!
  epsilon <- 1/sqrt(tau)
  # Diffuse normal priors for true x 
  for (i in 1:N){ 
    x[i]~dnorm(0,1e-3) 
    }
  for (i in 1:N){
    obsx[i] ~ dnorm(x[i], pow(errx[i], -2)) 
    obsy[i] ~ dnorm(y[i], pow(erry[i], -2)) 
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + beta*x[i]
  } 
}")
  
  
  ##ANCOVA (Main effects; multiple slopes)
  
  BLM2<-paste("model{
  tau0~",tau0BLM," \n ",
              "mu0~", mu0BLM2," \n ",
              "# Diffuse normal priors for predictors
  alpha ~", alphaBLM2," \n ",
              "for (j in 1:K) {
      beta[j] ~ dnorm(mu0, tau0)
    } 

  # Gamma prior for standard deviation
  tau ~ dgamma(1e-3, 1e-3) # precision
  sigma <- 1 / sqrt(tau) # standard deviation
  # Diffuse normal priors for true x
  for (i in 1:N){
    x1[i]~dnorm(0,1e-3) }
  # Likelihood function
  for (i in 1:N){
    obsy[i] ~ dnorm(y[i],pow(erry[i],-2))
    y[i] ~ dnorm(mu[i],tau)
    obsx1[i] ~ dnorm(x1[i],pow(errx1[i],-2))
    mu[i] <- alpha + beta[type[i]] * x1[i]
  } 
}")
  
  
  
  ##ANCOVA (interaction effects; multiple slopes and intercepts)
  BLM3<-paste("model{
  tau0~", tau0BLM3 ," \n ",
              "mu0~",mu0BLM3," \n ",
              "tau1~",tau1BLM3," \n ",
              "mu1~",mu0BLM3," \n ",
              "# Diffuse normal priors for predictors
    for (i in 1:K) {
      alpha[i] ~ dnorm(mu1, tau1)
      beta[i] ~ dnorm(mu0, tau0)
    } 

  # Gamma prior for standard deviation
  tau ~ dgamma(1e-3, 1e-3) # precision
  sigma <- 1 / sqrt(tau) # standard deviation
  # Diffuse normal priors for true x
  for (i in 1:N){
    x1[i]~dnorm(0,1e-3) }
  # Likelihood function
  for (i in 1:N){
    obsy[i] ~ dnorm(y[i],pow(erry[i],-2))
    y[i] ~ dnorm(mu[i],tau)
    obsx1[i] ~ dnorm(x1[i],pow(errx1[i],-2))
    mu[i] <- alpha[type[i]] + beta[type[i]] * x1[i]
  } 
}")
  
  
  ##Fit linear models
  if(hasMaterial == T){
    
    
    Y= IsoplotR::york(cbind(calibrationData[,c('T2','Temp_Error','D47','D47_SD')]))
    M0=lm(D47 ~ T2, calibrationData)
    M1=lm(D47 ~ T2+Material, calibrationData)
    M2=lm(D47 ~ T2*Material, calibrationData)
    
    
    ##Create the calibrationDatasets for Bayesian Models
    LM_Data <- list(obsx = calibrationData$T2 , obsy = calibrationData$D47 , 
                    errx = calibrationData$Temp_Error, erry = calibrationData$D47_SD, 
                    N=nrow(calibrationData))
    ANCOVA1_Data <- list(obsx1 = calibrationData$T2 , obsy = calibrationData$D47 ,
                         errx1 = calibrationData$Temp_Error, erry = calibrationData$D47_SD,
                         K=length(unique(calibrationData$Material)),
                         N=nrow(calibrationData),
                         type= as.numeric(calibrationData$Material))
    ANCOVA2_Data <- list(obsx1 = calibrationData$T2 , obsy = calibrationData$D47 , 
                         errx1 = calibrationData$Temp_Error, erry = calibrationData$D47_SD, 
                         K=length(unique(calibrationData$Material)),
                         N=nrow(calibrationData),
                         type= as.numeric(calibrationData$Material))
    
    ##Fit the models
    inits <- if(useInits==T){ function () {
      list(alpha = rnorm(1,0,.01),
           beta = rnorm(1,0,.01))
    }}else{NULL}
    
    BLM1_fit <- jags(data = LM_Data,inits = inits,
                     parameters = c("alpha","beta", "epsilon"),
                     model = textConnection(BLM1), n.chains = 3, 
                     n.iter = n.iter, n.burnin = n.iter*burninFrac)
    BLM1_fit <- update(BLM1_fit)
    BLM1_fit <- autojags(BLM1_fit) 
    
    ##ANCOVA 1
    inits <- if(useInits==T){ function () {
      list(alpha = rnorm(1,0,.01),
           beta = rnorm(ANCOVA1_Data$K, 0, 0.01))
    }}else{NULL}
    
    BLM2_fit <- jags(data = ANCOVA1_Data,inits = inits,
                     parameters = c("alpha","beta", "sigma"),model = textConnection(BLM2), n.chains = 3,
                     n.iter = n.iter, n.burnin = n.iter*burninFrac)
    BLM2_fit <- update(BLM2_fit)
    BLM2_fit <- autojags(BLM2_fit) 
    
    ##ANCOVA 2
    inits <- if(useInits==T){ function () {
      list(alpha = rnorm(ANCOVA2_Data$K,0,0.01),
           beta = rnorm(ANCOVA2_Data$K, 0, 0.01))
    }}else{NULL}
    
    BLM3_fit <- jags(data = ANCOVA2_Data,inits = inits,
                     parameters = c("alpha","beta", "sigma"), 
                     model = textConnection(BLM3), n.chains = 3,
                     n.iter = n.iter,  n.burnin = n.iter*burninFrac)
    BLM3_fit <- update(BLM3_fit)
    BLM3_fit <- autojags(BLM3_fit) 
    
    CompleteModelFit<-list('Y'=Y,"M0"=M0,"M1"=M1,"M2"=M2,"BLM1_fit"=BLM1_fit,"BLM2_fit"=BLM2_fit,"BLM3_fit"=BLM3_fit)
  }else{
    Y=IsoplotR::york(calibrationData[,c('T2','Temp_Error','D47','D47_SD')])
    M0=lm(D47 ~ T2, calibrationData)
    LM_Data <- list(obsx = calibrationData$T2 , obsy = calibrationData$D47 , 
                    errx = calibrationData$Temp_Error, erry = calibrationData$D47_SD, 
                    N=nrow(calibrationData))
    ##Fit the models
    inits <- if(useInits==T){ function () {
      list(alpha = rnorm(1,0,.01),
           beta = rnorm(1,0,.01))
    }}else{NULL}
    
    BLM1_fit <- jags(data = LM_Data,inits = inits,
                     parameters = c("alpha","beta", "epsilon"),
                     model = textConnection(BLM1), n.chains = 3, 
                     n.iter = n.iter,  n.burnin = n.iter*burninFrac)
    BLM1_fit <- update(BLM1_fit)
    BLM1_fit <- autojags(BLM1_fit) 
    CompleteModelFit<-list('Y'=Y,'M0'=M0,'BLM1_fit'=BLM1_fit)
  }
  attr(CompleteModelFit, 'data') <- calibrationData 
  return(CompleteModelFit)
}
