#' @param calibrationData A data.frame with at least the four following columns: 'T2','Temp_Error','D47','D47_SD'

fitYorkSingle<-function(calibrationData){
  Y= IsoplotR::york(cbind(calibrationData[,c('T2','Temp_Error','D47','D47_SD')]))
  CompleteModelFit<-list('Y'=Y)
  attr(CompleteModelFit, 'data') <- calibrationData 
  return(CompleteModelFit)
}

#' @param calibrationData A data.frame with at least the four following columns: 'T2','Temp_Error','D47','D47_SD'
#' @param n.iter number of iterations. The detault should be enough for most datasets
#' @param burninFrac Burnin fraction
#' @param alphaBLM1 Prior on the intercept
#' @param betaBLM1 Prior on the slope
#' @param useInits Use initial values in accordance to the prior or not

fitBayesianLinearSingle<-function(calibrationData, 
                                n.iter= 50000, burninFrac=0.1,
                                alphaBLM1='dnorm(0,1e-3)', betaBLM1= "dnorm(0,1e-3)",
                                useInits=T){
  
  
  ##Models
  BLM1<-paste("model{
  # Diffuse normal priors for predictors 
  alpha ~ ", alphaBLM1," \n ",
              "beta ~ ", betaBLM1," \n ", 
              "# Gamma prior for scatter 
  tau ~ dgamma(1e-3, 1e-3) ##NOTE TO SELF: IT DOESN'T DO THE ADAPTIVE PHASE WITH GAMMA!!
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
  
  LM_Data <- list(obsx = calibrationData$T2 , obsy = calibrationData$D47 , 
                  errx = calibrationData$Temp_Error, erry = calibrationData$D47_SD, 
                  N=nrow(calibrationData))
  
  inits <- if(useInits==T){ function () {
    #list(alpha = rnorm(1,0,.01),
    #    beta = rnorm(1,0,.01))
    
    list(alpha = rnorm(1,0.231,0.065),
         beta = rnorm(1,0.039,0.004))
    
  }}else{NULL}
  
  BLM1_fit <- jags(data = LM_Data,inits = inits,
                   parameters = c("alpha","beta", "tau"),
                   model = textConnection(BLM1), n.chains = 3, 
                   n.iter = n.iter, n.burnin = n.iter*burninFrac, n.thin = 10)
  #BLM1_fit <- update(BLM1_fit,n.iter = n.iter,  n.burnin = n.iter*burninFrac)
  BLM1_fit <- autojags(BLM1_fit,n.iter = n.iter,  n.burnin = n.iter*burninFrac, Rhat = 1.1, n.update = 2, n.thin = 4) 
  
  CompleteModelFit<-list("BLM1_fit"=BLM1_fit)
  attr(CompleteModelFit, 'data') <- calibrationData 
  return(CompleteModelFit)
  
}


#' @param calibrationData A data.frame with at least the four following columns: 'T2','Temp_Error','D47','D47_SD'
#' @param n.iter number of iterations. The detault should be enough for most datasets
#' @param burninFrac Burnin fraction
#' @param alphaBLM1 Prior on the intercept
#' @param betaBLM1 Prior on the slope
#' @param useInits Use initial values in accordance to the prior or not


fitBayesianLinearNoErrorsSingle<-function(calibrationData, 
                            n.iter= 50000, burninFrac=0.1,
                            alphaBLM1='dnorm(0,1e-3)', betaBLM1= "dnorm(0,1e-3)",
                            useInits=T){
  
  
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
  
LM_No_error_Data <- list(x = calibrationData$T2 , y = calibrationData$D47,
                         N=nrow(calibrationData))
  
  inits <- if(useInits==T){ function () {
    #list(alpha = rnorm(1,0,.01),
    #    beta = rnorm(1,0,.01))
    
    list(alpha = rnorm(1,0.231,0.065),
         beta = rnorm(1,0.039,0.004))
    
  }}else{NULL}
  
  BLM1_fit_NoErrors <- jags(data = LM_No_error_Data,inits = inits,
                            parameters = c("alpha","beta", "tau"),
                            model = textConnection(BLM1_NoErrors), n.chains = 3,
                            n.iter = n.iter,  n.burnin = n.iter*burninFrac, n.thin = 10)
  
  BLM1_fit_NoErrors <- autojags(BLM1_fit_NoErrors,n.iter = n.iter,  n.burnin = n.iter*burninFrac, Rhat = 1.1, n.update = 2, n.thin = 4)
  
  CompleteModelFit<-list('BLM1_fit_NoErrors'=BLM1_fit_NoErrors)
  attr(CompleteModelFit, 'data') <- calibrationData 
  return(CompleteModelFit)
  
}





#' @param calibrationData A data.frame with at least the four following columns: 'T2','Temp_Error','D47','D47_SD', and 'Material
#' @param n.iter number of iterations. The detault should be enough for most datasets
#' @param burninFrac Burnin fraction
#' @param alphaBLM2 Prior on the intercept
#' @param mu0BLM2 Prior on the mean slopes
#' @param tau0BLM Prior on the standard deviation for the prior for the slopes
#' @param useInits Use initial values in accordance to the prior or not


fitBayesianMainANCOVASimple<-function(calibrationData, 
                                          n.iter= 50000, burninFrac=0.1,
                                          alphaBLM2='dnorm(0,1e-3)', mu0BLM2='dnorm(0,1)', tau0BLM='dunif(1e-1,5)',
                                          useInits=T){
  
  
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
  
  inits <- if(useInits==T){ function () {
    #list(alpha = rnorm(1,0,.01),
    #     beta = rnorm(ANCOVA1_Data$K, 0, 0.01))
    
    list(alpha = rnorm(1,0.231,0.065),
         beta = rnorm(ANCOVA1_Data$K,0.039,0.004))
    
  }}else{NULL}
  
  ANCOVA1_Data <- list(obsx1 = calibrationData$T2 , obsy = calibrationData$D47 ,
                       errx1 = calibrationData$Temp_Error, erry = calibrationData$D47_SD,
                       K=length(unique(calibrationData$Material)),
                       N=nrow(calibrationData),
                       type= as.numeric(calibrationData$Material))
  
  BLM2_fit <- jags(data = ANCOVA1_Data,inits = inits,
                   parameters = c("alpha","beta", "tau"),model = textConnection(BLM2), n.chains = 3,
                   n.iter = n.iter, n.burnin = n.iter*burninFrac, n.thin = 10)
  #BLM2_fit <- update(BLM2_fit,n.iter = n.iter,  n.burnin = n.iter*burninFrac)
  BLM2_fit <- autojags(BLM2_fit,n.iter = n.iter,  n.burnin = n.iter*burninFrac, Rhat = 1.1, n.update = 2, n.thin = 4) 
  
  
  
  
  CompleteModelFit<-list("BLM2_fit"=BLM2_fit)
  attr(CompleteModelFit, 'data') <- calibrationData 
  return(CompleteModelFit)
  
}


#' @param calibrationData A data.frame with at least the four following columns: 'T2','Temp_Error','D47','D47_SD', and 'Material
#' @param n.iter number of iterations. The detault should be enough for most datasets
#' @param burninFrac Burnin fraction
#' @param mu0BLM3 Prior on the mean slopes
#' @param tau0BLM3 Prior on the standard deviation for the prior for the slopes
#' @param mu1BLM3 Prior on the mean intercept
#' @param tau1BLM3 Prior on the standard deviation for the prior for the intercepts
#' @param useInits Use initial values in accordance to the prior or not

fitBayesianInteractionANCOVASimple<-function(calibrationData, 
                                      n.iter= 50000, burninFrac=0.1,
                                      mu0BLM3='dnorm(0,1)', tau0BLM3='dunif(1e-1,5)', mu1BLM3='dnorm(0,1)', tau1BLM3='dunif(1e-1,5)',
                                      useInits=T){
  
  
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
  
  
  inits <- if(useInits==T){ function () {
    #list(alpha = rnorm(1,0,.01),
    #     beta = rnorm(ANCOVA1_Data$K, 0, 0.01))
    
    list(alpha = rnorm(1,0.231,0.065),
         beta = rnorm(ANCOVA1_Data$K,0.039,0.004))
    
  }}else{NULL}
  
  ANCOVA2_Data <- list(obsx1 = calibrationData$T2 , obsy = calibrationData$D47 , 
                       errx1 = calibrationData$Temp_Error, erry = calibrationData$D47_SD, 
                       K=length(unique(calibrationData$Material)),
                       N=nrow(calibrationData),
                       type= as.numeric(calibrationData$Material))
  
  inits <- if(useInits==T){ function () {
    #list(alpha = rnorm(ANCOVA2_Data$K,0,0.01),
    #     beta = rnorm(ANCOVA2_Data$K, 0, 0.01))
    
    list(alpha = rnorm(ANCOVA2_Data$K,0.231,0.065),
         beta = rnorm(ANCOVA2_Data$K,0.039,0.004))
    
  }}else{NULL}
  
  BLM3_fit <- jags(data = ANCOVA2_Data,inits = inits,
                   parameters = c("alpha","beta", "tau"), 
                   model = textConnection(BLM3), n.chains = 3,
                   n.iter = n.iter,  n.burnin = n.iter*burninFrac, n.thin = 10)
  #BLM3_fit <- update(BLM3_fit,n.iter = n.iter,  n.burnin = n.iter*burninFrac)
  BLM3_fit <- autojags(BLM3_fit,n.iter = n.iter,  n.burnin = n.iter*burninFrac, Rhat = 1.1, n.update = 2, n.thin = 4) 
  
  CompleteModelFit<-list("BLM3_fit"=BLM3_fit)
  attr(CompleteModelFit, 'data') <- calibrationData 
  return(CompleteModelFit)
  
}



##Example
require(R2jags)
library(IsoplotR)

Material
dput(Complete_Calibration_List[1:10,c('T2','Temp_Error','D47','D47_SD', 'Material')])



SampleData <-
  structure(
    list(
      T2 = c(
        13.5015404751374,
        13.5015404751374,
        13.5015404751374,
        13.5015404751374,
        13.2087270455852,
        13.2087270455852,
        12.8789929426016,
        12.8789929426016,
        13.5015404751374,
        13.2087270455852
      ),
      Temp_Error = c(2,
                     2, 2, 2, 2, 2, 1, 1, 2, 2),
      D47 = c(
        0.7679185,
        0.76607542,
        0.74889223,
        0.77348394,
        0.77808716,
        0.74983746,
        0.75385961,
        0.76021853,
        0.75197465,
        0.75198252
      ),
      D47_SD = c(
        0.01041691,
        0.00659912,
        0.02020052,
        0.00278529,
        0.00535896,
        0.0111864,
        0.0024896,
        0.01489061,
        0.01053873,
        0.01275677
      ),
      Material = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    ),
    row.names = c(NA,
                  10L),
    class = "data.frame"
  )

fitYorkSingle(SampleData)[[1]]
fitBayesianLinearSingle(SampleData)[[1]]
fitBayesianLinearNoErrorsSingle(SampleData)[[1]]
fitBayesianMainANCOVASimple(SampleData)[[1]]
fitBayesianInteractionANCOVASimple(SampleData)[[1]]






