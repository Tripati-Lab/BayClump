#" This function generate temperature predictions (in 10^6/T2) based on a 
#" calibration dataset and target D47. Note that this alternative function
#" propagates uncertainty around the target D47.
#" 
#" @param calibrationData The calibration dataset
#" @param hasMaterial Whether only a mixed model should be run
#" @param n.iter number of MCMC iterations
#" @param burninFrac burnin fraction (0-1)
#" @param D47Pred the target D47
#" @param D47Prederror error in the target D47
#" @param materialsPred Material of the target D47
#" @param priors Informative priors or not on the slope and intercept

BayesianPredictions <- function(bayeslincals, 
                                    D47Pred,
                                    D47Prederror,
                                    materialsPred,
                                    nsamp = NULL){
  

BLM1<-paste("model{
  for(i in 1:N){ 
    x[i] ~ dnorm(11, 0.0001)
    y[i] ~ dnorm(mu2[i], tau)
    x2[i] <- sqrt((beta * 10^6) / (y[i] - alpha)) - 273.15
    mu2[i] <- alpha  + beta * x[i]
    y2[i] ~ dnorm(y[i], pow(y2err[i],-2))
  }
}")
  
  BLMM<-paste("model{
    for(i in 1:N){ 
    x[i] ~ dnorm(11, 0.0001)
    y[i] ~ dnorm(mu2[i], tau)
    x2[i] <- sqrt((beta * 10^6) / (y[i] - alpha)) - 273.15
    mu2[i] <- alpha[typePred[i]]  + beta[typePred[i]] * x[i]
    y2[i] ~ dnorm(y[i], pow(y2err[i],-2))
  }
}")
  
  
  postBLM<- do.call(rbind, as.mcmc(bayeslincals$BLM1_fit_NoErrors))
  nsamp <- if(is.null(nsamp)){nrow(postBLM)}else{nsamp}
  
  postPredBLM1 <- do.call(rbind,lapply(sample(1:nrow(postBLM), nsamp), function(i){
    LM_No_error_Data <- list(
      N=length(D47Pred),
      y2=D47Pred,
      y2err=D47Prederror,
      alpha=postBLM[i,'alpha'],
      beta=postBLM[i,'beta'],
      tau=postBLM[i,'tau']
    )
    
    BLM1_fit_NoErrors <- jags(data = LM_No_error_Data,#inits = inits,
                              parameters = c("x2"),
                              model = textConnection(BLM1), n.chains = 3,
                              n.iter = 20000)
    unlist(BLM1_fit_NoErrors$BUGSoutput$mean[-1])
  }))
  
  postBLM<- do.call(rbind, as.mcmc(bayeslincals$BLM1_fit))
  postPredBLM2 <- do.call(rbind,lapply(sample(1:nrow(postBLM), nsamp), function(i){
    LM_No_error_Data <- list(
      N=length(D47Pred),
      y2=D47Pred,
      y2err=D47Prederror,
      alpha=postBLM[i,'alpha'],
      beta=postBLM[i,'beta'],
      tau=postBLM[i,'tau']
    )
    
    BLM1_fit_NoErrors <- jags(data = LM_No_error_Data,#inits = inits,
                              parameters = c("x2"),
                              model = textConnection(BLM1), n.chains = 3,
                              n.iter = 20000)
    unlist(BLM1_fit_NoErrors$BUGSoutput$mean[-1])
  }))
  
  postBLMM<- do.call(rbind, as.mcmc(bayeslincals$BLM3_fit))
  postPredBLMM <- do.call(rbind,lapply(sample(1:nrow(postBLMM), nsamp), function(i){
    
    LM_No_error_Data <- list(
      N=length(D47Pred),
      y2=D47Pred,
      y2err=D47Prederror,
      alpha=lapply(grep("alpha", colnames(postBLMM)), function(x) postBLMM[i,x]),
      beta=lapply(grep("beta", colnames(postBLMM)), function(x) postBLMM[i,x]),
      tau=postBLMM[i,'tau'],
      typePred=materialsPred
    )
    
    
    BLM1_fit_NoErrors <- jags(data = LM_No_error_Data,#inits = inits,
                              parameters = c("x2"),
                              model = textConnection(BLMM), n.chains = 3,
                              n.iter = 20000)
    unlist(BLM1_fit_NoErrors$BUGSoutput$mean[-1])
  }))
  
  CompleteModelFit<-list("BLM1_fit"=postPredBLM1,"BLM1_fit_NoErrors"=postPredBLM2, "BLM3"=postPredBLMM)
  
  return(CompleteModelFit)
}

#" Generate a dataset reflecting the priors used to run the analyses
#" 
#" @param prior Informative or not
#" @param n number of observations to simulate


generatePriorReconstructions <- function(prior, n=1000){
  
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
  x <- rnorm(n, 11,0.0001)
  xC <- sqrt(10^6/x)-273.15
  
  data <- cbind.data.frame(alpha=rnorm(n, params[1,2], params[1,3]), 
                           beta=rnorm(n, params[2,2], params[2,3]),
                           x,
                           xC)
  attr(data, "priors") <-  prior
  attr(data, "params") <-  params
  data
}


