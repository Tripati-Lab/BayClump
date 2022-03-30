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
                                    nsamp=500){
  

BLM1<-paste("model{
  for(i in 1:N){ 
    x[i] ~ dnorm(11, 0.394)
    y[i] ~ dnorm(mu2[i], tau)
    x2[i] <- sqrt((beta * 10^6) / (y[i] - alpha)) - 273.15
    mu2[i] <- alpha  + beta * x[i]
    y2[i] ~ dnorm(y[i], pow(y2err[i],-2))
  }
}")
  
  postBLM<- bayeslincals$BLM1_fit_NoErrors$BUGSoutput$sims.matrix
  postPredBLM1 <- do.call(rbind,lapply(sample(1:nrow(postBLM), nsamp), function(j){
    tryCatch({
      LM_No_error_Data <- list(
        N=length(D47Pred),
        y2=D47Pred,
        y2err=D47Prederror,
        alpha=postBLM[j,'alpha'],
        beta=postBLM[j,'beta'],
        tau=postBLM[j,'tau']
      )
    
    BLM1_fit_NoErrors <- jags(data = LM_No_error_Data,
                              parameters = c("x2"),
                              model = textConnection(BLM1), n.chains = 3,
                              n.iter =  20000)
    
    cbind.data.frame(D47Pred, 
                     D47Prederror, 
                     BLM1_fit_NoErrors$BUGSoutput$mean[-1],
                     BLM1_fit_NoErrors$BUGSoutput$sd[-1]
                     )
    
  }, error=function(e){c(NA,NA)})
  }))
  postPredBLM1 <-aggregate(postPredBLM1[, 3:4], list(postPredBLM1$D47Pred, postPredBLM1$D47Prederror), mean)
  
  
  postBLM<- bayeslincals$BLM1_fit$BUGSoutput$sims.matrix
  postPredBLM2 <- do.call(rbind,lapply(sample(1:nrow(postBLM), nsamp), function(j){
    tryCatch({
      LM_No_error_Data <- list(
        N=length(D47Pred),
        y2=D47Pred,
        y2err=D47Prederror,
        alpha=postBLM[j,'alpha'],
        beta=postBLM[j,'beta'],
        tau=postBLM[j,'tau']
      )
      
      BLM1_fit_NoErrors <- jags(data = LM_No_error_Data,
                                parameters = c("x2"),
                                model = textConnection(BLM1), n.chains = 3,
                                n.iter =  20000)
      
      cbind.data.frame(D47Pred, 
                       D47Prederror, 
                       BLM1_fit_NoErrors$BUGSoutput$mean[-1],
                       BLM1_fit_NoErrors$BUGSoutput$sd[-1]
      )
      
    }, error=function(e){c(NA,NA)})
  }))
  postPredBLM2 <-aggregate(postPredBLM2[, 3:4], list(postPredBLM2$D47Pred, postPredBLM2$D47Prederror), mean)
  
  postBLMM<- bayeslincals$BLM3_fit$BUGSoutput$sims.matrix
  postPredBLMM <- do.call(rbind,lapply(sample(1:nrow(postBLM), nsamp), function(i){
  tryCatch({
    
    alphas=lapply(grep("alpha", colnames(postBLMM)), function(x) postBLMM[,x])
    betas=lapply(grep("beta", colnames(postBLMM)), function(x) postBLMM[,x])
    
      LM_No_error_Data <- list(
        N=length(D47Pred),
        y2=D47Pred,
        y2err=D47Prederror,
        alpha= alphas[[1]][j], #Currently working for a single material
        beta= betas[[1]][j], #Currently working for a single material
        tau=postBLMM[j,'tau']
      )
    
      BLM1_fit_NoErrors <- jags(data = LM_No_error_Data,#inits = inits,
                                parameters = c("x2"),
                                model = textConnection(BLM1), n.chains = 3,
                                n.iter =  10000)
      
      cbind.data.frame(D47Pred, 
                       D47Prederror, 
                       BLM1_fit_NoErrors$BUGSoutput$mean[-1],
                       BLM1_fit_NoErrors$BUGSoutput$sd[-1]
      )
      
      }, error=function(e){c(NA,NA)})
  }))
  postPredBLMM <-aggregate(postPredBLMM[, 3:4], list(postPredBLMM$D47Pred, postPredBLMM$D47Prederror), mean)

  
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


