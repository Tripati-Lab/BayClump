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


BayesianPredictions <- function(calibrationData, 
                                n.iter = 5000, 
                                priors = "Informative",
                                samples=NULL,
                                init.values = FALSE, 
                                D47Pred,
                                materialsPred){
  
  
  
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
  
  
  BLM1 <- paste("model{
                # Diffuse normal priors for predictors
                alpha ~ ", alphaBLM1," \n ",
                "beta ~ ", betaBLM1," \n ",
                "
    tau <- pow(sigma, -2) 
    sigma ~ dunif(0, 100)                             
    taux ~ dunif(40, 2700) #based on SD used for simulations
    
  # calibration
  for(i in 1:N){   
    truex[i] ~ dnorm(11,0.01)
    x[i] ~ dnorm(truex[i], taux)
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + beta * truex[i]
  }
  
  # Inversion
  for(i in 1:NPred){
  yTemp[i] <- sqrt((beta * 10 ^ 6) / (yPred[i] - alpha)) - 273.15
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
  
  # Inversion
  for(i in 1:NPred){
  yTemp[i] <- sqrt((beta * 10 ^ 6) / (yPred[i] - alpha)) - 273.15
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
    taux ~ dunif(40, 2700) #based on SD used for simulations


    # Likelihood function
    for (i in 1:N){
    truex[i] ~ dnorm(11,0.01)
    x[i] ~ dnorm(truex[i], taux)
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha[type[i]] + beta[type[i]] * truex[i]
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
     yPredFixed[j] <-  sum(alpha[type[j]] + beta[type[j]] * truex[j])
  }

  varFixed <- (sd(yPredFixed))^2
  varResidual <- sigma^2 # get the variance of residuals
  varRandom <- subjEffSigma^2  # get the variance of random plot effect
  # calculate marginal R^2
  marginalR2 <- varFixed / (varFixed + varRandom + varResidual)
  # calculate conditional R^2
  conditionalR2 <- (varRandom + varFixed) / (varFixed + varRandom + varResidual) 

  # Inversion
  for(i in 1:NPred){
  yTemp[i] <- sqrt((beta[typePred[i]] * 10 ^ 6) / (yPred[i] - alpha[typePred[i]])) - 273.15
  }

}")
  
  #Data
  
  LMM_Data <- list(x = calibrationData$Temperature, 
                   y = calibrationData$D47, 
                   K = length(unique(calibrationData$Material)),
                   N = nrow(calibrationData),
                   type = as.numeric(calibrationData$Material),
                   NPred = length(D47Pred),
                   yPred = D47Pred,
                   typePred = materialsPred)
  
  LM_Data <- list(x = calibrationData$Temperature, 
                  y = calibrationData$D47,
                  N = nrow(calibrationData),
                  NPred = length(D47Pred),
                  yPred = D47Pred)
  
  
  #Inits
  initsMixed <- function () {
    list(alpha = rnorm(LMM_Data$K,0.231,0.065),
         beta = rnorm(LMM_Data$K,0.039,0.004))
    
  }
  
  initsSimple <- function () {
    list(alpha = rnorm(1,0.231,0.065),
         beta = rnorm(1,0.039,0.004))
    
  }
  
  #Fit models
  BLM3_fit <- jags(data = LMM_Data, 
                   inits = if(init.values){initsMixed}else{NULL},
                   parameters = c("yTemp"), 
                   model = textConnection(BLM3), 
                   n.chains = 3,
                   n.iter = n.iter)
  BLM3_fit <- autojags(BLM3_fit)
  
  BLM1_fit <- jags(data = LM_Data, 
                   inits = if(init.values){initsSimple}else{NULL},
                   parameters = c("yTemp"),
                   model = textConnection(BLM1), 
                   n.chains = 3, 
                   n.iter = n.iter)
  BLM1_fit <- autojags(BLM1_fit)
  
  BLM1_fit_NoErrors <- jags(data = LM_Data, 
                            inits = if(init.values){initsSimple}else{NULL},
                            parameters = c("yTemp"),
                            model = textConnection(BLM1_NoErrors), 
                            n.chains = 3,
                            n.iter = n.iter)
  BLM1_fit_NoErrors <- autojags(BLM1_fit_NoErrors)
  
  dsSumBLM3 <- summary(as.mcmc(BLM3_fit))$quantiles[-1,]
  dsSumBLM1 <- summary(as.mcmc(BLM1_fit))$quantiles[-1,]
  dsSumBLM1_NE <- summary(as.mcmc(BLM1_fit_NoErrors))$quantiles[-1,]
  
  postPredBLMM <-     cbind.data.frame(D47Pred, 
                                       Material = materialsPred,
                                       Tc = BLM3_fit$BUGSoutput$mean$yTemp,
                                       sd = BLM3_fit$BUGSoutput$sd$yTemp)
  
  
  postPredBLM2 <-     cbind.data.frame(D47Pred, 
                                       Tc = BLM1_fit_NoErrors$BUGSoutput$mean$yTemp,
                                       sd = BLM1_fit_NoErrors$BUGSoutput$sd$yTemp)
  
  
  postPredBLM1 <-     cbind.data.frame(D47Pred, 
                                       Tc = BLM1_fit$BUGSoutput$mean$yTemp,
                                       sd = BLM1_fit$BUGSoutput$sd$yTemp)
  
  
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


