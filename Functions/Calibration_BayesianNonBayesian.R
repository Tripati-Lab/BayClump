#" This function generate temperature predictions (in 10^6/T2) based on a 
#" calibration dataset and target D47. Note that this alternative function
#" propagates uncertainty around the target D47.
#" 
#" @param calibrationData The calibration dataset
#" @param hasMaterial Whether only a mixed model should be run
#" @param n.iter number of MCMC iterations
#" @param burninFrac burnin fraction (0-1)
#" @param priors Informative priors or not on the slope and intercept
#" @param useInits whether we should use inits to the run
#" @param D47error the column in calibrationData containing the uncertainty in D47


fitClumpedRegressions <<- function(calibrationData, 
                                hasMaterial = FALSE, 
                                n.iter = 5000, 
                                burninFrac = 0.5,
                                priors = "Informative",
                                useInits = TRUE, 
                                D47error = "D47error"){
  
  
  if(priors == "Informative"){
    alphaBLM1 = "dnorm(0.231,0.065)" 
    betaBLM1 = "dnorm(0.039,0.004)"}else{
      alphaBLM1 = "dnorm(0, 0.01)" 
      betaBLM1 = "dnorm(0, 0.01)"
    }
  
  
  ##Models
  BLM1<-paste(" model{
    # Diffuse normal priors for predictors
    alpha ~ ", alphaBLM1," \n ",
              "beta ~ ", betaBLM1," \n ", 
              "
    sigma <- 1/sqrt(tau)                              
    tau ~ dgamma(0.1, 0.1)                                
    
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
    
    ##Log-likelihood
  for(i in 1:N){ 
   regression_residual[i] <- y[i] - mu[i]
   zloglik[i] <- logdensity.norm(y[i], mu[i], tau)
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
  
  ##Log-likelihood
  for(i in 1:N){ 
   regression_residual[i] <- y[i] - mu[i]
   zloglik[i] <- logdensity.norm(y[i], mu[i], tau)
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

##Log-likelihood
  for(i in 1:N){ 
   regression_residual[i] <- y[i] - mu[i]
   zloglik[i] <- logdensity.norm(y[i], mu[i], tau)
  }

}")
  
  
  
  
  ##If mixed
  if(hasMaterial == T){
    
    Y= NULL#IsoplotR::york(cbind(calibrationData[,c("Temperature","TempError","D47", "D47error")]))
    M0=NULL#lm(D47 ~ Temperature, calibrationData)
    M1=NULL#lm(D47 ~ Temperature+Material, calibrationData)
    M2=NULL#lm(D47 ~ Temperature*Material, calibrationData)
    
    ANCOVA2_Data <- list(obsx1 = calibrationData$Temperature , obsy = calibrationData$D47 , 
                         errx1 = abs(calibrationData$TempError), erry = calibrationData[,D47error], 
                         K=length(unique(calibrationData$Material)),
                         N=nrow(calibrationData),
                         type= as.numeric(calibrationData$Material))
    
    inits <- if(useInits==T){ function () {
      list(alpha = rnorm(ANCOVA2_Data$K,0.231,0.065),
           beta = rnorm(ANCOVA2_Data$K,0.039,0.004))
      
    }}else{NULL}
    
    BLM3_fit <- jags(data = ANCOVA2_Data, #inits = inits,
                     parameters = c("alpha","beta","conditionalR2", "marginalR2","zloglik"), 
                     model = textConnection(BLM3), n.chains = 3,
                     n.iter = n.iter,  n.burnin = n.iter*burninFrac)
    
    tmatrix <- BLM3_fit$BUGSoutput$sims.matrix
    tmatrix <- tmatrix[grep("zloglik", colnames(tmatrix)),]
    aM<-waic(tmatrix)
    
    #Avoid running the other models when running the mixed model
    BLM1_fit<- BLM3_fit
    BLM1_fit_NoErrors <- BLM3_fit
    
    R2sComplete<-rbind.data.frame(getR2Bayesian(BLM1_fit, calibrationData=calibrationData, hasMaterial = T),
                                  getR2Bayesian(BLM1_fit_NoErrors, calibrationData=calibrationData, hasMaterial = T))
    
    BLM3_fitR2<-t(as.data.frame(BLM3_fit$BUGSoutput$summary[grep("conditional",row.names(BLM3_fit$BUGSoutput$summary)),c(5,3,7)]))
    colnames(BLM3_fitR2)<-names(R2sComplete)
    R2sComplete<-rbind(R2sComplete,BLM3_fitR2)
    row.names(R2sComplete)<-NULL
    BLM3_fitR2M<-t(as.data.frame(BLM3_fit$BUGSoutput$summary[grep("marginal",row.names(BLM3_fit$BUGSoutput$summary)),c(5,3,7)]))
    colnames(BLM3_fitR2M)<-names(R2sComplete)
    row.names(BLM3_fitR2M)<-NULL
    R2sComplete<-rbind(R2sComplete,BLM3_fitR2M)
    R2sComplete$model<-c("BLM1_fit", "BLM1_fit_NoErrors", "BLM3_fit","BLM3_fit")
    R2sComplete$class<-c("Conditional", "Conditional", "Conditional","Marginal")
    
    DICs<-c(aM$estimates[3,1], aM$estimates[3,1],  aM$estimates[3,1])
    names(DICs)<-c("BLM1_fit", "BLM1_fit_NoErrors", "BLM3_fit")
    
    CompleteModelFit<-list("Y"=Y,"M0"=M0,"M1"=M1,"M2"=M2,"BLM1_fit"=BLM1_fit,"BLM1_fit_NoErrors"=BLM1_fit_NoErrors, "BLM3_fit"=BLM3_fit)
  }else{
    
    
    ##Non-mixed models
    
    LM_No_error_Data <- list(x = calibrationData$Temperature , y = calibrationData$D47,
                             N=nrow(calibrationData))
    
    ##Create the calibrationDatasets for Bayesian Models
    LM_Data <- list(obsx = calibrationData$Temperature , obsy = calibrationData$D47 , 
                    errx = abs(calibrationData$TempError), erry = calibrationData[,D47error], 
                    N=nrow(calibrationData))
    
    
    ##Fit the models
    inits <- if(useInits==T){ function () {
      list(alpha = rnorm(1,0.231,0.065),
           beta = rnorm(1,0.039,0.004))
      
    }}else{NULL}
    
    BLM1_fit <- jags(data = LM_Data, #inits = inits,
                     parameters = c("alpha","beta", "tau","zloglik"),
                     model = textConnection(BLM1), n.chains = 3, 
                     n.iter = n.iter, n.burnin = n.iter*burninFrac)
    
    BLM1_fit_NoErrors <- jags(data = LM_No_error_Data,#inits = inits,
                              parameters = c("alpha","beta", "tau","zloglik"),
                              model = textConnection(BLM1_NoErrors), n.chains = 3,
                              n.iter = n.iter,  n.burnin = n.iter*burninFrac)
    
    tmatrix <- BLM1_fit$BUGSoutput$sims.matrix
    tmatrix <- tmatrix[grep("zloglik", colnames(tmatrix)),]
    aMErrors <-waic(tmatrix)
    
    tmatrix <- BLM1_fit_NoErrors$BUGSoutput$sims.matrix
    tmatrix <- tmatrix[grep("zloglik", colnames(tmatrix)),]
    aMNoErrors<-waic(tmatrix)
    
    Y=NULL#IsoplotR::york(calibrationData[,c("Temperature","TempError","D47","D47error")])
    M0=NULL#lm(D47 ~ Temperature, calibrationData)
    
    R2sComplete<-rbind.data.frame(getR2Bayesian(BLM1_fit, calibrationData=calibrationData),
                                  getR2Bayesian(BLM1_fit_NoErrors, calibrationData=calibrationData))
    R2sComplete$model<-c("BLM1_fit", "BLM1_fit_NoErrors")
    DICs<-c(aMErrors$estimates[3,1], aMNoErrors$estimates[3,1])
    names(DICs)<-c("BLM1_fit", "BLM1_fit_NoErrors")
    
    CompleteModelFit<-list("Y"=Y,"M0"=M0,"BLM1_fit"=BLM1_fit,"BLM1_fit_NoErrors"=BLM1_fit_NoErrors)
  }
  
  attr(CompleteModelFit, "data") <- calibrationData 
  attr(CompleteModelFit, "R2s") <- R2sComplete 
  attr(CompleteModelFit, "DICs") <- DICs 
  
  return(CompleteModelFit)
}


#" Bootstrap York regression models from a calibration dataset
#" 
#" @param data The calibration dataset
#" @param replicates Number of bootstrap replicates
#" @param samples Number of samples per replicate
#" @param D47error The column in data containing the errors in D47

simulateYork_measured <<- function(data, 
                                   replicates, 
                                   samples = NULL, 
                                   D47error = "D47error"){
  do.call(rbind,lapply(1:replicates, function(x){
    dataSub <- data[sample(seq_along(data[,1]), if(is.null(samples)){nrow(data)}else{samples}, replace = T),]
    dataSub$y_SE <- dataSub[,D47error]
    dataSub$x_SE <- abs(dataSub$TempError)
    Reg <- york(cbind.data.frame(dataSub$Temperature, dataSub$x_SE, dataSub$D47, dataSub$y_SE))
    cbind.data.frame("intercept"=Reg$a[1],"slope"=Reg$b[1])
  }))
}


#" Bootstrap an OLS regression models from a calibration dataset
#" 
#" @param data The calibration dataset
#" @param replicates Number of bootstrap replicates
#" @param samples Number of samples per replicate
#" @param D47error The column in data containing the errors in D47

simulateLM_measured <<- function(data, 
                               replicates, 
                               samples = NULL, 
                               D47error="D47error"){
  
  a<-lapply(1:replicates, function(x){
    dataSub<-data[sample(seq_along(data[,1]), if(is.null(samples)){nrow(data)}else{samples}, replace = T),]
    dataSub$y_SE<-dataSub[,D47error]
    dataSub$x_SE<-dataSub$TempError
    Reg<-summary(lm(D47 ~ Temperature,  dataSub))
    res<-cbind.data.frame("intercept"=Reg$coefficients[1,1],"slope"=Reg$coefficients[2,1])
    attr(res, "R2") <- Reg$r.squared
    res
  })
  
  R2s<-unlist(lapply(a, function(x) attributes(x)$R2))
  R2s<-data.frame(median=median(R2s), lwr=quantile(R2s, 0.025), upr=quantile(R2s, 0.975))
  a<-do.call(rbind,a)
  attr(a, "R2") <- R2s
  return(a)
}


#" Bootstrap a weighted OLS regression models from a calibration dataset
#" 
#" @param data The calibration dataset
#" @param replicates Number of bootstrap replicates
#" @param samples Number of samples per replicate
#" @param D47error The column in data containing the errors in D47


simulateLM_inverseweights <<- function(data, 
                                     replicates, 
                                     samples = NULL, 
                                     D47error="D47error"){
  a<-lapply(1:replicates, function(x){
    dataSub<-data[sample(seq_along(data[,1]), if(is.null(samples)){nrow(data)}else{samples}, replace = T),]
    dataSub$y_SE<-dataSub[,D47error]
    dataSub$x_SE<-abs(dataSub$TempError)
    Reg0<-lm(D47 ~ Temperature,  dataSub)
    wt <- 1 / lm(abs(Reg0$residuals) ~ Reg0$fitted.values)$fitted.values^2
    Reg<-summary(lm(D47 ~ Temperature,  dataSub, weights = wt))
    res<-cbind.data.frame("intercept"=Reg$coefficients[1,1],"slope"=Reg$coefficients[2,1])
    attr(res, "R2") <- Reg$r.squared
    res
  })
  
  R2s<-unlist(lapply(a, function(x) attributes(x)$R2))
  R2s<-data.frame(median=median(R2s), lwr=quantile(R2s, 0.025), upr=quantile(R2s, 0.975))
  a<-do.call(rbind,a)
  attr(a, "R2") <- R2s
  return(a)
}


#" Bootstrap Deming regression models from a calibration dataset
#" 
#" @param data The calibration dataset
#" @param replicates Number of bootstrap replicates
#" @param samples Number of samples per replicate
#" @param D47error The column in data containing the errors in D47


simulateDeming <<- function(data, 
                          replicates, 
                          samples = NULL, 
                          D47error="D47error", 
                          multicore=TRUE){
  
  if(multicore){
  
  ncores = parallel::detectCores()
  do.call(rbind,mclapply(1:replicates, mc.cores = ncores, function(x){
    dataSub<-data[sample(seq_along(data[,1]), if(is.null(samples)){nrow(data)}else{samples}, replace = T),]
    dataSub$y_SE<-abs(dataSub[,D47error])/sqrt(nrow(dataSub))
    dataSub$x_SE<-abs(dataSub$TempError)/sqrt(nrow(dataSub))
    Reg<-deming(D47 ~ Temperature, dataSub, xstd= 1/x_SE^2, ystd= 1/y_SE^2)
    cbind.data.frame("intercept"=Reg$coefficients[1],"slope"=Reg$coefficients[2])
  }))
  
  }else{
    do.call(rbind,lapply(1:replicates, mc.cores = ncores, function(x){
      dataSub<-data[sample(seq_along(data[,1]), if(is.null(samples)){nrow(data)}else{samples}, replace = T),]
      dataSub$y_SE<-abs(dataSub[,D47error])/sqrt(nrow(dataSub))
      dataSub$x_SE<-abs(dataSub$TempError)/sqrt(nrow(dataSub))
      Reg<-deming(D47 ~ Temperature, dataSub, xstd= 1/x_SE^2, ystd= 1/y_SE^2)
      cbind.data.frame("intercept"=Reg$coefficients[1],"slope"=Reg$coefficients[2])
    }))
  }
  
}


#" Bootstrap Bayesian regression models from a calibration dataset
#" 
#" @param data The calibration dataset
#" @param replicates Number of bootstrap replicates
#" @param samples Number of samples per replicate
#" @param generations number of MCMC generations
#" @param isMixed should we run only the mixed model?
#" @param priors Informative priors or not?
#" @param multicore run the analyses using multiple cores?
#" @param D47error The column in data containing the errors in D47


simulateBLM_measuredMaterial <<- function(data, 
                                        replicates, 
                                        samples = NULL, 
                                        generations=20000, 
                                        isMixed=FALSE,
                                        priors = "Informative", 
                                        multicore=TRUE){
  
  data_BR_Measured<-data
  
  single_rep<-function(i){
    
    dataSub <- data_BR_Measured
    
    if(length(unique(data_BR_Measured$Material)) == 1 ){
      dataSub<-data_BR_Measured[sample(seq_along(data_BR_Measured[,1]), if(is.null(samples)){nrow(data_BR_Measured)}else{samples}, replace = T),]

    }else{
      dataSub<- do.call(rbind,lapply(unique(data_BR_Measured$Material), function(x){
        material1<-data_BR_Measured[data_BR_Measured$Material ==x,]
        dataSub1<-material1[sample(seq_along(material1[,1]), round((if(is.null(samples)){nrow(data_BR_Measured)}else{samples}) /length(unique(data_BR_Measured$Material))), replace = T),]
        dataSub1
      } ))
    }
    
    Reg<-fitClumpedRegressions(calibrationData=dataSub,
                               hasMaterial = isMixed, n.iter = generations,
                               priors = priors)
    
    if(isMixed == FALSE){
      to_ret<-list(
        cbind.data.frame("intercept"=Reg$BLM1_fit$BUGSoutput$summary[1,1],"slope"=Reg$BLM1_fit$BUGSoutput$summary[2,1]),
        cbind.data.frame("intercept"=Reg$BLM1_fit_NoErrors$BUGSoutput$summary[1,1],"slope"=Reg$BLM1_fit_NoErrors$BUGSoutput$summary[2,1]))
      attr(to_ret, "R2s") <- attr(Reg, "R2s") 
      attr(to_ret, "DICs") <- attr(Reg, "DICs") 
      attr(to_ret, "Conv") <- list(BLM1_fit=Reg$BLM1_fit$BUGSoutput$summary, 
                                   BLM1_fit_NoErrors= Reg$BLM1_fit_NoErrors$BUGSoutput$summary)
      attr(to_ret, "PosteriorOne") <- list(BLM1_fit=Reg$BLM1_fit$BUGSoutput$sims.matrix, BLM1_fit_NoErrors=Reg$BLM1_fit_NoErrors$BUGSoutput$sims.matrix)
      
      to_ret
    }else{
      
      nmaterials<-length(unique(data_BR_Measured$Material))
      
      to_ret<-list(
        NA,
        NA,
        cbind.data.frame("intercept"=Reg$BLM3_fit$BUGSoutput$summary[c(1:nmaterials),1],"slope"=Reg$BLM3_fit$BUGSoutput$summary[c((nmaterials+1):c(nmaterials+nmaterials)),1], 
                         "material"=unique(dataSub$Material )))
      attr(to_ret, "R2s") <- attr(Reg, "R2s") 
      attr(to_ret, "DICs") <- attr(Reg, "DICs") 
      attr(to_ret, "Conv") <- list(#BLM1_fit=Reg$BLM1_fit$BUGSoutput$summary, 
                                   #BLM1_fit_NoErrors= Reg$BLM1_fit_NoErrors$BUGSoutput$summary,
                                   BLM3_fit=Reg$BLM3_fit$BUGSoutput$summary
      )
      attr(to_ret, "PosteriorOne") <- Reg$BLM3_fit$BUGSoutput$sims.matrix
      
      to_ret
    }
    
  }
  
  # Find out how many cores there are
  
  if(multicore){
  ncores = parallel::detectCores()
  tot = mclapply(1:replicates, mc.cores = ncores, single_rep)
  }else{
  tot = lapply(1:replicates, single_rep)
  }
  
  if(isMixed == FALSE){
    
    to_ret<-list("BLM_Measured_errors"=
           do.call(rbind,lapply(tot, function(x) x[[1]])),
         "BLM_Measured_no_errors"=do.call(rbind,lapply(tot, function(x) x[[2]]))
    )
    
    rs2<-do.call(rbind,lapply(tot, function(x) attr(x, "R2s") ))
    rs2<-aggregate(rs2[, 1:3], list(rs2$model), median)
    attr(to_ret, "R2s") <- rs2
    
    DICs<-lapply(tot, function(x) attr(x, "DICs") )
    
    DICs<-do.call(rbind,lapply(1:2 , function(x){ 
      a<-unlist(lapply(seq_along(DICs), function(y){ DICs[[y]][x] }))
      cbind.data.frame(median=median(a), lwr=quantile(a, 0.025), upr=quantile(a, 0.975), model=names(a[1]))
      }))

    attr(to_ret, "DICs") <- DICs
    
    Conv<-lapply(tot, function(x) attr(x, "Conv") )
    attr(to_ret, "Conv") <- Conv
    attr(to_ret, "PosteriorOne") <- attr(tot[[1]], "PosteriorOne")
    
    to_ret
  }else{
    
    targetlist<-lapply(tot, function(x) x[[3]])
    if(any(sapply(targetlist, is.null))){targetlist<-targetlist[-which(sapply(targetlist, is.null))]}
    BLMMFin<-do.call(rbind,targetlist)
    #BLMMFin<-BLMMFin[grep("[",row.names(BLMMFin), fixed = T),]
    
    to_ret<-list("BLM_Measured_errors"=
           do.call(rbind,lapply(tot, function(x) x[[1]])),
         "BLM_Measured_no_errors"=do.call(rbind,lapply(tot, function(x) x[[2]])),
         "BLMM_Measured_errors"= BLMMFin
    )
    rs2<-do.call(rbind,lapply(tot, function(x) attr(x, "R2s") ))
    rs21<-aggregate(rs2[, 1:3], list(rs2$class,rs2$model), median)
    rs22<-aggregate(rs2[, 1:3], list(rs2$class,rs2$model), FUN=quantile, probs=0.025)[,3]
    rs23<-aggregate(rs2[, 1:3], list(rs2$class,rs2$model), FUN=quantile, probs=0.975)[,3]
    rs21$lwr<-rs22
    rs21$upr<-rs23
    
    attr(to_ret, "R2s") <- rs21
    
    DICs<-lapply(tot, function(x) attr(x, "DICs") )
    
    DICs<-do.call(rbind,lapply(1:3 , function(x){ 
      a<-unlist(lapply(seq_along(DICs), function(y){ DICs[[y]][x] }))
      cbind.data.frame(median=median(a), lwr=quantile(a, 0.025), upr=quantile(a, 0.975), model=names(a[1]))
    }))
    attr(to_ret, "DICs") <- DICs
    
    Conv<-lapply(tot, function(x) attr(x, "Conv") )
    attr(to_ret, "Conv") <- Conv
    attr(to_ret, "PosteriorOne") <- attr(tot[[1]], "PosteriorOne")
    
    to_ret
  }
  
}


#" Estimate the R2 for Bayesian models
#" 
#" @param model The model (from R2jags)
#" @param calibrationData the calibration dataset used to fit the model
#" @param hasMaterial Is it the mixed model?

getR2Bayesian <<- function(model, 
                           calibrationData, 
                           hasMaterial=FALSE) {
  if(hasMaterial==F){
    mcmc <- model$BUGSoutput$sims.matrix
    Xmat = model.matrix( ~ Temperature, calibrationData)
    coefs = mcmc[, c("alpha", "beta")]
    fit = coefs %*% t(Xmat)
    resid = sweep(fit, 2, calibrationData$D47, "-")
    var_f = apply(fit, 1, var)
    var_e = apply(resid, 1, var)
    R2 = var_f / (var_f + ifelse(is.na(var_e), 0,var_e))
    cbind.data.frame(
      mean = mean(R2),
      lwr = quantile(R2, 0.025),
      upr = quantile(R2, 0.975)
    )
  }else{
    mcmc <- model$BUGSoutput$sims.matrix
    Xmat = model.matrix( ~ Temperature, calibrationData)
    coefs = as.data.frame(mcmc[,-ncol(mcmc) ])
    new <- rowSums(coefs[,grep("alpha", names(coefs)), drop=FALSE])
    new2 <- rowSums(coefs[,grep("beta", names(coefs)), drop=FALSE])
    coefs<-cbind(new,new2)
    fit = coefs %*% t(Xmat)
    resid = sweep(fit, 2, calibrationData$D47, "-")
    var_f = apply(fit, 1, var)
    var_e = apply(resid, 1, var)
    R2 = var_f / (var_f + var_e)
    cbind.data.frame(
      mean = mean(R2),
      lwr = quantile(R2, 0.025),
      upr = quantile(R2, 0.975)
    )
    
  }
}


#" This function is used to generate CI estimates at given intervals. It is currently
#" used for plotting in BayClump.
#" 
#" @param data A data.frame with two columns labeled slope and intercept. 
#" This should be the result of bootstrapping a given calibration set.
#" @param from the lower limit in x
#" @param to the upper limit in x
#" @param length.out the number of breaks


RegressionSingleCI<-function(data, from, to, length.out=100){
  
  sampleDataReplicates<- as.data.frame(data)
  
  Theta_mat <- sampleDataReplicates
  
  # Points where to evaluate the model
  x_eval <- seq(from, to, length.out = length.out)
  
  # Matrix with the predictions
  fun <- function(x, theta) as.numeric(theta["intercept"]) + (x * as.numeric(theta["slope"]))
  Pred_mat <- apply(Theta_mat, 1, function(theta) fun(x_eval, theta))
  
  
  # Pack the estimates for plotting
  uncertaintyModels <- cbind(
    x = x_eval, 
    as.data.frame(t(apply(Pred_mat, 1, function(y_est) c(
      mean_est = mean(y_est), 
      ci_lower_est = quantile(y_est, probs = 0.025, names = FALSE), 
      ci_upper_est = quantile(y_est, probs = 0.975, names = FALSE)
    ))))
  )
  
  return(list(uncertaintyModels))
  
}

#" Generate a dataset reflecting the priors used to run the analyses
#" 
#" @param prior Informative or not
#" @param n number of observations to simulate

generatePriorDistCalibration <- function(prior, n=1000){
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
  
  data <- cbind.data.frame(alpha=rnorm(n, params[1,2], params[1,3]), 
                   beta=rnorm(n, params[2,2], params[2,3]))
  attr(data, "priors") <-  prior
  attr(data, "params") <-  params
  data
}



