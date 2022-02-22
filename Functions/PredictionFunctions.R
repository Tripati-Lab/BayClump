predictTcBayes <- function(calibrationData, 
                           data, 
                           generations, 
                           hasMaterial=T, 
                           bootDataset=T, 
                           onlyMedian=F, 
                           replicates = 1000, 
                           multicore= TRUE, 
                           priors='informative',
                           errorsD47=T){
  
  single_rep <<- function(i){
  
  errors<-data
  if(ncol(errors) < 3 ){ errors<-cbind(errors,Material=1)}
  
  if(bootDataset){calibrationData <- calibrationData[sample(1:nrow(calibrationData),nrow(calibrationData), replace = T),] }
  
  if(errorsD47){
    
    predictionsWithinBayesian<-fitClumpedRegressionsPredictions_D47errors(calibrationData=calibrationData, 
                                                                    useInits=T, 
                                                                    hasMaterial = hasMaterial,
                                                                    D47Pred=errors[,1],
                                                                    D47Prederror=errors[,2],
                                                                    materialsPred=errors[,3],
                                                                    n.iter= generations,
                                                                    priors=priors)
  }else{
  
  predictionsWithinBayesian<-fitClumpedRegressionsPredictions_D47(calibrationData=calibrationData, 
                                                              useInits=T, 
                                                              hasMaterial = hasMaterial,
                                                              D47Pred=errors[,1],
                                                              D47Prederror=errors[,2],
                                                              materialsPred=errors[,3],
                                                              n.iter= generations,
                                                              priors=priors)
  }
  predsComplete<-if(hasMaterial){
    
    fullProp<-
      cbind.data.frame(model='BLM3_fit',errors,median=predictionsWithinBayesian$BLM3_fit$BUGSoutput$summary[c(1:nrow(errors)+1),c(5)],
                       lwr=predictionsWithinBayesian$BLM3_fit$BUGSoutput$summary[c(1:nrow(errors)+1),c(3)],
                       upr=predictionsWithinBayesian$BLM3_fit$BUGSoutput$summary[c(1:nrow(errors)+1),c(7)])
    
  

      colnames(fullProp)<-c('model', 'D47', 'D47error','Material' ,'Tc', 'lwr', 'upr')
      fullProp

  }else{
    
    fullProp<-  rbind(
      cbind.data.frame(model='BLM1_fit', errors, median=predictionsWithinBayesian$BLM1_fit$BUGSoutput$summary[c(1:nrow(errors)+1),c(5)],
                       lwr=predictionsWithinBayesian$BLM1_fit$BUGSoutput$summary[c(1:nrow(errors)+1),c(3)],
                       upr=predictionsWithinBayesian$BLM1_fit$BUGSoutput$summary[c(1:nrow(errors)+1),c(7)]
      ),
      cbind.data.frame(model='BLM1_fit_NoErrors',errors,median=predictionsWithinBayesian$BLM1_fit_NoErrors$BUGSoutput$summary[c(1:nrow(errors)+1),c(5)],
                       lwr=predictionsWithinBayesian$BLM1_fit_NoErrors$BUGSoutput$summary[c(1:nrow(errors)+1),c(3)],
                       upr=predictionsWithinBayesian$BLM1_fit_NoErrors$BUGSoutput$summary[c(1:nrow(errors)+1),c(7)])  ) 
    
      colnames(fullProp)<-c('model', 'D47', 'D47error',"Material", 'Tc', 'lwr', 'upr')
      fullProp
  }

  row.names(predsComplete) <-NULL
  predsComplete
  }
  
  if(bootDataset){
    
    # Find out how many cores there are
    ncores = parallel::detectCores()
    
    # Use all available cores
    tot =if( multicore ){ 
      pbmclapply(1:replicates, mc.cores = ncores, single_rep) } else {
              lapply(1:replicates, single_rep)
        }
    tot <- do.call(rbind,tot)
    
    if(onlyMedian){
     dat <- ddply(tot,~model+D47+D47error+Material,summarise, median=median(Tc), lwr= quantile(Tc, 0.025), upr= quantile(Tc, 0.975))
     names(dat)[5] <- "Tc"
     dat$Tc <- sqrt(10^6/dat$Tc)-273.15
     dat$lwr <- sqrt(10^6/dat$lwr)-273.15
     dat$upr <- sqrt(10^6/dat$upr)-273.15
     dat$se <- (dat$lwr - dat$upr) / 3.92
     dat$sd <- dat$se * sqrt(replicates)
     dat[,c(1:5,9)]
    }else{
      dat <- ddply(tot,~model+D47+D47error+Material,summarise, median=median(Tc), lwr= median(lwr), upr= median(upr) )
      names(dat)[5] <- "Tc"
      dat$Tc <- sqrt(10^6/dat$Tc)-273.15
      dat$lwr <- sqrt(10^6/dat$lwr)-273.15
      dat$upr <- sqrt(10^6/dat$upr)-273.15
      dat$se <- (dat$lwr - dat$upr) / 3.92
      dat$sd <- dat$se * sqrt(replicates)
      dat[,c(1:5,9)]
    }
    
  }else{
    single_rep()
  }
  
}
