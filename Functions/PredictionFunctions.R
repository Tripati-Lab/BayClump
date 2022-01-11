predictTcBayes <- function(calibrationData, data, generations, hasMaterial=T, bootDataset=T, onlyMedian=T, replicates = 1000, multicore= TRUE, priors='informative'){
  
  single_rep <<- function(i){
  
  errors<-data
  if(ncol(errors) < 3 ){ errors<-cbind(errors,Material=1)}
  
  if(bootDataset){calibrationData <- calibrationData[sample(1:nrow(calibrationData),nrow(calibrationData), replace = T),] }
  
  predictionsWithinBayesian<-fitClumpedRegressionsPredictions(calibrationData=calibrationData, 
                                                              useInits=T, 
                                                              hasMaterial = hasMaterial,
                                                              D47Pred=errors[,1],
                                                              D47Prederror=errors[,2],
                                                              materialsPred=errors[,3],
                                                              n.iter= generations,
                                                              priors=priors
                                                              )
  
  predsComplete<-if(hasMaterial){
    
    fullProp<-
      cbind.data.frame(model='BLM3_fit',errors,median=predictionsWithinBayesian$BLM3_fit$BUGSoutput$summary[c(1:nrow(errors)),c(5)],
                       lwr=predictionsWithinBayesian$BLM3_fit$BUGSoutput$summary[c(1:nrow(errors)),c(7)],
                       upr=predictionsWithinBayesian$BLM3_fit$BUGSoutput$summary[c(1:nrow(errors)),c(3)])
    
    
    fullProp$type<-"Parameter Uncertainty"
    fullProp$BayesianPredictions<-'Yes'
    

      colnames(fullProp)<-c('model', 'D47', 'D47error','Material' ,'Tc', 'lwr', 'upr', 'type', 'BayesianPredictions')
      fullProp

  }else{
    
    fullProp<-  rbind(
      cbind.data.frame(model='BLM1_fit', errors, median=predictionsWithinBayesian$BLM1_fit$BUGSoutput$summary[c(1:nrow(errors)),c(5)],
                       lwr=predictionsWithinBayesian$BLM1_fit$BUGSoutput$summary[c(1:nrow(errors)),c(7)],
                       upr=predictionsWithinBayesian$BLM1_fit$BUGSoutput$summary[c(1:nrow(errors)),c(3)]
      ),
      cbind.data.frame(model='BLM1_fit_NoErrors',errors,median=predictionsWithinBayesian$BLM1_fit_NoErrors$BUGSoutput$summary[c(1:nrow(errors)),c(5)],
                       lwr=predictionsWithinBayesian$BLM1_fit_NoErrors$BUGSoutput$summary[c(1:nrow(errors)),c(7)],
                       upr=predictionsWithinBayesian$BLM1_fit_NoErrors$BUGSoutput$summary[c(1:nrow(errors)),c(3)])  ) 
    
    fullProp$type<-"Parameter Uncertainty"
    fullProp$BayesianPredictions<-'Yes'
  
      colnames(fullProp)<-c('model', 'D47', 'D47error',"Material", 'Tc', 'lwr', 'upr', 'type', 'BayesianPredictions')
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
              lapply(1:replicates, single_rep)}
    tot <- do.call(rbind,tot)
    
    if(onlyMedian){
     dat <- ddply(tot,~model+D47+D47error+Material,summarise, median=mean(Tc), se=(sd(Tc)/sqrt(length(Tc))))
     names(dat)[5] <- "Tc"
     dat
    }else{
      tot %>% 
        group_by(model, D47, D47error, Material) %>%
        summarise(across(-c(type,BayesianPredictions), median, na.rm = TRUE)) %>% 
        as.data.frame()
      names(dat)[5] <- "Tc"
      dat
    }
    
  }else{
    single_rep()[,c(1:7)]
  }
  
}
