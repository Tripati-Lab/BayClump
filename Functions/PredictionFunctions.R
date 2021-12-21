predictTcBayes<-function(calibrationData, data, generations,hasMaterial=T){
  
  errors<-data
  if(ncol(errors) < 3 ){ errors<-cbind(errors,Material=1)}
  
  predictionsWithinBayesian<-fitClumpedRegressionsPredictions(calibrationData=calibrationData, 
                                                              useInits=T, 
                                                              hasMaterial = hasMaterial,
                                                              D47Pred=errors[,1],
                                                              D47Prederror=errors[,2],
                                                              materialsPred=errors[,3],
                                                              n.iter= generations)
  
  predsComplete<-if(hasMaterial){
    
    fullProp<-rbind(
      cbind.data.frame(model='BLM1_fit', errors, median=predictionsWithinBayesian$BLM1_fit$BUGSoutput$summary[c(1:nrow(errors)),c(5)],
                       lwr=predictionsWithinBayesian$BLM1_fit$BUGSoutput$summary[c(1:nrow(errors)),c(7)],
                       upr=predictionsWithinBayesian$BLM1_fit$BUGSoutput$summary[c(1:nrow(errors)),c(3)]
      ),
      cbind.data.frame(model='BLM1_fit_NoErrors',errors,median=predictionsWithinBayesian$BLM1_fit_NoErrors$BUGSoutput$summary[c(1:nrow(errors)),c(5)],
                       lwr=predictionsWithinBayesian$BLM1_fit_NoErrors$BUGSoutput$summary[c(1:nrow(errors)),c(7)],
                       upr=predictionsWithinBayesian$BLM1_fit_NoErrors$BUGSoutput$summary[c(1:nrow(errors)),c(3)]),
      cbind.data.frame(model='BLM3_fit',errors,median=predictionsWithinBayesian$BLM3_fit$BUGSoutput$summary[c(1:nrow(errors)),c(5)],
                       lwr=predictionsWithinBayesian$BLM3_fit$BUGSoutput$summary[c(1:nrow(errors)),c(7)],
                       upr=predictionsWithinBayesian$BLM3_fit$BUGSoutput$summary[c(1:nrow(errors)),c(3)])
    )
    
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
