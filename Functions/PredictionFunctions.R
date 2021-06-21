CItoSE<-function(upper, lower){
  (upper-lower)/3.92
}
predictTcNonBayes<-function(data, slope, slpcnf, intercept, intcnf ){
  
  errors<-as.data.frame(data)
  ##Ignore uncertainty in slope and intercept (the usual)
  
  resnoParUn<-do.call(rbind,lapply(1:nrow(errors), function(x){
    resnoParUn<-clumpedr::revcal(D47=c(errors[x,1],errors[x,1]+(errors[x,2]*1.96), errors[x,1]-(errors[x,2]*1.96)),
                                 slope = slope, intercept = intercept,
                                 slpcnf =  slpcnf, intcnf =  intcnf,
                                 ignorecnf = T)
    
    data.frame(D47=errors[x,1], D47error=errors[x,2] ,'Tc'=resnoParUn[1], 'lwr'=resnoParUn[2], 'upr'=resnoParUn[3])
  }))
  
  ##Account for uncertainty in parameters
  
  resParUn<-do.call(rbind,lapply(1:nrow(errors), function(x){
    resParUn<-clumpedr::revcal(D47=c(errors[x,1],errors[x,1]+(errors[x,2]*1.96), errors[x,1]-(errors[x,2]*1.96)),
                               slope = as.numeric(slope), intercept = as.numeric(intercept),
                               slpcnf =  as.numeric(slpcnf), intcnf =  as.numeric(intcnf),
                               ignorecnf = F)
    
    data.frame(D47=errors[x,1], D47error=errors[x,2] ,'Tc'=as.numeric(resParUn[1,1]), 'lwr'=as.numeric(resParUn[2,4]), 'upr'=as.numeric(resParUn[3,3]))
  }))
  
  
  rbind(cbind.data.frame(type='No parameter Uncertainty',resnoParUn),
        cbind.data.frame(type='Parameter Uncertainty', resParUn))
  
}
predictTcBayes<-function(calibrationData, data, generations,hasMaterial=F, onlyWithinBayesian=F){
  
  errors<-as.data.frame(data)
  
  if(isFALSE(onlyWithinBayesian)){
  
  SingleRep<-fitClumpedRegressions(calibrationData=calibrationData,
                                   hasMaterial = hasMaterial, n.iter = generations)
  
  if(hasMaterial){
    
    a<-predictTcNonBayes(data=errors, slope=SingleRep$BLM1_fit$BUGSoutput$summary[2,1], 
                 slpcnf=CItoSE(SingleRep$BLM1_fit$BUGSoutput$summary[2,7], SingleRep$BLM1_fit$BUGSoutput$summary[2,3]), 
                 intercept=SingleRep$BLM1_fit$BUGSoutput$summary[1,1], 
                 intcnf=CItoSE(SingleRep$BLM1_fit$BUGSoutput$summary[1,7], SingleRep$BLM1_fit$BUGSoutput$summary[1,3]) )
    
    b<-predictTcNonBayes(data=errors, slope=SingleRep$BLM1_fit_NoErrors$BUGSoutput$summary[2,1], 
                 slpcnf=CItoSE(SingleRep$BLM1_fit_NoErrors$BUGSoutput$summary[2,7], SingleRep$BLM1_fit_NoErrors$BUGSoutput$summary[2,3]), 
                 intercept=SingleRep$BLM1_fit_NoErrors$BUGSoutput$summary[1,1], 
                 intcnf=CItoSE(SingleRep$BLM1_fit_NoErrors$BUGSoutput$summary[1,7], SingleRep$BLM1_fit_NoErrors$BUGSoutput$summary[1,3]) )
    
    c<-predictTcNonBayes(data=errors, slope=SingleRep$BLM3_fit$BUGSoutput$summary[2,1], 
                 slpcnf=CItoSE(SingleRep$BLM3_fit$BUGSoutput$summary[2,7], SingleRep$BLM3_fit$BUGSoutput$summary[2,3]), 
                 intercept=SingleRep$BLM3_fit$BUGSoutput$summary[1,1], 
                 intcnf=CItoSE(SingleRep$BLM3_fit$BUGSoutput$summary[1,7], SingleRep$BLM3_fit$BUGSoutput$summary[1,3]) )
  }else{
    
    a<-predictTcNonBayes(data=errors, slope=SingleRep$BLM1_fit$BUGSoutput$summary[2,1], 
                 slpcnf=CItoSE(SingleRep$BLM1_fit$BUGSoutput$summary[2,7], SingleRep$BLM1_fit$BUGSoutput$summary[2,3]), 
                 intercept=SingleRep$BLM1_fit$BUGSoutput$summary[1,1], 
                 intcnf=CItoSE(SingleRep$BLM1_fit$BUGSoutput$summary[1,7], SingleRep$BLM1_fit$BUGSoutput$summary[1,3]) )
    
    b<-predictTcNonBayes(data=errors, slope=SingleRep$BLM1_fit_NoErrors$BUGSoutput$summary[2,1], 
                 slpcnf=CItoSE(SingleRep$BLM1_fit_NoErrors$BUGSoutput$summary[2,7], SingleRep$BLM1_fit_NoErrors$BUGSoutput$summary[2,3]), 
                 intercept=SingleRep$BLM1_fit_NoErrors$BUGSoutput$summary[1,1], 
                 intcnf=CItoSE(SingleRep$BLM1_fit_NoErrors$BUGSoutput$summary[1,7], SingleRep$BLM1_fit_NoErrors$BUGSoutput$summary[1,3]) )
    
  }
  
  
  }
  
  # Bayesian predictions inside the Bayesian framework (full propagation)
  
  predictionsWithinBayesian<-fitClumpedRegressionsPredictions(calibrationData=calibrationData, 
                                                              useInits=F, 
                                                              hasMaterial = hasMaterial,
                                                              D47Prederror=errors[,2],
                                                              D47Pred=errors[,1],
                                                              n.iter= generations)
  
  
  predsComplete<-if(hasMaterial){
    
    fullProp<- if(nrow(errors)==1 ){
      
    rbind(
      cbind.data.frame(model='BLM1_fit', errors, t(matrix(predictionsWithinBayesian$BLM1_fit$BUGSoutput$summary[c(1:nrow(errors)),c(5,3,7)]))),
      cbind.data.frame(model='BLM1_fit_NoErrors',errors,t(matrix(predictionsWithinBayesian$BLM1_fit_NoErrors$BUGSoutput$summary[c(1:nrow(errors)),c(5,3,7)]))),
      cbind.data.frame(model='BLM3_fit',errors,t(matrix(predictionsWithinBayesian$BLM3_fit$BUGSoutput$summary[c(1:nrow(errors)),c(5,3,7)])))
    )
    }else{
      rbind(
        cbind.data.frame(model='BLM1_fit', errors, predictionsWithinBayesian$BLM1_fit$BUGSoutput$summary[c(1:nrow(errors)),c(5,3,7)]),
        cbind.data.frame(model='BLM1_fit_NoErrors',errors,predictionsWithinBayesian$BLM1_fit_NoErrors$BUGSoutput$summary[c(1:nrow(errors)),c(5,3,7)]),
        cbind.data.frame(model='BLM3_fit',errors,predictionsWithinBayesian$BLM3_fit$BUGSoutput$summary[c(1:nrow(errors)),c(5,3,7)])
      ) 
      }
      
    fullProp$type<-"Parameter Uncertainty"
    fullProp$BayesianPredictions<-'Yes'
    
    if(isFALSE(onlyWithinBayesian)){
    a$model<-'BLM1_fit'
    b$model<-'BLM1_fit_NoErrors'
    c$model<-'BLM3_fit'
    colnames(fullProp)<-c('model', 'D47', 'D47error', 'Tc', 'lwr', 'upr', 'type', 'BayesianPredictions')
    rbindlist(list(fullProp,a,b,c),fill=TRUE)
    }else{
      colnames(fullProp)<-c('model', 'D47', 'D47error', 'Tc', 'lwr', 'upr', 'type', 'BayesianPredictions')
      fullProp
      }
    
  }else{
    
      fullProp<- 
        if(nrow(errors)==1 ){
        rbind(
        cbind(model='BLM1_fit', errors, t(matrix(predictionsWithinBayesian$BLM1_fit$BUGSoutput$summary[c(1:nrow(errors)),c(5,3,7)]))),
        cbind.data.frame(model='BLM1_fit_NoErrors',errors,t(matrix(predictionsWithinBayesian$BLM1_fit_NoErrors$BUGSoutput$summary[c(1:nrow(errors)),c(5,3,7)])))  )
        }else{
          rbind(
            cbind(model='BLM1_fit', errors, predictionsWithinBayesian$BLM1_fit$BUGSoutput$summary[c(1:nrow(errors)),c(5,3,7)]),
            cbind.data.frame(model='BLM1_fit_NoErrors',errors,predictionsWithinBayesian$BLM1_fit_NoErrors$BUGSoutput$summary[c(1:nrow(errors)),c(5,3,7)]))
        }
      
      
      fullProp$type<-"Parameter Uncertainty"
      fullProp$BayesianPredictions<-'Yes'
      if(isFALSE(onlyWithinBayesian)){
      a$model<-'BLM1_fit'
      b$model<-'BLM1_fit_NoErrors'
      
      colnames(fullProp)<-c('model', 'D47', 'D47error', 'Tc', 'lwr', 'upr', 'type', 'BayesianPredictions')
      rbindlist(list(fullProp,a,b),fill=TRUE)
    }else{
      colnames(fullProp)<-c('model', 'D47', 'D47error', 'Tc', 'lwr', 'upr', 'type', 'BayesianPredictions')
      fullProp
    }
    

    
    
  }
  
  
  predsComplete
}
