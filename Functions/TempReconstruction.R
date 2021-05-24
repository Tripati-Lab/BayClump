TempReconstructionCoefficients<-function(
                             dataSlopeIntercept, 
                             targety=c(0.6,0.7), 
                             errory=c(0.005, 0.01), 
                             nmeasurements=c(2,3)){


modelUncertaintyPredictions<-do.call(rbind,lapply(seq_along(errory), function(d){
    Theta_mat <- dataSlopeIntercept
    
    errors<-cbind.data.frame(targety=targety, errory=c(errory), 
                             nmeasurements=nmeasurements)
    
    reps_complete<- do.call(rbind,lapply(1:100, function(i){
      # Points where to evaluate the model
      y_eval<-sapply(targety, function(w) rnorm(1, mean=w, sd= errory[d]* sqrt(errors$nmeasurements)))
      
      
      # Matrix with the predictions
      fun <- function(x, theta) (x-as.numeric(theta["intercept"])) / ( as.numeric(theta["slope"]))
      Pred_mat <- apply(Theta_mat, 1, function(theta) fun(y_eval, theta))
      
      
      # Pack the estimates for plotting
      Estims_plot <- cbind(
        replicate=i,
        y = targety, 
        errory=errory[d],
        as.data.frame(t(apply(Pred_mat, 1, function(x_est) c(
          median_est = median(x_est), 
          ci_lower_est = quantile(x_est, probs = 0.025, names = FALSE), 
          ci_upper_est = quantile(x_est, probs = 0.975, names = FALSE)
        ))))      )
      
    }))
  }))

sum<-cbind.data.frame(ddply(modelUncertaintyPredictions,.(y),function(x) median(x$median_est)),
                      lowerCI=ddply(modelUncertaintyPredictions,.(y),function(x) quantile(x$median_est, 0.025))[,2],
                      upperCI=ddply(modelUncertaintyPredictions,.(y),function(x) quantile(x$median_est, 0.975))[,2])

colnames(sum)[c(1:2)]<-c('targetD47', 'MedianPrediction')

return(list('summary'=sum, 'replicates'=modelUncertaintyPredictions))

}

TempReconstructionBayesian<-function(dataCalibration, 
                             targety=c(0.6,0.7), 
                             errory=c(0.005, 0.01), 
                             materialPred=c(1,2), 
                             nmeasurements=c(2,3),
                             hasMaterial=T){
  
  errors<-cbind.data.frame(targety=targety, errory=c(errory), 
                           materialPred=materialPred, 
                           nmeasurements=nmeasurements)
  
  perfBayesianModelComplete<-do.call(rbind,lapply(1:20, function(i){
    tryCatch({
      print(i)
      dataCal<-dataCalibration
      compsingleRep<- lapply(1:nrow(errors), function(y){
        
        calibrationDataSub<- if(hasMaterial == T ){
          ddply(dataCal,.(Material),function(x) x[sample(nrow(x),25, replace = T),])
        }else{
          dataCal[sample(1:nrow(dataCal), 50, replace = T), ]
        }
        
        k<-fitClumpedRegressionsPredictions(calibrationData=calibrationDataSub, 
                                            useInits=F, 
                                            hasMaterial = hasMaterial,
                                            D47Prederror=errors$errory[y],
                                            D47Pred=errors$targety[y],
                                            materialPred=errors$materialPred[y],
                                            n.iter= 20000)
        
        pu<-list('BLM1_fit'=as.data.frame(t(as.data.frame(k$BLM1_fit$BUGSoutput$summary[1,c(1,2,3,5,7)]))),
                 'BLM1_fit_NoErrors'=as.data.frame(t(as.data.frame(k$BLM1_fit_NoErrors$BUGSoutput$summary[1,c(1,2,3,5,7)]))),
                 'BLM3_fit'=as.data.frame(t(as.data.frame(k$BLM3_fit$BUGSoutput$summary[1,c(1,2,3,5,7)]))))
        pu2<-cbind(rbindlist(pu, idcol = T), errors[y,])
        pu2
      })
      
      compsingleRep2 <- rbindlist(compsingleRep, idcol = T)
      
      
      return(compsingleRep2)
    }, error=function(e){})
  }))
  
  colnames(perfBayesianModelComplete)[c(1,2,8,9,6,5,7)]<-c('replicate','model', 'y', 'errory', 'median_est','ci_lower_est','ci_upper_est')
  
  
  sum<-cbind.data.frame(ddply(perfBayesianModelComplete,.(y),function(x) median(x$median_est)),
                   lowerCI=ddply(perfBayesianModelComplete,.(y),function(x) quantile(x$median_est, 0.025))[,2],
                   upperCI=ddply(perfBayesianModelComplete,.(y),function(x) quantile(x$median_est, 0.975))[,2])
  
  colnames(sum)[c(1:2)]<-c('targetD47', 'MedianPrediction')
  
  return(list('summary'=sum, 'replicates'=perfBayesianModelComplete))
}


