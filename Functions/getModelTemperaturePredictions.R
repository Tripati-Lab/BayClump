##Predictions table
getModelTemperaturePredictions<-function(CompleteModelFit,predictionData ){
  
  if(length(CompleteModelFit) ==7){
    df_preds<-cbind.data.frame(predictionData, predictYork(CompleteModelFit$Y, predictionData$D47, predictionData$D47_SD, predictionData$n_measurements),
                               predictLm(CompleteModelFit$M0, predictionData$D47, predictionData$D47_SD, predictionData$n_measurements),
                               predictLmCovariate(attr(CompleteModelFit, "data"),CompleteModelFit$M1, predictionData$D47, predictionData$D47_SD, predictionData$n_measurements, predictionData$Material),
                               predictLmCovariate(attr(CompleteModelFit, "data"),CompleteModelFit$M2, predictionData$D47, predictionData$D47_SD, predictionData$n_measurements, predictionData$Material),
                               predictLmJags(CompleteModelFit$BLM1_fit, predictionData$D47, predictionData$D47_SD, predictionData$n_measurements),
                               predictANCOVA1Jags(CompleteModelFit$BLM2_fit, predictionData$D47, predictionData$D47_SD, predictionData$n_measurements, predictionData$Material),
                               predictANCOVA2Jags(CompleteModelFit$BLM3_fit, predictionData$D47_SD, predictionData$n_measurements, predictionData$D47, predictionData$Material))
  }else{
    df_preds<-cbind.data.frame(predictionData,predictYork(CompleteModelFit$Y, predictionData$D47, predictionData$D47_SD,  predictionData$n_measurements),
                               predictLm(CompleteModelFit$M0, predictionData$D47, predictionData$D47_SD, predictionData$n_measurements),
                               predictLmJags(CompleteModelFit$BLM1_fit, predictionData$D47, predictionData$D47_SD, predictionData$n_measurements))
  }
  names(df_preds)<-make.unique(names(df_preds), sep="_")
  return(df_preds)
  
}
