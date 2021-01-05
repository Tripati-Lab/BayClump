##Predictions table
getModelD47Predictions<-function(CompleteModelFit,predictionData){
  if(length(CompleteModelFit) == 7){
    df_preds<-cbind.data.frame(predictionData, CalibratePredictYork(CompleteModelFit$Y, predictionData),
                               CalibratePredictLm(CompleteModelFit$M0, predictionData$T2),
                               CalibratePredictLmCovariate(attr(CompleteModelFit, "data"),CompleteModelFit$M1, predictionData$T2, predictionData$Material),
                               CalibratePredictLmCovariate(attr(CompleteModelFit, "data"), CompleteModelFit$M2, predictionData$T2, predictionData$Material),
                               CalibratePredictLmJags(CompleteModelFit$BLM1_fit, predictionData$T2),
                               CalibratePredictANCOVA1Jags(CompleteModelFit$BLM2_fit, predictionData$T2, predictionData$Material),
                               CalibratePredictANCOVA2Jags(CompleteModelFit$BLM3_fit, predictionData$T2, predictionData$Material))
  }else{
    df_preds<-cbind.data.frame(predictionData, CalibratePredictYork(CompleteModelFit$Y, predictionData),
                               CalibratePredictLm(CompleteModelFit$M0, predictionData$T2),
                               CalibratePredictLmJags(CompleteModelFit$BLM1_fit, predictionData$T2))
  }
  names(df_preds)<-make.unique(names(df_preds), sep="_")
  return(df_preds)
  
}
