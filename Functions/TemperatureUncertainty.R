Temperature95Uncertainty<-function(CompleteModelFit, dataset){
  
  tp<-getModelTemperaturePredictions(CompleteModelFit=RF_improper, predictionData=dataset)
  tp<-tp[,grep('er95',colnames(tp))]
  tp<-as.data.frame(sapply(tp, as.numeric))
  do.call(rbind.data.frame,lapply(seq(1,ncol(tp), by=2), function(x){
    
    cbind.data.frame(Model=strsplit(colnames(tp[x]), '_')[[1]][1],
                     mean_95endpoints_distance=mean(tp[,x+1]-tp[,x], na.rm=T),
                     max_95endpoints_distance=max(tp[,x+1]-tp[,x], na.rm=T),
                     min_95endpoints_distance=min(tp[,x+1]-tp[,x], na.rm=T))
  }))
  
}