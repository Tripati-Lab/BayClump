# Dataset must include T2 and D47.

plotObservedPredictedDistributions<-function(CompleteModelFit, dataset, plotOnly=T){
  ##Observed and predicted Temperature
  checkModelTemps<-getModelTemperaturePredictions(CompleteModelFit, dataset)
  id.vars<-grep('T2|mean', colnames(checkModelTemps))
  names(checkModelTemps) <- paste0("value.", names(checkModelTemps))
  checkModelTemps<-reshape(checkModelTemps[,id.vars],                    # data
          direction = "long",    # long or wide
          varying = 1:length(id.vars),         # the columns that should be stacked
          timevar = "condition"  # name of "time" variable, basically groups
  )
  colnames(checkModelTemps)[1]<-c('Model')
  
  p<-ggplot(checkModelTemps, aes(x = value, y = Model)) +
    geom_density_ridges(aes(fill = Model))
  
  
  ##Observed and predicted D47
  checkModelD47<-getModelD47Predictions(CompleteModelFit, dataset)
  
  id.varsModels<-grep('mean', colnames(checkModelD47))
  id.vars<-grep('\\D47\\b', colnames(checkModelD47))
  id.vars<-c(id.vars,id.varsModels)
  
  names(checkModelD47) <- paste0("value.", names(checkModelD47))
  
  checkModelD47<-reshape(checkModelD47[,id.vars],                    # data
                           direction = "long",    # long or wide
                           varying = 1:length(id.vars),         # the columns that should be stacked
                           timevar = "condition"  # name of "time" variable, basically groups
  )
  colnames(checkModelD47)[1]<-c('Model')
  
  q<-ggplot(checkModelD47, aes(x = value, y = Model)) +
    geom_density_ridges(aes(fill = Model))
  
  if(plotOnly==T){
    ggarrange(p, q)
  }else{
    list('plot'=ggarrange(p, q),'checkModelTemps'=checkModelTemps,'checkModelD47'=checkModelD47)
  }
}
