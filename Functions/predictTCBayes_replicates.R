predictTcBayes_replicates<-function(calData, 
                                    targetD47, 
                                    error_targetD47, 
                                    material, 
                                    nrep=1000, 
                                    hasMaterial=F, 
                                    generations=20000){
  
    singleRep<-function(i) {predictTcBayes(calibrationData=calData, 
                                           data=cbind.data.frame(targetD47,error_targetD47, material),
                                           generations=generations, 
                                           hasMaterial=hasMaterial)
    }
    
    ncores = parallel::detectCores()
    
    totalRep<-pbmclapply(1:nrep,singleRep, mc.cores = ncores)
    
    uncertaintyPredictionWithinBayesianModels<-do.call(rbind,totalRep)
    
    uncertaintyPredictionWithinBayesianModelsComplete<- lapply(unique(uncertaintyPredictionWithinBayesianModels$model), function(y){
    uncertaintyPredictionWithinBayesianModels<-uncertaintyPredictionWithinBayesianModels[uncertaintyPredictionWithinBayesianModels$model== y,]
    
    uncertaintyPredictionWithinBayesianModels<- do.call(rbind,lapply(1:(nrow(uncertaintyPredictionWithinBayesianModels)/nrep), function(x){
      targetRows<- uncertaintyPredictionWithinBayesianModels[x,c('D47','D47error')]
      subData<-uncertaintyPredictionWithinBayesianModels[which(uncertaintyPredictionWithinBayesianModels$D47 == targetRows$D47 &
                                                                 uncertaintyPredictionWithinBayesianModels$D47error ==targetRows$D47error ) ,]
      cbind.data.frame(model=uncertaintyPredictionWithinBayesianModels$model[1],
                       D47=targetRows$D47,
                       D47error=targetRows$D47error,
                       Tc=median(subData$Tc),
                       lwr=median(subData$lwr),
                       upr=median(subData$upr) )
    }))
    
    colnames(uncertaintyPredictionWithinBayesianModels)<-c("model","D47","D47error","Tc", "lwr","upr")
    
    
    '%nin%'<-Negate('%in%')
    for(x in 1:nrow(uncertaintyPredictionWithinBayesianModels)){
      ds<-uncertaintyPredictionWithinBayesianModels[x,c("Tc","lwr",'upr')]
      up<-ds[which.max(ds)]
      lw<-ds[which.min(ds)]
      Tc<-ds[which(ds %nin% c(up,lw))]
      uncertaintyPredictionWithinBayesianModels[x,c("Tc","lwr",'upr')]<-c(Tc, lw,up)
    }
    
    uncertaintyPredictionWithinBayesianModels
    })
    
    return(uncertaintyPredictionWithinBayesianModelsComplete)
    
 
  
}
