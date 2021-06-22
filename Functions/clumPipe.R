lw<-function(x) quantile(x, 0.025)
up<-function(x) quantile(x, 0.975)

##Closest functions
dist <- function(df2, D47error, TError, bestCase=T){
  dt <- data.table((df2$Dataset_D47Error-D47error)^2+(df2$Dataset_TError-TError)^2)
  return(if(bestCase){which.min(dt$V1)}else{which.max(dt$V1)})
}

distD47 <- function(df2, TargetD47error, Target_D47, bestCase=T){
  dt <- data.table((df2$Target_D47Error-TargetD47error)^2+(df2$Target_D47-Target_D47)^2)
  return(if(bestCase){which.min(dt$V1)}else{which.max(dt$V1)})
}


clumpipe<-function(calData, PipCriteria, targetD47, error_targetD47, nrep=1000, BayesianOnly=F){
  
  ##Find overall error scenario for the dataset
  
  targetScenario<-unlist(PipCriteria[dist(df2=PipCriteria, 
                                          D47error=mean(calData$D47error), 
                                          TError=mean(calData$TempError), bestCase = T),1])
  
  
  subPipCriteria<-PipCriteria[PipCriteria$Error_scenario ==targetScenario, ]
  
  ##Find model for each target D47
  
  compilationModels<-cbind.data.frame(targetD47, error_targetD47, do.call(rbind,lapply(1:length(targetD47), function(x){
    #Don't change the bestCase argument
    selRow<-distD47(subPipCriteria, TargetD47error=error_targetD47[x], Target_D47=targetD47[x], bestCase = T)
    subPipCriteria[selRow,]
  })))
  
  
  
  if(BayesianOnly==T){
    
    
    singleRep<-function(i) {predictTcBayes(calibrationData=calData, 
                                           data=cbind(targetD47,error_targetD47),
                                           generations=10000, 
                                           hasMaterial=F, onlyWithinBayesian=T)
    }
    
    totalRep<-pbmclapply(1:nrep,singleRep, mc.cores = 4)
    uncertaintyPredictionWithinBayesianModels<-do.call(rbind,totalRep)
    uncertaintyPredictionWithinBayesianModels<-uncertaintyPredictionWithinBayesianModels[uncertaintyPredictionWithinBayesianModels$model== 'BLM1_fit',]
    
    uncertaintyPredictionWithinBayesianModels<- lapply(1:(nrow(uncertaintyPredictionWithinBayesianModels)/nrep), function(x){
     targetRows<- uncertaintyPredictionWithinBayesianModels[x,c('D47','D47error')]
     
     subData<-uncertaintyPredictionWithinBayesianModels[which(uncertaintyPredictionWithinBayesianModels$D47 == targetRows$D47 &
             uncertaintyPredictionWithinBayesianModels$D47error ==targetRows$D47error ) ,]
     
     cbind.data.frame(model=uncertaintyPredictionWithinBayesianModels$model[1],
                      D47=targetRows$D47,
                      D47error=targetRows$D47error,
                      Tc=median(subData$Tc),
                      lwr=lw(subData$Tc),
                      upr=up(subData$Tc) )

    })
    uncertaintyPredictionWithinBayesianModels<-do.call(rbind,uncertaintyPredictionWithinBayesianModels)
    

    colnames(compilationModels)[c(1,2)]<-c('D47','D47error')
    completeDf<-cbind.data.frame(compilationModels,uncertaintyPredictionWithinBayesianModels)
    
    return(completeDf)

  }else{
    
    
    ##Bayesian predictions
    
    BP<-compilationModels[compilationModels$PredictionType == 'Bayesian predictions',]
    NBP<-compilationModels[compilationModels$PredictionType != 'Bayesian predictions',]
    
    
    
    
    ##Perform non-Bayesian predictions
    if(nrow(NBP)>0){
      prediction_NBP<-do.call(rbind,lapply(1:nrow(NBP), function(x){
        SM<-NBP[x,'ModelSelected']
        
        ##Run the model
        if(SM == 'LM'){
          valsReps<- simulateLM_measured(data=calData, replicates=nrep, samples=50)
          
        }
        if(SM=='York'){
          valsReps<- simulateYork_measured(data=calData, replicates=nrep, samples=50)
          
        }
        if(SM=='Deming'){
          valsReps<-  simulateDeming(data=calData, replicates=nrep, samples=50)
          
        }
        
        if(SM=='invweLM'){
          valsReps<- simulateLM_inverseweights(data=calData, replicates=nrep, samples=50)
          
        }
        
        if(length(grep('BLM',SM))>0){
          ##Just do the one with errors for now
          valsReps<-simulateBLM_measuredMaterial(data=calData, replicates=nrep, samples=50, generations=10000, isMixed=F)[[1]]
        }
        
        ##Get the predictions
        
        predictTcNonBayes(data=cbind(NBP$targetD47[x],NBP$error_targetD47[x]), 
                          slope=median(valsReps$slope), 
                          slpcnf=CItoSE(quantile(valsReps$slope, 0.975), quantile(valsReps$slope, 0.025)), 
                          intercept=median(valsReps$intercept), 
                          intcnf=CItoSE(quantile(valsReps$intercept, 0.975), quantile(valsReps$intercept, 0.025)))[1,]
        
      }))
      if(nrow(BP)==0){
        return(cbind.data.frame(compilationModels,prediction_NBP))
      }
    }
    
    
    if(nrow(BP)>0){
      
      ##Perform Bayesian predictions
      
      singleRep<-function(i) {predictTcBayes(calibrationData=calData, 
                                             data=cbind(BP$targetD47,BP$error_targetD47),
                                             generations=10000, 
                                             hasMaterial=F, onlyWithinBayesian=T)[1,]
      }
      
      totalRep<-pbmclapply(1:nrep,singleRep, mc.cores = 4)
      uncertaintyPredictionWithinBayesianModels<-do.call(rbind,totalRep)
      uncertaintyPredictionWithinBayesianModels<-summaryBy(Tc ~ model+D47+D47error, data = uncertaintyPredictionWithinBayesianModels, 
                                                           FUN = list(median,lw,up))
      colnames(uncertaintyPredictionWithinBayesianModels)<-c("model","D47","D47error","Tc", "lwr","upr")
      
      
      resPreds<-rbindlist(list(prediction_NBP,uncertaintyPredictionWithinBayesianModels))
      
      colnames(compilationModels)[c(1,2)]<-c('D47','D47error')
      completeDf<-merge(compilationModels, resPreds, by=c('D47','D47error'))
      return(completeDf)
      
    }
  }
  
}
