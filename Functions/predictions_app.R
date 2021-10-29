
#' @param targetD47 should have two columns (sample ID and D47 value)
#' @param sampleID Sample ID column in the targetD47 df

predictions_app<-function(calData, targetD47, sampleID, materials, nrep=1000, BayesianOnly=F, hasMaterial=F,generations=20000){
 
  bds <- aggregate(as.formula(paste0("D47 ~", sampleID)), targetD47, function(x) c(mean = mean(x), sd = sd(x)))
  bds <- bds[complete.cases(bds),]
  bds$Material <- sapply(bds$Sample.Name, function(x) targetD47[targetD47$Sample.Name == x,'Material'][1])
  nbds <- targetD47[,c("D47", sampleID)]
  
  if(BayesianOnly==T){
    
    singleRep<-function(i) {predictTcBayes(calibrationData=calData, 
                                           data=cbind.data.frame(bds$D47[,1],bds$D47[,2], bds$Material),
                                           generations=generations, 
                                           hasMaterial=T, onlyWithinBayesian=T)
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
                         lwr=max(subData$lwr),
                         upr=min(subData$upr) )
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
                                             generations=generations, 
                                             hasMaterial=F, onlyWithinBayesian=T)
      }
      
      totalRep<-pbmclapply(1:nrep,singleRep, mc.cores = 4)
      uncertaintyPredictionWithinBayesianModels<-do.call(rbind,totalRep)
      uncertaintyPredictionWithinBayesianModels<-uncertaintyPredictionWithinBayesianModels[uncertaintyPredictionWithinBayesianModels$model== 'BLM1_fit',]
      uncertaintyPredictionWithinBayesianModels<- do.call(rbind,lapply(1:(nrow(uncertaintyPredictionWithinBayesianModels)/nrep), function(x){
        targetRows<- uncertaintyPredictionWithinBayesianModels[x,c('D47','D47error')]
        subData<-uncertaintyPredictionWithinBayesianModels[which(uncertaintyPredictionWithinBayesianModels$D47 == targetRows$D47 &
                                                                   uncertaintyPredictionWithinBayesianModels$D47error ==targetRows$D47error ) ,]
        cbind.data.frame(model=uncertaintyPredictionWithinBayesianModels$model[1],
                         D47=targetRows$D47,
                         D47error=targetRows$D47error,
                         Tc=median(subData$Tc),
                         lwr=max(subData$lwr),
                         upr=min(subData$upr) )
      }))
      
      colnames(uncertaintyPredictionWithinBayesianModels)<-c("model","D47","D47error","Tc", "lwr","upr")
      resPreds<-rbindlist(list(prediction_NBP,uncertaintyPredictionWithinBayesianModels))
      
      colnames(compilationModels)[c(1,2)]<-c('D47','D47error')
      completeDf<-cbind.data.frame(compilationModels, resPreds)
      return(completeDf)
      
    }
  }
  
}
