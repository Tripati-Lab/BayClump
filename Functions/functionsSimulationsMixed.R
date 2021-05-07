
simulateYork_measured<<-function(data, replicates=100, samples=10){
  do.call(rbind,pblapply(1:replicates, function(x){
    dataSub<-data[sample(seq_along(data[,1]), samples, replace = T),]
    Reg<-york(cbind.data.frame(dataSub$T2, dataSub$Temp_Error, dataSub$D47, dataSub$D47_SD))
    cbind.data.frame('intercept'=Reg$a[1],'slope'=Reg$b[1])
  }))
}


simulateLM_measured<<-function(data, replicates=100, samples=30){
  do.call(rbind,pblapply(1:replicates, function(x){
    dataSub<-data[sample(seq_along(data[,1]), samples, replace = T),]
    dataSub$y_SE<-dataSub$D47_SD/sqrt(2)
    Reg<-summary(lm(D47~ T2,  dataSub, weights = 1/y_SE^2))
    cbind.data.frame('intercept'=Reg$coefficients[1,1],'slope'=Reg$coefficients[2,1])
  }))
}


simulateLM_inverseweights<<-function(data, replicates=100, samples=30){
  do.call(rbind,pblapply(1:replicates, function(x){
    dataSub<-data[sample(seq_along(data[,1]), samples, replace = T),]
    Reg<-summary(lm(D47~ T2,  dataSub))
    cbind.data.frame('intercept'=Reg$coefficients[1,1],'slope'=Reg$coefficients[2,1])
  }))
}


simulateDeming<<-function(data, replicates=100, samples=30){
  do.call(rbind,pblapply(1:replicates, function(x){
    dataSub<-data[sample(seq_along(data[,1]), samples, replace = T),]
    dataSub$y_SE<-abs(dataSub$D47_SD/sqrt(20))
    dataSub$x_SE<-abs(dataSub$Temp_Error/sqrt(20))
    Reg<-deming(D47 ~ T2, dataSub, xstd= 1/dataSub$x_SE, ystd= dataSub$y_SE)
    cbind.data.frame('intercept'=Reg$coefficients[1],'slope'=Reg$coefficients[2])
  }))
}


simulateBLM_measuredMaterial<<-function(data, replicates=100, samples=30, generations=20000, isMixed=F){
  
  data_BR_Measured<-data
  
  single_rep<-function(i){
    
    
    if(isMixed == F){
      dataSub<-data_BR_Measured[sample(seq_along(data_BR_Measured[,1]), samples, replace = T),]
      
    }else{
      material1<-data_BR_Measured[data_BR_Measured$Material ==1,]
      dataSub1<-material1[sample(seq_along(material1[,1]), round(samples/2), replace = T),]
      
      material2<-data_BR_Measured[data_BR_Measured$Material ==2,]
      dataSub2<-material2[sample(seq_along(material2[,1]), round(samples/2), replace = T),]
      
      dataSub<-rbind.data.frame(dataSub1,dataSub2)
    }
    
    Reg<-fitClumpedRegressions(calibrationData=dataSub,
                               hasMaterial = T, n.iter = generations)
    
    
    if(length(unique(dataSub$Material))==1 ){
      list(
        cbind.data.frame('intercept'=Reg$BLM1_fit$BUGSoutput$summary[1,1],'slope'=Reg$BLM1_fit$BUGSoutput$summary[2,1]),
        cbind.data.frame('intercept'=Reg$BLM1_fit_NoErrors$BUGSoutput$summary[1,1],'slope'=Reg$BLM1_fit_NoErrors$BUGSoutput$summary[2,1]),
        cbind.data.frame('intercept'=Reg$BLM3_fit$BUGSoutput$summary[c(1),1],'slope'=Reg$BLM3_fit$BUGSoutput$summary[c(2),1]))
      
    }else{
      
      list(
        cbind.data.frame('intercept'=Reg$BLM1_fit$BUGSoutput$summary[1,1],'slope'=Reg$BLM1_fit$BUGSoutput$summary[2,1]),
        cbind.data.frame('intercept'=Reg$BLM1_fit_NoErrors$BUGSoutput$summary[1,1],'slope'=Reg$BLM1_fit_NoErrors$BUGSoutput$summary[2,1]),
        cbind.data.frame('intercept'=Reg$BLM3_fit$BUGSoutput$summary[c(1,2),1],'slope'=Reg$BLM3_fit$BUGSoutput$summary[c(3,4),1], 
                         'material'=c(1,2)))
      
    }
    
  }
  
  tot = mclapply(1:replicates, mc.cores = 4, single_rep)
  
  
  if(isMixed == F){
    
    list('BLM_Measured_errors'=
           do.call(rbind,lapply(tot, function(x) x[[1]])),
         'BLM_Measured_no_errors'=do.call(rbind,lapply(tot, function(x) x[[2]]))
    )
    
  }else{
    
    targetlist<-lapply(tot, function(x) x[[3]])
    if(any(sapply(targetlist, is.null))){targetlist<-targetlist[-which(sapply(targetlist, is.null))]}
    BLMMFin<-do.call(rbind,targetlist[unlist(lapply(targetlist, function(x) nrow(x)))==2 ])
    #BLMMFin<-BLMMFin[grep("[",row.names(BLMMFin), fixed = T),]
    
    list('BLM_Measured_errors'=
           do.call(rbind,lapply(tot, function(x) x[[1]])),
         'BLM_Measured_no_errors'=do.call(rbind,lapply(tot, function(x) x[[2]])),
         'BLMM_Measured_errors'= BLMMFin
    )
    
  }
  
  
  
}



#Analize all at once


simulateAll<-function(data, error, replicates, samples, generations, isMixed=F){
  York_measured_reps<-simulateYork_measured(data=data, replicates=replicates, samples=samples)
  LM_measured_reps<-simulateLM_measured(data=data, replicates=replicates, samples=samples)
  Deming<-simulateDeming(data=data, replicates=replicates, samples=samples)
  invweLM<-simulateLM_inverseweights(data=data, replicates=replicates, samples=samples)
  
  ##Bayesian
  BLM_measured<-simulateBLM_measuredMaterial(data=data, replicates=replicates, samples=samples, generations=generations, isMixed=isMixed)
  
  
  sal<-list('York'=York_measured_reps,
            'LM'=LM_measured_reps,
            'BLM_measured'=BLM_measured,
            'Deming'=Deming,
            'invweLM'=invweLM)
    
  full_replicates<- if(isTRUE(isMixed)){
  rbind.data.frame( 
    cbind.data.frame(rbindlist(sal[-c(3)], idcol=TRUE),material=NA),
    cbind.data.frame(rbindlist(sal[[3]][1:2], idcol=TRUE), material=NA),
    cbind.data.frame(.id= 'BLMM_Measured_errors', sal[[3]][[3]]))
}else{
  rbind.data.frame( 
    rbindlist(sal[-c(3)], idcol=TRUE),
    rbindlist(sal[[3]], idcol=TRUE))
}
  
  row.names(full_replicates)<-NULL
      
  return(full_replicates)
  
}