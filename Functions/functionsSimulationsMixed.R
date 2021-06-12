
simulateYork_measured<<-function(data, replicates, samples=NULL, D47error="D47error"){
  do.call(rbind,pblapply(1:replicates, function(x){
    dataSub<-data[sample(seq_along(data[,1]), if(is.null(samples)){nrow(data)}else{nrow(data)*samples}, replace = T),]
    dataSub$y_SE<-dataSub[,D47error]
    dataSub$x_SE<-dataSub$Temp_Error
    Reg<-york(cbind.data.frame(dataSub$T2, dataSub$x_SE, dataSub$D47, dataSub$y_SE))
    cbind.data.frame('intercept'=Reg$a[1],'slope'=Reg$b[1])
  }))
}


simulateLM_measured<<-function(data, replicates, samples=NULL, D47error="D47error"){
  do.call(rbind,pblapply(1:replicates, function(x){
    dataSub<-data[sample(seq_along(data[,1]), if(is.null(samples)){nrow(data)}else{nrow(data)*samples}, replace = T),]
    dataSub$y_SE<-dataSub[,D47error]
    dataSub$x_SE<-dataSub$Temp_Error
    Reg<-summary(lm(D47~ T2,  dataSub, weights = y_SE))
    cbind.data.frame('intercept'=Reg$coefficients[1,1],'slope'=Reg$coefficients[2,1])
  }))
}


simulateLM_inverseweights<<-function(data, replicates, samples=NULL, D47error="D47error"){
  do.call(rbind,pblapply(1:replicates, function(x){
    dataSub<-data[sample(seq_along(data[,1]), if(is.null(samples)){nrow(data)}else{nrow(data)*samples}, replace = T),]
    Reg<-summary(lm(D47~ T2,  dataSub))
    cbind.data.frame('intercept'=Reg$coefficients[1,1],'slope'=Reg$coefficients[2,1])
  }))
}


simulateDeming<<-function(data, replicates, samples=NULL, D47error="D47error"){
  do.call(rbind,pblapply(1:replicates, function(x){
    dataSub<-data[sample(seq_along(data[,1]), if(is.null(samples)){nrow(data)}else{nrow(data)*samples}, replace = T),]
    dataSub$y_SE<-dataSub[,D47error]
    dataSub$x_SE<-dataSub$Temp_Error
    Reg<-deming(D47 ~ T2, dataSub, xstd= dataSub$x_SE, ystd= dataSub$y_SE)
    cbind.data.frame('intercept'=Reg$coefficients[1],'slope'=Reg$coefficients[2])
  }))
}


simulateBLM_measuredMaterial<<-function(data, replicates, samples=NULL, generations=20000, isMixed=F){
  
  data_BR_Measured<-data
  
  single_rep<-function(i){
    
    
    if(isMixed == F){
      dataSub<-data_BR_Measured[sample(seq_along(data_BR_Measured[,1]), if(is.null(samples)){nrow(data)}else{nrow(data)*samples}, replace = T),]
      
    }else{
      material1<-data_BR_Measured[data_BR_Measured$Material ==1,]
      dataSub1<-material1[sample(seq_along(material1[,1]), round(if(is.null(samples)){nrow(data)}else{nrow(data)*samples}/2), replace = T),]
      
      material2<-data_BR_Measured[data_BR_Measured$Material ==2,]
      dataSub2<-material2[sample(seq_along(material2[,1]), round(if(is.null(samples)){nrow(data)}else{nrow(data)*samples}/2), replace = T),]
      
      dataSub<-rbind.data.frame(dataSub1,dataSub2)
    }
    
    Reg<-fitClumpedRegressions(calibrationData=dataSub,
                               hasMaterial = isMixed, n.iter = generations)
    
    
    if(isMixed==F ){
      list(
        cbind.data.frame('intercept'=Reg$BLM1_fit$BUGSoutput$summary[1,1],'slope'=Reg$BLM1_fit$BUGSoutput$summary[2,1]),
        cbind.data.frame('intercept'=Reg$BLM1_fit_NoErrors$BUGSoutput$summary[1,1],'slope'=Reg$BLM1_fit_NoErrors$BUGSoutput$summary[2,1]))
      
    }else{
      
      list(
        cbind.data.frame('intercept'=Reg$BLM1_fit$BUGSoutput$summary[1,1],'slope'=Reg$BLM1_fit$BUGSoutput$summary[2,1]),
        cbind.data.frame('intercept'=Reg$BLM1_fit_NoErrors$BUGSoutput$summary[1,1],'slope'=Reg$BLM1_fit_NoErrors$BUGSoutput$summary[2,1]),
        cbind.data.frame('intercept'=Reg$BLM3_fit$BUGSoutput$summary[c(1,2),1],'slope'=Reg$BLM3_fit$BUGSoutput$summary[c(3,4),1], 
                         'material'=unique(dataSub$Material )))
      
    }
    
  }
  
  tot = lapply(1:replicates, single_rep)
  
  
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