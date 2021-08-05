
simulateYork_measured<<-function(data, replicates, samples=NULL, D47error="D47error"){
  do.call(rbind,pblapply(1:replicates, function(x){
    dataSub<-data[sample(seq_along(data[,1]), if(is.null(samples)){nrow(data)}else{samples}, replace = T),]
    dataSub$y_SE<-dataSub[,D47error]
    dataSub$x_SE<-dataSub$TempError
    Reg<-york(cbind.data.frame(dataSub$Temperature, dataSub$x_SE, dataSub$D47, dataSub$y_SE))
    cbind.data.frame('intercept'=Reg$a[1],'slope'=Reg$b[1])
  }))
}


simulateLM_measured<<-function(data, replicates, samples=NULL, D47error="D47error"){
  do.call(rbind,pblapply(1:replicates, function(x){
    dataSub<-data[sample(seq_along(data[,1]), if(is.null(samples)){nrow(data)}else{samples}, replace = T),]
    dataSub$y_SE<-dataSub[,D47error]
    dataSub$x_SE<-dataSub$TempError
    Reg<-summary(lm(D47 ~ Temperature,  dataSub))
    cbind.data.frame('intercept'=Reg$coefficients[1,1],'slope'=Reg$coefficients[2,1])
  }))
}


simulateLM_inverseweights<<-function(data, replicates, samples=NULL, D47error="D47error"){
  do.call(rbind,pblapply(1:replicates, function(x){
    dataSub<-data[sample(seq_along(data[,1]), if(is.null(samples)){nrow(data)}else{nrow(data)*samples}, replace = T),]
    dataSub$y_SE<-dataSub[,D47error]
    dataSub$x_SE<-dataSub$TempError
    Reg<-summary(lm(D47 ~ Temperature,  dataSub, weights = y_SE))
    cbind.data.frame('intercept'=Reg$coefficients[1,1],'slope'=Reg$coefficients[2,1])
  }))
}


simulateDeming<<-function(data, replicates, samples=NULL, D47error="D47error"){
  do.call(rbind,pblapply(1:replicates, function(x){
    dataSub<-data[sample(seq_along(data[,1]), if(is.null(samples)){nrow(data)}else{nrow(data)*samples}, replace = T),]
    dataSub$y_SE<-dataSub[,D47error]
    dataSub$x_SE<-dataSub$TempError
    Reg<-deming(D47 ~ Temperature, dataSub, xstd= 1/dataSub$x_SE, ystd= 1/dataSub$y_SE)
    cbind.data.frame('intercept'=Reg$coefficients[1],'slope'=Reg$coefficients[2])
  }))
}


simulateBLM_measuredMaterial<<-function(data, replicates, samples=NULL, generations=20000, isMixed=F){
  
  data_BR_Measured<-data
  
  single_rep<-function(i){
    
    
    if(isMixed == F){
      dataSub<-data_BR_Measured[sample(seq_along(data_BR_Measured[,1]), if(is.null(samples)){nrow(data)}else{samples}, replace = T),]
      
    }else{
      dataSub<- do.call(rbind,lapply(unique(data_BR_Measured$Material), function(x){
        material1<-data_BR_Measured[data_BR_Measured$Material ==x,]
        dataSub1<-material1[sample(seq_along(material1[,1]), round(if(is.null(samples)){nrow(data)}else{nrow(data)*samples}/length(unique(data_BR_Measured$Material))), replace = T),]
        dataSub1
      } ))
      
    }
    
    Reg<-fitClumpedRegressions(calibrationData=dataSub,
                               hasMaterial = isMixed, n.iter = generations)
    
    
    
    if(isMixed==F ){
      to_ret<-list(
        cbind.data.frame('intercept'=Reg$BLM1_fit$BUGSoutput$summary[1,1],'slope'=Reg$BLM1_fit$BUGSoutput$summary[2,1]),
        cbind.data.frame('intercept'=Reg$BLM1_fit_NoErrors$BUGSoutput$summary[1,1],'slope'=Reg$BLM1_fit_NoErrors$BUGSoutput$summary[2,1]))
      attr(to_ret, 'R2s') <- attr(Reg, "R2s") 
      attr(to_ret, 'DICs') <- attr(Reg, "DICs") 
      to_ret
    }else{
      
      nmaterials<-length(unique(data_BR_Measured$Material))
      
      to_ret<-list(
        cbind.data.frame('intercept'=Reg$BLM1_fit$BUGSoutput$summary[1,1],'slope'=Reg$BLM1_fit$BUGSoutput$summary[2,1]),
        cbind.data.frame('intercept'=Reg$BLM1_fit_NoErrors$BUGSoutput$summary[1,1],'slope'=Reg$BLM1_fit_NoErrors$BUGSoutput$summary[2,1]),
        cbind.data.frame('intercept'=Reg$BLM3_fit$BUGSoutput$summary[c(1:nmaterials),1],'slope'=Reg$BLM3_fit$BUGSoutput$summary[c((nmaterials+1):c(nmaterials+nmaterials)),1], 
                         'material'=unique(dataSub$Material )))
      attr(to_ret, 'R2s') <- attr(Reg, "R2s") 
      attr(to_ret, 'DICs') <- attr(Reg, "DICs") 
      to_ret
    }
    
  }
  
  # Find out how many cores there are
  ncores = parallel::detectCores()
  
  # Use all available cores
  tot = pbmclapply(1:replicates, mc.cores = ncores, single_rep)
  #tot = pblapply(1:replicates, single_rep)
  
  if(isMixed == F){
    
    to_ret<-list('BLM_Measured_errors'=
           do.call(rbind,lapply(tot, function(x) x[[1]])),
         'BLM_Measured_no_errors'=do.call(rbind,lapply(tot, function(x) x[[2]]))
    )
    
    rs2<-do.call(rbind,lapply(tot, function(x) attr(x, "R2s") ))
    rs2<-aggregate(rs2[, 1:3], list(rs2$model), mean)
    attr(to_ret, 'R2s') <- rs2
    
    DICs<-lapply(tot, function(x) attr(x, "DICs") )
    
    DICs<-do.call(rbind,lapply(1:2 , function(x){ 
      a<-unlist(lapply(seq_along(DICs), function(y){ DICs[[y]][x] }))
      cbind.data.frame(median=median(a), lwr=quantile(a, 0.025), upr=quantile(a, 0.975), model=names(a[1]))
      }))

    attr(to_ret, 'DICs') <- DICs
    to_ret
  }else{
    
    targetlist<-lapply(tot, function(x) x[[3]])
    if(any(sapply(targetlist, is.null))){targetlist<-targetlist[-which(sapply(targetlist, is.null))]}
    BLMMFin<-do.call(rbind,targetlist)
    #BLMMFin<-BLMMFin[grep("[",row.names(BLMMFin), fixed = T),]
    
    to_ret<-list('BLM_Measured_errors'=
           do.call(rbind,lapply(tot, function(x) x[[1]])),
         'BLM_Measured_no_errors'=do.call(rbind,lapply(tot, function(x) x[[2]])),
         'BLMM_Measured_errors'= BLMMFin
    )
    rs2<-do.call(rbind,lapply(tot, function(x) attr(x, "R2s") ))
    rs2<-aggregate(rs2[, 1:3], list(rs2$model), mean)
    attr(to_ret, 'R2s') <- rs2
    DICs<-lapply(tot, function(x) attr(x, "DICs") )
    
    DICs<-do.call(rbind,lapply(1:3 , function(x){ 
      a<-unlist(lapply(seq_along(DICs), function(y){ DICs[[y]][x] }))
      cbind.data.frame(median=median(a), lwr=quantile(a, 0.025), upr=quantile(a, 0.975), model=names(a[1]))
    }))
    attr(to_ret, 'DICs') <- DICs
    to_ret
  }
  
  
  
}