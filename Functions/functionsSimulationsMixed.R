#' @param  n number of simulated observations
#' @param  sdx standard deviation on the simulated true x in Centigrades
#' @param  meanX Mean value for the simulated true X in 10^6/T2
#' @param  sdobs standard deviation on the measurement of X in 10^6/T2
#' @param  alpha True intercept
#' @param  beta True slope
#' @param  p_sdy Process error in Y
#' @param  obs_sdy standard deviation on the measurement of Y


sim_slr <- function(nobs, sdobsx, meanX, alpha, beta, sdy, sdobsy) {
  
  eqsdX<-data.frame(sdC=c(0.25,0.5,1,2,3,5,10), sd10K=c(.019,.038,.077,0.155,.234,.394,.808))
  sdobsx<-eqsdX[which(eqsdX[,1] >=  sdobsx ),2][1]

  truex <- rnorm(nobs,meanX,0.8759)    
  errx <- rnorm(nobs, 0, sdobsx)
  obsx <- truex + errx
  
  erry <- rnorm(nobs, 0, sdobsy)
  truey <- rnorm(nobs,alpha + beta*truex,sdy)
  obsy <- truey + erry

  res<-data.frame(x_TRUE = truex, x_measured= obsx, x_SD=errx, 
             y_TRUE=truey, y_SD=erry, y_measured=obsy, material=1)
  
    return(res)
}


#' Same as above but can handle multiple slopes and intercepts

sim_slr_mult <- function(nobs, sdobsx, meanX=12.02585352, alpha, beta, sdy, sdobsy) {
  
  nobs<-if(length(nobs)==1){c(nobs,nobs)}else{nobs}
  
  res<-do.call(rbind.data.frame,lapply(seq_along(alpha), function(x){
    eqsdX<-data.frame(sdC=c(0.25,0.5,1,2,3,5,10), sd10K=c(.019,.038,.077,0.155,.234,.394,.808))
    sdobsx<-eqsdX[which(eqsdX[,1] >=  sdobsx ),2][1]
    
    truex <- rnorm(nobs[x],meanX,0.8759)    
    errx <- rnorm(nobs[x], 0, sdobsx)
    obsx <- truex + errx
    
    erry <- rnorm(nobs[x], 0, sdobsy)
    truey <- rnorm(nobs[x],alpha[x] + beta[x]*truex,sdy)
    obsy <- truey + erry
    
    res<-data.frame(x_TRUE = truex, x_measured= obsx, x_SD=errx, 
                    y_TRUE=truey, y_SD=erry, y_measured=obsy, material=x)
  }))
  
  return(res)
}


sim_slr_givenX <- function(x, x_un, alpha, beta, sdy) {
  
  n=length(x)
  errorx <- x_un
  obsx <- x
  
  # simulate response data
  errory <-  sdy
  obsy <- alpha + beta*obsx + errory
  
  data.frame(x_measured= obsx, x_SD=errorx, 
             y_measured=obsy, y_SD=errory)
}


getBST_CI<-function(x, samples){
  sort(x)[c(round(samples*0.025), round(samples*0.975))]
}

simulateYork_true<<-function(data,error, replicates=100, samples=10){
  do.call(rbind,pblapply(1:replicates, function(x){
    dataSub<-data[sample(seq_along(data[,1]), samples, replace = T),]
    Reg<-york(cbind.data.frame(dataSub$x_TRUE, error, dataSub$y_measured, dataSub$y_SD))
    cbind.data.frame('intercept'=Reg$a[1],'slope'=Reg$b[1])
  }))
}



simulateYork_measured<<-function(data, replicates=100, samples=10){
  do.call(rbind,pblapply(1:replicates, function(x){
    dataSub<-data[sample(seq_along(data[,1]), samples, replace = T),]
    Reg<-york(cbind.data.frame(dataSub$x_measured, dataSub$x_SD, dataSub$y_measured, dataSub$y_SD))
    cbind.data.frame('intercept'=Reg$a[1],'slope'=Reg$b[1])
  }))
}

simulateLM_true<<-function(data, replicates=100, samples=10){
  do.call(rbind,pblapply(1:replicates, function(x){
    dataSub<-data[sample(seq_along(data[,1]), samples, replace = T),]
    Reg<-summary(lm(y_measured~ x_TRUE,  dataSub))
    cbind.data.frame('intercept'=Reg$coefficients[1,1],'slope'=Reg$coefficients[2,1])
  }))
}



simulateLM_measured<<-function(data, replicates=100, samples=30){
  do.call(rbind,pblapply(1:replicates, function(x){
    dataSub<-data[sample(seq_along(data[,1]), samples, replace = T),]
    dataSub$y_SE<-dataSub$y_SD/sqrt(2)
    Reg<-summary(lm(y_measured~ x_measured,  dataSub, weights = 1/y_SE^2))
    cbind.data.frame('intercept'=Reg$coefficients[1,1],'slope'=Reg$coefficients[2,1])
  }))
}


simulateLM_inverseweights<<-function(data, replicates=100, samples=30){
  do.call(rbind,pblapply(1:replicates, function(x){
    dataSub<-data[sample(seq_along(data[,1]), samples, replace = T),]
    Reg<-summary(lm(y_measured~ x_measured,  dataSub))
    cbind.data.frame('intercept'=Reg$coefficients[1,1],'slope'=Reg$coefficients[2,1])
  }))
}


simulateDeming<<-function(data, replicates=100, samples=30){
  do.call(rbind,pblapply(1:replicates, function(x){
    dataSub<-data[sample(seq_along(data[,1]), samples, replace = T),]
    Reg<-deming(y_measured ~ x_measured, dataSub)
    cbind.data.frame('intercept'=Reg$coefficients[1],'slope'=Reg$coefficients[2])
  }))
}


simulateBLM_measured<<-function(data, replicates=100, samples=30, generations=20000, dataGiven=F){
  
  data_BR_Measured<-if(dataGiven == F){data[,-c(1,4)] }else{
    data
  }
  
  colnames(data_BR_Measured)<-c('T2', 'Temp_Error','D47_SD','D47')
  
  single_rep<-function(i){
    dataSub<-data_BR_Measured[sample(seq_along(data_BR_Measured[,1]), samples, replace = T),]
    Reg<-fitClumpedRegressions(calibrationData=dataSub,
                               hasMaterial = F, n.iter = generations)
    
    list(
      cbind.data.frame('intercept'=Reg$BLM1_fit$BUGSoutput$summary[1,1],'slope'=Reg$BLM1_fit$BUGSoutput$summary[2,1]),
      cbind.data.frame('intercept'=Reg$BLM1_fit_NoErrors$BUGSoutput$summary[1,1],'slope'=Reg$BLM1_fit_NoErrors$BUGSoutput$summary[2,1])
    )
  }
  
  tot = pbmclapply(1:replicates, mc.cores = 4, single_rep)
  
  
  list('BLM_Measured_errors'=
         do.call(rbind,lapply(tot, function(x) x[[1]])),
       'BLM_Measured_no_errors'=do.call(rbind,lapply(tot, function(x) x[[2]]))
  )
  
}



simulateBLM_measuredMaterial<<-function(data, replicates=100, samples=30, generations=20000, dataGiven=F, isMixed=F){
  
  data_BR_Measured<-if(dataGiven == F){data[,-c(1,4)] }else{
    data
  }
  
  colnames(data_BR_Measured)<-c('T2', 'Temp_Error','D47_SD','D47', 'Material')
  
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
  
  tot = pbmclapply(1:replicates, mc.cores = 4, single_rep)
  
  
  if(isMixed == F){
    
    list('BLM_Measured_errors'=
           do.call(rbind,lapply(tot, function(x) x[[1]])),
         'BLM_Measured_no_errors'=do.call(rbind,lapply(tot, function(x) x[[2]])),
         'BLMM_Measured_errors'= do.call(rbind,lapply(tot, function(x) x[[3]]))
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


simulateBLM_true<<-function(data, replicates=100, samples=30, generations=20000){
  
  data_BR_Measured<-data[,-c(2,6)]
  colnames(data_BR_Measured)<-c('T2', 'Temp_Error','D47', 'D47_SD')
  
  single_rep<-function(i){
    dataSub<-data_BR_Measured[sample(seq_along(data_BR_Measured[,1]), samples, replace = T),]
    Reg<-fitClumpedRegressions(calibrationData=dataSub,
                               hasMaterial = F, n.iter = generations)
    list(
      cbind.data.frame('intercept'=Reg$BLM1_fit$BUGSoutput$summary[1,1],'slope'=Reg$BLM1_fit$BUGSoutput$summary[2,1]),
      cbind.data.frame('intercept'=Reg$BLM1_fit_NoErrors$BUGSoutput$summary[1,1],'slope'=Reg$BLM1_fit_NoErrors$BUGSoutput$summary[2,1])
    )
  }
  
  tot = pbmclapply(1:replicates, mc.cores = 4, single_rep)
  
  
  list('BLM_True_errors'=
         do.call(rbind,lapply(tot, function(x) x[[1]])),
       'BLM_True_no_errors'=do.call(rbind,lapply(tot, function(x) x[[2]]))
  )
  
}




#Analize all at once


simulateAll<-function(data, error, replicates, samples, generations, name='Simulations', isMixed=F){
  York_measured_reps<-simulateYork_measured(data=data, replicates=replicates, samples=samples)
  LM_measured_reps<-simulateLM_measured(data=data, replicates=replicates, samples=samples)
  Deming<-simulateDeming(data=data, replicates=replicates, samples=samples)
  invweLM<-simulateLM_inverseweights(data=data, replicates=replicates, samples=samples)
  
  ##Bayesian
  BLM_measured<-simulateBLM_measuredMaterial(data=data, replicates=replicates, samples=samples, generations=generations, isMixed=isMixed)

  
  
  
  return(list('York'=York_measured_reps,
              'LM'=LM_measured_reps,
              'BLM_measured'=BLM_measured,
              'Deming'=Deming,
              'invweLM'=invweLM
              ))
  
}

simulateAll_dataGiven<-function(data, error, replicates, samples, generations, name='Simulations'){
  York_measured_reps<-simulateYork_measured(data=data, replicates=replicates, samples=samples)
  LM_measured_reps<-simulateLM_measured(data=data, replicates=replicates, samples=samples)
  
  ##Bayesian
  BLM_measured<-simulateBLM_measured(data=data, replicates=replicates, samples=samples, generations=generations, dataGiven = T)

  
  ints<-c(York_measured_reps[,1], LM_measured_reps[,1],BLM_measured[[1]][,1],BLM_measured[[2]][,1])
  slopes<-c( York_measured_reps[,2], LM_measured_reps[,2],BLM_measured[[1]][,2],BLM_measured[[2]][,2])
  
  ran_in<-c(min(ints), max(ints))
  ran_sl<-c(min(slopes), max(slopes))
  
  
  pdf(paste0(name,'_replicates=',replicates,'_samples=',samples, '_generations=',generations ,'Standarized.pdf' ), 7,7)
  par(mfrow=c(4,2))
  
  quantiles1<-getBST_CI(York_measured_reps[,1],replicates);quantiles2<-getBST_CI(York_measured_reps[,2],replicates)
  hist(York_measured_reps[,1], col="lightblue", main = 'Intercept', xlim=ran_in);abline(v = 0.268, col="red", lwd=3, lty=2);abline(v = quantiles1, col="blue", lwd=3)
  mtext(paste0('True=',0.268, ", 95% CI= ", paste(round(quantiles1,4), collapse = '-')), 3, line=0)
  hist(York_measured_reps[,2], col="lightblue", main = 'Slope', xlim=ran_sl);abline(v = 0.0369, col="red", lwd=3, lty=2);abline(v = quantiles2, col="blue", lwd=3)
  mtext(paste0('True=',0.0369, ", 95% CI= ", paste(round(quantiles2,4), collapse = '-')), 3, line=0)

  quantiles1<-getBST_CI(LM_measured_reps[,1],replicates);quantiles2<-getBST_CI(LM_measured_reps[,2],replicates)
  hist(LM_measured_reps[,1], col="lightblue", main = 'Intercept', xlim=ran_in);abline(v = 0.268, col="red", lwd=3, lty=2);abline(v = quantiles1, col="blue", lwd=3)
  mtext(paste0('True=',0.268, ", 95% CI= ", paste(round(quantiles1,4), collapse = '-')), 3, line=0)
  hist(LM_measured_reps[,2], col="lightblue", main = 'Slope', xlim=ran_sl);abline(v = 0.0369, col="red", lwd=3, lty=2);abline(v = quantiles2, col="blue", lwd=3)
  mtext(paste0('True=',0.0369, ", 95% CI= ", paste(round(quantiles2,4), collapse = '-')), 3, line=0)
  
  quantiles1<-getBST_CI(BLM_measured[[1]][,1],replicates);quantiles2<-getBST_CI(BLM_measured[[1]][,2],replicates)
  hist(BLM_measured[[1]][,1], col="lightblue", main = 'Intercept', xlab='Measured BLM with errors', xlim=ran_in);abline(v = 0.268, col="red", lwd=3, lty=2);abline(v = quantiles1, col="blue", lwd=3)
  mtext(paste0('True=',0.268, ", 95% CI= ", paste(round(quantiles1,4), collapse = '-')), 3, line=0)
  hist(BLM_measured[[1]][,2], col="lightblue", main = 'Slope', xlab='BLM with errors', xlim=ran_sl);abline(v = 0.0369, col="red", lwd=3, lty=2);abline(v = quantiles2, col="blue", lwd=3)
  mtext(paste0('True=',0.0369, ", 95% CI= ", paste(round(quantiles2,4), collapse = '-')), 3, line=0)
  
  quantiles1<-getBST_CI(BLM_measured[[2]][,1],replicates);quantiles2<-getBST_CI(BLM_measured[[2]][,2],replicates)
  hist(BLM_measured[[2]][,1], col="lightblue", main = 'Intercept', xlab='Measured BLM no errors', xlim=ran_in);abline(v = 0.268, col="red", lwd=3, lty=2);abline(v = quantiles1, col="blue", lwd=3)
  mtext(paste0('True=',0.268, ", 95% CI= ", paste(round(quantiles1,4), collapse = '-')), 3, line=0)
  hist(BLM_measured[[2]][,2], col="lightblue", main = 'Slope', xlab='BLM no errors', xlim=ran_sl);abline(v = 0.0369, col="red", lwd=3, lty=2);abline(v = quantiles2, col="blue", lwd=3)
  mtext(paste0('True=',0.0369, ", 95% CI= ", paste(round(quantiles2,4), collapse = '-')), 3, line=0)
  
  dev.off()
  
  
  
  pdf(paste0(name,'_replicates=',replicates,'_samples=',samples, '_generations=',generations ,'Raw.pdf' ), 7,7)
  par(mfrow=c(4,2))
  
  quantiles1<-getBST_CI(York_measured_reps[,1],replicates);quantiles2<-getBST_CI(York_measured_reps[,2],replicates)
  hist(York_measured_reps[,1], col="lightblue", main = 'Intercept');abline(v = 0.268, col="red", lwd=3, lty=2);abline(v = quantiles1, col="blue", lwd=3)
  mtext(paste0('True=',0.268, ", 95% CI= ", paste(round(quantiles1,4), collapse = '-')), 3, line=0)
  hist(York_measured_reps[,2], col="lightblue", main = 'Slope');abline(v = 0.0369, col="red", lwd=3, lty=2);abline(v = quantiles2, col="blue", lwd=3)
  mtext(paste0('True=',0.0369, ", 95% CI= ", paste(round(quantiles2,4), collapse = '-')), 3, line=0)
  

  quantiles1<-getBST_CI(LM_measured_reps[,1],replicates);quantiles2<-getBST_CI(LM_measured_reps[,2],replicates)
  hist(LM_measured_reps[,1], col="lightblue", main = 'Intercept');abline(v = 0.268, col="red", lwd=3, lty=2);abline(v = quantiles1, col="blue", lwd=3)
  mtext(paste0('True=',0.268, ", 95% CI= ", paste(round(quantiles1,4), collapse = '-')), 3, line=0)
  hist(LM_measured_reps[,2], col="lightblue", main = 'Slope');abline(v = 0.0369, col="red", lwd=3, lty=2);abline(v = quantiles2, col="blue", lwd=3)
  mtext(paste0('True=',0.0369, ", 95% CI= ", paste(round(quantiles2,4), collapse = '-')), 3, line=0)
  
  quantiles1<-getBST_CI(BLM_measured[[1]][,1],replicates);quantiles2<-getBST_CI(BLM_measured[[1]][,2],replicates)
  hist(BLM_measured[[1]][,1], col="lightblue", main = 'Intercept', xlab='Measured BLM with errors');abline(v = 0.268, col="red", lwd=3, lty=2);abline(v = quantiles1, col="blue", lwd=3)
  mtext(paste0('True=',0.268, ", 95% CI= ", paste(round(quantiles1,4), collapse = '-')), 3, line=0)
  hist(BLM_measured[[1]][,2], col="lightblue", main = 'Slope', xlab='Measured BLM with errors');abline(v = 0.0369, col="red", lwd=3, lty=2);abline(v = quantiles2, col="blue", lwd=3)
  mtext(paste0('True=',0.0369, ", 95% CI= ", paste(round(quantiles2,4), collapse = '-')), 3, line=0)
  
  quantiles1<-getBST_CI(BLM_measured[[2]][,1],replicates);quantiles2<-getBST_CI(BLM_measured[[2]][,2],replicates)
  hist(BLM_measured[[2]][,1], col="lightblue", main = 'Intercept', xlab='Measured BLM no errors');abline(v = 0.268, col="red", lwd=3, lty=2);abline(v = quantiles1, col="blue", lwd=3)
  mtext(paste0('True=',0.268, ", 95% CI= ", paste(round(quantiles1,4), collapse = '-')), 3, line=0)
  hist(BLM_measured[[2]][,2], col="lightblue", main = 'Slope', xlab='Measured BLM no errors');abline(v = 0.0369, col="red", lwd=3, lty=2);abline(v = quantiles2, col="blue", lwd=3)
  mtext(paste0('True=',0.0369, ", 95% CI= ", paste(round(quantiles2,4), collapse = '-')), 3, line=0)
  
 
  
  dev.off()
  
  return(list('York_measured_reps'=York_measured_reps,
              'LM_measured_reps'=LM_measured_reps,
              'BLM_measured'=BLM_measured
  ))
  
}



SimulateRegressions<-function(universe,nobs=1000 ,wd=getwd(), folderName=Sys.Date()){
  
  
  dir.create(file.path(wd, folderName), showWarnings = FALSE)
  setwd(file.path(wd, folderName))
  
  fullDataSimulations<-pblapply(seq_along(universe[,1]), function(x){
    
    ##Generate name
    name = paste0("y_error=",universe$y_error[x], 
                  '_y_measurement_error=', universe$y_measurement_error[x],
                  "_x_error=",universe$x_error[x], 
                  "_sampled=",universe$samples[x],
                  "_replicates=", universe$replicates[x], 
                  "_generations=", universe$generations[x])
    
    dir.create(file.path(wd,folderName, name), showWarnings = FALSE)
    setwd(file.path(wd, folderName, name))
    
    #Generate data
    
    data<-sim_slr(nobs=1000, sdobsx = universe$x_error[x], meanX = 12.02585352, alpha=0.268, beta=0.0369, sdy=universe$y_error[x], sdobsy = universe$y_measurement_error[x])
    
    
    q<-ggplot(data = data)+ 
      geom_errorbar(aes(x = x_measured, ymin = y_measured-y_SD,ymax = y_measured+y_SD), color='grey')+
      geom_errorbarh(aes(y=y_measured,xmin = x_measured-x_SD,xmax = x_measured+x_SD), color='grey')+
      geom_point(aes(x = x_measured,y = y_measured))+
      xlab('Temperature (10^6/T^2)')+ylab("D47")+ggtitle('Measurement error in X')+ 
      geom_abline(intercept=0.268, slope=0.0369)+ 
      geom_abline(intercept=0.268, slope=0.0369, color='blue')
    
    pdf(paste0(name, '_Data.pdf'), 7,5)
    print(q)
    dev.off()
    
    #Run the models
    sal<-simulateAll(data, error=1e-5, replicates=universe$replicates[x], samples=universe$samples[x], generations=universe$generations[x], name=name)

    full_replicates<-rbind.data.frame( 
      rbindlist(sal[-c(3)], idcol=TRUE),
      rbindlist(sal[[3]], idcol=TRUE))
    
    write.csv(full_replicates, paste0(name, 'replicates.csv'))
    write.csv(data, paste0(name, 'data.csv'))

    setwd(file.path(wd, folderName))
    
    
  })
  setwd(wd)
  
  names(fullDataSimulations)<-unlist(lapply(1:nrow(universe), function(x){ paste0("y_error=",universe$y_error[x], 
                                                                                  "_x_error=",universe$x_error[x], 
                                                                                  "_sampled=",universe$samples[x],
                                                                                  "_replicates=", universe$replicates[x], 
                                                                                  "_generations=", universe$generations[x])}))
  
  
  fullDataSimulations
  
}


SimulateRegressions_dataGiven<-function(universe, X, X_un){
  
  fullDataSimulations<-pblapply(seq_along(universe[,1]), function(x){
    
    ##Generate name
    name = paste0("y_error=",universe$y_error[x], 
                  "_sampled=",universe$samples[x],
                  "_replicates=", universe$replicates[x], 
                  "_generations=", universe$generations[x],
                  '_dataGiven')
    
    #Generate data
    
    data=sim_slr_givenX(x=X, x_un=X_un, alpha=0.268, beta=0.0369, sdy=universe$y_error[x])
    
    
    #Plot the data
   
    q<-ggplot(data = data)+ 
      geom_errorbar(aes(x = x_measured, ymin = y_measured-abs(y_SD),ymax = y_measured+abs(y_SD)), color='grey')+
      geom_errorbarh(aes(y=y_measured,xmin = x_measured-x_SD,xmax = x_measured+x_SD), color='grey')+
      geom_point(aes(x = x_measured,y = y_measured))+
      xlab('Temperature (10^2/T^2)')+ylab("D47")+ggtitle('Measurement error in X')+ 
      geom_abline(intercept=0.268, slope=0.0369) 

    pdf(paste0(name, '_Data.pdf'), 7,5)
    print(q)
    dev.off()
    
    #Run the models
    sal<-simulateAll_dataGiven(data, error=1e-5, replicates=universe$replicates[x], samples=universe$samples[x], generations=universe$generations[x], name=name)
    
    full_replicates<-rbind.data.frame( 
      rbindlist(sal[-c(3)], idcol=TRUE),
      rbindlist(sal[[3]], idcol=TRUE))
    write.csv(full_replicates, paste0(name, '.csv'))
    
  })
  
  
  names(fullDataSimulations)<-unlist(lapply(1:nrow(universe), function(x){ paste0("y_error=",universe$y_error[x], 
                                                                                  "_x_error=",universe$x_error[x], 
                                                                                  "_sampled=",universe$samples[x],
                                                                                  "_replicates=", universe$replicates[x], 
                                                                                  "_generations=", universe$generations[x])}))
  
  
  fullDataSimulations
  
}



SimulateRegressionsMixed<-function(universe, nobs=100 ,wd=getwd(), folderName=Sys.Date(), alpha=c(0.268,0.268), beta=c(0.0369,0.0369)){
  
  
  dir.create(file.path(wd, folderName), showWarnings = FALSE)
  setwd(file.path(wd, folderName))
  
  fullDataSimulations<-pblapply(seq_along(universe[,1]), function(x){
    tryCatch({
    ##Generate name
    name = paste0("y_error=",universe$y_error[x], 
                  '_y_measurement_error=', universe$y_measurement_error[x],
                  "_x_error=",universe$x_error[x], 
                  "_sampled=",universe$samples[x],
                  "_replicates=", universe$replicates[x], 
                  "_generations=", universe$generations[x])
    
    dir.create(file.path(wd,folderName, name), showWarnings = FALSE)
    setwd(file.path(wd, folderName, name))
    
    #Generate data
        data<-sim_slr_mult(nobs=nobs, sdobsx = universe$x_error[x], meanX = 12.02585352, alpha=alpha , beta=beta, sdy=universe$y_error[x], sdobsy = universe$y_measurement_error[x])
        #Plot the data
 
        q<-ggplot(data = data)+ 
          geom_errorbar(aes(x = x_measured, ymin = y_measured-y_SD,ymax = y_measured+y_SD), color='grey')+
          geom_errorbarh(aes(y=y_measured,xmin = x_measured-x_SD,xmax = x_measured+x_SD), color='grey')+
          geom_point(aes(x = x_TRUE,y = y_TRUE),color='red')+
          geom_point(aes(x = x_measured,y = y_measured)) + 
          xlab('Temperature (10^6/T^2)')+ylab("D47")+ggtitle('Measurement error in X')+ 
          geom_abline(intercept=alpha[1], slope=beta[1], color='blue')+
          geom_abline(intercept=alpha[2], slope=beta[2], color='red')+theme_bw()
        
        
        
        
        pdf(paste0(name, '_Data.pdf'), 7,5)
        print(q)
        dev.off()
        
        sal<-simulateAll(data, error=1e-5, 
                              replicates=universe$replicates[x], 
                              samples=universe$samples[x], 
                              generations=universe$generations[x],isMixed=T)
  
    #Run the models

        full_replicates<-rbind.data.frame( 
          cbind.data.frame(rbindlist(sal[-c(3)], idcol=TRUE),material=NA),
          cbind.data.frame(rbindlist(sal[[3]][1:2], idcol=TRUE), material=NA),
          cbind.data.frame(.id= 'BLMM_Measured_errors', sal[[3]][[3]]))
      
        
    write.csv(full_replicates, paste0(name, 'replicates.csv'))
    write.csv(data, paste0(name, 'data.csv'))
    
    setwd(file.path(wd, folderName))
    }, error=function(e){})
  })
  setwd(wd)
  
  names(fullDataSimulations)<-unlist(lapply(1:nrow(universe), function(x){ paste0("y_error=",universe$y_error[x], 
                                                                              "_generations=", universe$generations[x])}))
  
  
  fullDataSimulations
  
}



