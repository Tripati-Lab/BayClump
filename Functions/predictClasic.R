library(investr)

#' This function performs temp reconstruction (10^6/T^2 with T in K) for
#' multiple replicates of the same target.
#' 
#' @param targety independent measurements of D47
#' @param calData Calibration data with the following columns: D47, T2, 
#'                TempError, D47error
#' @param model string: lm, wlm, Deming, York


predictTclassic <<- function(calData, targety, model='lm', replicates=1000, DegreeC=T, bootDataset=T, onlyMedian=T){
  calData$T2 <- calData$Temperature
  
  if(model == 'lm'){
    
    if(bootDataset){
    res <<- do.call(rbind, lapply(1:replicates, function(x){
    mod <<- lm(D47 ~ T2, calData[sample(1:nrow(calData), nrow(calData), replace = T),])
    estimate <<- investr::invest(mod, y0 = targety,
           interval = "percentile", 
           nsim = if(bootDataset){1}else{replicates}, seed = 3, 
           extendInt="yes", progress=F, 
           lower=-100,
           upper=100)
    }))
    
    }else{
      
      mod <<- lm(D47 ~ T2, calData)
      estimate <<- investr::invest(mod, y0 = targety,
                                   interval = "percentile", 
                                   nsim = if(bootDataset){1}else{replicates}, seed = 3, 
                                   extendInt="yes", progress=F, 
                                   lower=-100,
                                   upper=100)
      
    }
    
  }
  
  if(model == 'wlm'){
    
    if(bootDataset){
      res <- do.call(rbind, lapply(1:replicates, function(x){
        ds <<- calData[sample(1:nrow(calData), nrow(calData), replace = T),]
        Reg0 <- lm(D47 ~ T2,  ds)
        wt <- 1 / lm(abs(Reg0$residuals) ~ Reg0$fitted.values)$fitted.values^2
        mod<<-lm(D47 ~ T2,  ds, weights = wt)
        estimate <<- investr::invest(mod, y0 = targety,
                                     interval = "percentile", 
                                     nsim = if(bootDataset){1}else{replicates}, seed = 3, 
                                     extendInt="yes", progress=F, 
                                     lower=-100,
                                     upper=100)
      }))
      
    }else{
    
    Reg0 <<- lm(D47 ~ T2,  calData)
    wt <<- 1 / lm(abs(Reg0$residuals) ~ Reg0$fitted.values)$fitted.values^2
    mod <<-lm(D47 ~ T2,  calData, weights = wt)
    estimate <<- investr::invest(mod, y0 = targety,
                       interval = "percentile", 
                       nsim = if(bootDataset){1}else{replicates}, seed = 3, 
                       extendInt="yes", progress=F, 
                       lower=-100,
                       upper=100)
    }
  }
  
  if(model == 'York'){
    
    if(bootDataset){
      res <<- do.call(rbind, lapply(1:replicates, function(x){
        
        dat <<- calData[sample(1:nrow(calData), nrow(calData), replace = T),]
        
        fit<<-york(cbind.data.frame(dat$T2, dat$TempError, dat$D47, dat$D47error))
        fit1 <<- stats::nls(formula = D47  ~ a + b1*T2,
                            data = dat, 
                            start = list(a = fit$a[1],b1 = fit$b[1]),
                            lower = c(a = (fit$a[1]-fit$a[2]*1.96),b1 = (fit$b[1]-fit$b[2]*1.96)),
                            upper = c(a = (fit$a[1]+fit$a[2]*1.96),b1 = (fit$b[1]+fit$b[2]*1.96)),
                            algorithm = "port") 
        estimate <<- investr::invest(fit1, y0 = targety,
                                     interval = "percentile", 
                                     nsim = if(bootDataset){1}else{replicates}, seed = 3, 
                                     extendInt="yes", progress=F, 
                                     lower=-100,
                                     upper=100)
        
      }))
      
    }else{
    
    fit<<-york(cbind.data.frame(calData$T2, calData$TempError, calData$D47, calData$D47error))
    fit1 <<- stats::nls(formula = D47  ~ a + b1*T2,
                data = calData, 
                start = list(a = fit$a[1],b1 = fit$b[1]),
                lower = c(a = (fit$a[1]-fit$a[2]*1.96),b1 = (fit$b[1]-fit$b[2]*1.96)),
                upper = c(a = (fit$a[1]+fit$a[2]*1.96),b1 = (fit$b[1]+fit$b[2]*1.96)),
                algorithm = "port") 
    estimate <<- investr::invest(fit1, y0 = targety,
                       interval = "percentile", 
                       nsim = if(bootDataset){1}else{replicates}, seed = 3, 
                       extendInt="yes", progress=F, 
                       lower=-100,
                       upper=100)
    }
  }
  
  if(model == 'Deming'){
    
    
    
    if(bootDataset){
      res <<- do.call(rbind, lapply(1:replicates, function(x){
        
        dat <<- calData[sample(1:nrow(calData), nrow(calData), replace = T),]
        
        fit<<-deming(D47 ~ T2, dat, xstd= 1/TempError^2, ystd= 1/D47error^2)
        fit1 <<- nls(formula = D47  ~ a + b1*T2,
                     data = dat, 
                     start = list(a = fit$coefficients[1],b1 = fit$coefficients[2]),
                     lower = c(a = (fit$ci[1,1]),b1 = (fit$ci[2,1])),
                     upper = c(a = (fit$ci[1,2]),b1 = (fit$ci[2,2])),
                     algorithm = "port") 
        
        estimate <<- investr::invest(fit1, y0 = targety,
                                     interval = "percentile", 
                                     nsim = if(bootDataset){1}else{replicates}, seed = 3, 
                                     extendInt="yes", progress=F, 
                                     lower=-100,
                                     upper=100)
      }))
      
    }else{
    
    fit<<-deming(D47 ~ T2, calData, xstd= 1/TempError^2, ystd= 1/D47error^2)
    fit1 <<- nls(formula = D47  ~ a + b1*T2,
                data = calData, 
                start = list(a = fit$coefficients[1],b1 = fit$coefficients[2]),
                lower = c(a = (fit$ci[1,1]),b1 = (fit$ci[2,1])),
                upper = c(a = (fit$ci[1,2]),b1 = (fit$ci[2,2])),
                algorithm = "port") 
    
    estimate <<- investr::invest(fit1, y0 = targety,
                       interval = "percentile", 
                       nsim = if(bootDataset){1}else{replicates}, seed = 3, 
                       extendInt="yes", progress=F, 
                       lower=-100,
                       upper=100)
    }
  }
  
 
  
  
  if(bootDataset){
    res<-as.data.frame(res)
    res[,1:5] <- lapply(res[,1:5], as.numeric)
    se <- sd(res$estimate)/sqrt(length(res$estimate))

    if(onlyMedian){
    preds <- cbind.data.frame(D47= mean( targety ), 
                              D47se=sd(targety)/sqrt(length(targety)),
                              temp=median(res$estimate), se=se, lwr=mean(res$estimate)-1.96*se, upr=mean(res$estimate)+1.96*se
    )
    }else{
      preds <- cbind.data.frame(D47= mean( targety ), 
                                D47se=sd(targety)/sqrt(length(targety)),
                                temp=median(res$estimate), se=median(res$se), lwr=median(res$lower), upr=mean(res$lower)
      )
    }
    
    
  }else{
    preds <- cbind.data.frame(D47= mean( targety ), 
                              D47se=sd(targety)/sqrt(length(targety)),
                              temp=estimate$estimate, se=estimate$se, lwr=estimate$lower, upr=estimate$upper
    )
  }
    if(DegreeC){
      
      predsUp<- sqrt(10^6/(preds$temp+preds$se))-273.15
      predsPoint <- sqrt(10^6/preds$temp)-273.15
      
      preds <- cbind.data.frame(D47= preds$D47, 
                                D47se= preds$D47se, 
                                temp=predsPoint, 
                                se=predsPoint-predsUp, 
                                lwr=predsPoint-((predsPoint-predsUp)*1.96), 
                                upr=predsPoint+((predsPoint-predsUp)*1.96)
      )
      
      return(preds)
    }else{
      return(preds)
    }
    
}

##Predictions based on the actual replicates

classicCalibration <- function(reps, targetD47, error_targetD47, material,  mixed=F) {
  
  if(mixed){
   do.call(rbind, lapply(1:length(targetD47), function(x){
    point <- sqrt((median(reps$slope[reps$material == material[x]]) * 10 ^ 6) / 
                    (targetD47[x] -  median(reps$intercept[reps$material == material[x]]))) - 273.15
      
    error_point <- point - (sqrt((median(reps$slope[reps$material == material[x]]) * 10 ^ 6) / (targetD47[x] + error_targetD47[x] - median(reps$intercept[reps$material == material[x]]))) - 273.15)
    cbind.data.frame(targetD47=targetD47[x],error_targetD47=error_targetD47[x], material=material[x], Tc=point, se=error_point)
    
   
    })
   )
    
  }else{
    point <- sqrt((median(reps$slope) * 10^6) / (targetD47 -  median(reps$intercept))) - 273.15
    error_point <- point - (sqrt((median(reps$slope) * 10 ^ 6) / (targetD47 + error_targetD47  - median(reps$intercept))) - 273.15)
    cbind.data.frame(targetD47=targetD47,error_targetD47=error_targetD47, Tc=point, se=error_point)
  }
}