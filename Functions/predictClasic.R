library(investr)

#' This function performs temp reconstruction (10^6/T^2 with T in K) for
#' multiple replicates of the same target.
#' 
#' @param targety independent measurements of D47
#' @param calData Calibration data with the following columns: D47, T2, 
#'                TempError, D47error
#' @param model string: lm, wlm, Deming, York


predictTclassic<-function(calData, targety, model='lm'){
  
  if("T2" %in% colnames(calData) ){}else{
  calData$T2 <- calData$Temperature
  }
  
  if(model == 'lm'){
    mod <- lm(D47 ~ T2, calData)
    estimate <- invest(mod, y0 = targety,
           interval = "percentile", 
           nsim = 1000, seed = 3, 
           extendInt="yes", progress=T, 
           lower=-100,
           upper=100)
  }
  
  if(model == 'wlm'){
    mod<-lm(D47 ~ T2,  calData, weights = D47error)
    estimate <- invest(mod, y0 = targety,
                       interval = "percentile", 
                       nsim = 1000, seed = 3, 
                       extendInt="yes", progress=T, 
                       lower=-100,
                       upper=100)
  }
  
  if(model == 'York'){
    fit<-york(cbind.data.frame(calData$T2, calData$TempError, calData$D47, calData$D47error))
    fit1 <- nls(formula = D47  ~ a + b1*T2,
                data = calData, 
                start = list(a = fit$a[1],b1 = fit$b[1]),
                lower = c(a = (fit$a[1]-fit$a[2]*1.96),b1 = (fit$b[1]-fit$b[2]*1.96)),
                upper = c(a = (fit$a[1]+fit$a[2]*1.96),b1 = (fit$b[1]+fit$b[2]*1.96)),
                algorithm = "port") 

    estimate <- invest(fit1, y0 = targety,
                       interval = "percentile", 
                       nsim = 1000, seed = 3, 
                       extendInt="yes", progress=T, 
                       lower=-100,
                       upper=100)
  }
  
  if(model == 'Deming'){
    fit<-deming(D47 ~ Temperature, calData, xstd= calData$TempError, ystd= calData$D47error)
    fit1 <- nls(formula = D47  ~ a + b1*T2,
                data = calData, 
                start = list(a = fit$coefficients[1],b1 = fit$coefficients[2]),
                lower = c(a = (fit$ci[1,1]),b1 = (fit$ci[2,1])),
                upper = c(a = (fit$ci[1,2]),b1 = (fit$ci[2,2])),
                algorithm = "port") 
    
    estimate <- invest(fit1, y0 = targety,
                       interval = "percentile", 
                       nsim = 1000, seed = 3, 
                       extendInt="yes", progress=T, 
                       lower=-100,
                       upper=100)
  }
  
  cbind.data.frame(D47= paste( targety, collapse = "," ), 
                   temp=estimate$estimate, lwr=estimate$lower, upr=estimate$upper
                   )

}


