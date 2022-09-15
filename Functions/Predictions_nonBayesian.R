#' This function performs temp reconstruction (10^6/T^2 with T in K) for
#' multiple replicates of the same target.
#' 
#' @param calData Calibration dataset
#' @param recData Reconstruction dataset
#' @param obCal Data.frame summarizing the distribution of slopes and intercepts
#' @clumpedClassic Whether the classic calibration or simple inversion should be conducted

predictTc <<- function(calData,
                       recData,
                       obCal,
                       clumpedClassic = TRUE){
  
  
  if(clumpedClassic){
    
    temp <- sqrt((mean(obCal$beta) * 10 ^ 6) / 
                   (recData$D47 - mean(obCal$alpha))) - 273.15
    temp_SE <- sqrt((mean(obCal$beta) * 10 ^ 6) / 
                      (recData$D47 + recData$D47error - mean(obCal$alpha))) - 273.15
    se = temp-temp_SE
    recTempS <- cbind.data.frame(Sample = recData$Sample,
                                 D47 = recData$D47,
                                 D47error = recData$D47error,
                                 meanTemp = temp, 
                                 Temp_L = temp - 1.96 * se, 
                                 Temp_H = temp + 1.96 * se)
    
  }else{
    
    #McClelland, H. L., Halevy, I., Wolf‐Gladrow, D. A., Evans, D., & Bradley, A. S. (2021). Statistical uncertainty in paleoclimate proxy reconstructions. Geophysical Research Letters, 48(15), e2021GL092773.
    
    X_new_hat  <- (recData$D47 - mean(obCal$alpha)) / mean(obCal$beta)
    
    Cqts   <- c(0.025, 0.975)
    Tval   <- qt(Cqts, df = (nrow(calData)-2)) 
    Yhat   <-  mean(obCal$alpha) + calData$Temperature*mean(obCal$beta) 
    ysd    <- sd(calData$D47 - Yhat)
    Xe1    <- (ysd * Tval) / mean(obCal$beta)
    Xe1L   <- sqrt(10^6/ (X_new_hat + Xe1[1])) - 273.15
    Xe1U   <- sqrt(10^6/(X_new_hat + Xe1[2])) - 273.15 
    
    cbind.data.frame(Sample = recData$Sample, 
                     D47 = recData$D47, 
                     D47error = recData$D47error, 
                     meanTemp = sqrt(10^6/X_new_hat)-273.15, 
                     Temp_L = Xe1U, 
                     Temp_H = Xe1L)
  }
  
}


