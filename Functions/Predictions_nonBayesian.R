#' This function performs temp reconstruction (10^6/T^2 with T in K) for
#' multiple replicates of the same target.
#' 
#' @param calData Calibration dataset
#' @param recData Reconstruction dataset
#' @param obCal Data.frame summarizing the distribution of slopes and intercepts
#' @clumpedClassic Whether the classic calibration or simple inversion should be conducted

predictTc <<- function(calData,
                       recData,
                       obCal){
  
    temp <- sqrt((mean(obCal$beta) * 10 ^ 6) / 
                   (recData$D47 - mean(obCal$alpha))) - 273.15
    temp_SE <- sqrt((mean(obCal$beta) * 10 ^ 6) / 
                      (recData$D47 + recData$D47error - mean(obCal$alpha))) - 273.15
    error = (temp-temp_SE)

    recTempS <- cbind.data.frame(Sample = recData$Sample,
                                 D47 = recData$D47,
                                 D47error = recData$D47error,
                                 meanTemp = temp, 
                                 error = error)
    
  
}


