#root mean squared error
RMSE <- function(predicted, observed){
  sqrt(mean((predicted - observed)^2, na.rm=T))
}

std <- function(x) sd(x, na.rm = T)/sqrt(length(x))

#Error range overlap
REO <- function(predictedlow, predictedupper ,observedMean, observedSE){
  
 a<- do.call(rbind.data.frame,
  lapply(seq_along(predictedlow), function(x){
    
    Pint<-c(predictedlow,predictedupper)
    Oint<-c(observedMean[x]-1.96*observedSE[x], observedMean[x]+1.96*observedSE[x])
    
    ##If they are nested
    ON<-Oint[1] > Pint[1] & Oint[2] < Pint[2] #Observed is nested
    PN<-Pint[1] > Oint[1] & Pint[2] < Oint[2] #Predicted is nested
    Inside<-if(ON | PN){ifelse(isTRUE(ON) , 'ObservedRange','PredictedRange')}else{ NA }
    if(!is.na(Inside) ){ PartialOverlap<-NA; DistanceNoOverlap<-NA }
    
    ##If they do not overlap
    if(is.na(Inside) ){
      
      PHi<-Oint[2] <  Pint[1] #Predicted is higher
      ObsHi<-Pint[2] <  Oint[1] #Observed is higher
      
      if(any( c(PHi,ObsHi))){
        DistanceNoOverlap<-if(PHi){ Pint[1] - Oint[2] }else{ Pint[2]-Oint[1] }
        PartialOverlap<-NA
        Inside<-NA
      }else{
        ##If they partially overlap
        PpartiallyHi<- Pint[1] < Oint[2]   #Predicted is partially higher
        ObspartiallyHi<- Oint[1] < Pint[2]   #Observed is partially higher
        
        DistanceNoOverlap<-NA
        Inside<-NA
        PartialOverlap<-if(PpartiallyHi){ Oint[2] - Pint[1] }else{ Pint[2]-Oint[1] }
        
      }  
    }
    
    data.frame(Inside, PartialOverlap, DistanceNoOverlap)

  }))
  
 tg<-which(names(table(a$Inside)) == 'PredictedRange' )
 tg2<- length(which(names(table(a$Inside)) == 'ObservedRange' ))
 
 cbind.data.frame(
 #% of predictions within the observed range
 '%PredWithinObsRange' =as.data.frame(table(a$Inside)/nrow(a))[tg,2]*100,
 #% of observations within the predicted range
 '%ObsWithinPredRange' = if( tg2!=0 ){ as.data.frame(table(a$Inside)/nrow(a))[1,2]*100}else{NA },
 #% of predictions partially overlapped
 '%PartOverlap' = (length(na.omit(a$PartialOverlap))/nrow(a))*100,
 #% of predictions no overlap
 '%NoOverlap' = (length(na.omit(a$DistanceNoOverlap))/nrow(a))*100,
 
 #mean partial overlap
 '%MeanPartialOverlap' = (mean(na.omit(a$PartialOverlap))),
 '%SEPartialOverlap' =(std(na.omit(a$PartialOverlap))),
 '%MeanNoOverlap' =(mean(na.omit(a$DistanceNoOverlap))),
 '%SENoOverlap' = (std(na.omit(a$DistanceNoOverlap)))
 )


}

#mean absolute error
MAE <- function(predicted, observed)
{
  mean(abs(predicted - observed), na.rm=T)
}


confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector, na.rm = T)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector, na.rm = T)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- cbind.data.frame("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}
