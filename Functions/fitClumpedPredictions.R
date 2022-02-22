fitClumpedPredictions<-function(calibrationData, 
                                n.iter= 5000, 
                                burninFrac=0.5,
                                priors = "informative",
                                D47error='D47error', 
                                D47Pred,
                                D47Prederror){
  
  
  if(priors == 'informative'){
    alphaBLM1='dnorm(0.231,0.065)' 
    betaBLM1= 'dnorm(0.039,0.004)'}else{
      alphaBLM1='dnorm(0, 0.01)' 
      betaBLM1= 'dnorm(0, 0.01)'
    }
  
  
  ##Models
  BLM1<-paste(" model{
    # Diffuse normal priors for predictors
    alpha ~ ", alphaBLM1," \n ",
              "beta ~ ", betaBLM1," \n ", 
              "
    sigma <- 1/sqrt(tauy)                              
    tauy ~ dgamma(0.1, 0.1)                                
    
    for (i in 1:N){
        x[i] ~ dnorm(11,0.01)
    }
    # Likelihood
    for (i in 1:N){
        obsy[i] ~ dnorm(y[i],pow(erry[i],-2))
        y[i] ~ dnorm(mu[i],tauy)
        obsx[i] ~ dnorm(x[i],pow(errx[i],-2))
        mu[i] <- alpha + beta*x[i]
    }
}")
  
  LM_Data <- list(obsx = calibrationData$Temperature , obsy = calibrationData$D47 , 
                  errx = calibrationData$TempError, erry = calibrationData[,D47error], 
                  N=nrow(calibrationData))
  BLM1_fit <- jags(data = LM_Data, inits = NULL,
                   parameters = c("alpha","beta", "tauy"),
                   model = textConnection(BLM1), n.chains = 3, 
                   n.iter = n.iter, n.burnin = n.iter*burninFrac)
  
  
  BLM1Pred<-paste("model{

for(i in 1:N)
{
obsy1[i] ~ dnorm(y1[i],erry[i])
y1[i]~dnorm(mup[i],tauy)
mup[i] <-  beta*Tcpropagated[i] + alpha
Tcpropagated[i]~dnorm(sqrt(10^6/T0[i])-273.15, 35)
}

for(i in 1:N)
{
obsy2[i] ~ dnorm(y2[i],erry[i])
y2[i]~dnorm(mup[i],tauy)
mup2[i] <-  beta*Tcpropagated2[i] + alpha
Tcpropagated2[i]~dnorm(sqrt(10^6/T1[i])-273.15, 35)
}

Terr <- Tcpropagated-Tcpropagated2

}")
  
  sa <- nrow(BLM1_fit$BUGSoutput$sims.matrix)
  
  valsPreds <- lapply(1:ifelse(sa>500, 500), function(i){
    
    T0=(D47Pred- BLM1_fit$BUGSoutput$sims.matrix[i,1])/
      BLM1_fit$BUGSoutput$sims.matrix[i,2]
    T0 <- ifelse(T0 == 0, 11, T0)
    
    T1=(D47Pred+D47Prederror- BLM1_fit$BUGSoutput$sims.matrix[i,1])/
      BLM1_fit$BUGSoutput$sims.matrix[i,2]
    T1 <- ifelse(T1 == 0, 11, T1)
    
  LM_Data_pred <- list(N=length(D47Pred), obsy1 = D47Pred, obsy2=D47Pred+D47Prederror,
                       erry=D47Prederror,
                       alpha= BLM1_fit$BUGSoutput$sims.matrix[i,1],
                       beta=BLM1_fit$BUGSoutput$sims.matrix[i,2],
                       tauy=BLM1_fit$BUGSoutput$sims.matrix[i,4],
                       T0=T0,
                       T1=T1
  )
  
   BLM1_fit_pred <- jags(data = LM_Data_pred,inits = NULL,
                        parameters = c("Tcpropagated","Terr"),
                        model = textConnection(BLM1Pred), n.chains = 3, 
                        n.iter = n.iter,  n.burnin = n.iter*burninFrac)
  cbind.data.frame(D47Pred,D47Prederror, Temp=BLM1_fit_pred$BUGSoutput$mean$Tcpropagated,
  SE=BLM1_fit_pred$BUGSoutput$mean$Terr)
  })
  
  valsPreds <- do.call(rbind, valsPreds)
  
  Preds <- valsPreds %>% 
    group_by(D47Pred, D47Prederror) %>% 
    summarise(Temp = mean(Temp),
              SE = mean(SE), .groups = 'drop')
  
  BLM1_NoErrors<-paste("model{
                # Diffuse normal priors for predictors
                alpha ~ ", alphaBLM1," \n ",
                       "beta ~ ", betaBLM1," \n ",
                       "
    sigma <- 1/sqrt(tau)                              
    tau ~ dgamma(0.1, 0.1)  
                for (i in 1:N){
                y[i] ~ dnorm(mu[i],tau)
                mu[i]<- eta[i]
                eta[i] <- alpha + inprod(x[i],beta)
  }
}")
 
    LM_No_error_Data <- list(x = calibrationData$Temperature , y = calibrationData$D47,
                             N=nrow(calibrationData))
   
    BLM1_fit_NoErrors <- jags(data = LM_No_error_Data,inits = NULL,
                              parameters = c("alpha","beta", "tau"),
                              model = textConnection(BLM1_NoErrors), n.chains = 3,
                              n.iter = n.iter,  n.burnin = n.iter*burninFrac)

    
    
    sa <- nrow(BLM1_fit_NoErrors$BUGSoutput$sims.matrix)
    valsPreds <- lapply(sample(1:ifelse(sa>500,500,sa),ifelse(sa>500,500,sa)), function(i){
      
      
      T0=(D47Pred- BLM1_fit_NoErrors$BUGSoutput$sims.matrix[i,1])/
        BLM1_fit_NoErrors$BUGSoutput$sims.matrix[i,2]
      T0 = ifelse(T0==0, 11, T0)
      
      T1=(D47Pred+D47Prederror- BLM1_fit_NoErrors$BUGSoutput$sims.matrix[i,1])/
        BLM1_fit_NoErrors$BUGSoutput$sims.matrix[i,2]
      T1 = ifelse(T1==0, 11, T1)
      
      LM_Data_pred <- list(N=length(D47Pred), obsy1 = D47Pred, obsy2=D47Pred+D47Prederror,
                           erry=D47Prederror,
                           alpha= BLM1_fit_NoErrors$BUGSoutput$sims.matrix[i,1],
                           beta=BLM1_fit_NoErrors$BUGSoutput$sims.matrix[i,2],
                           tauy=BLM1_fit_NoErrors$BUGSoutput$sims.matrix[i,4],
                           T0=T0,
                           T1=T1
      )
      
      BLM1_fit_pred <- jags(data = LM_Data_pred,inits = NULL,
                            parameters = c("Tcpropagated","Terr"),
                            model = textConnection(BLM1Pred), n.chains = 3, 
                            n.iter = n.iter,  n.burnin = n.iter*burninFrac)
      cbind.data.frame(D47Pred,D47Prederror, Temp=BLM1_fit_pred$BUGSoutput$mean$Tcpropagated,
                       SE=BLM1_fit_pred$BUGSoutput$mean$Terr)
    })
    
    valsPreds <- do.call(rbind, valsPreds)
    
    Preds_NE <- valsPreds %>% 
      group_by(D47Pred, D47Prederror) %>% 
      summarise(Temp = mean(Temp),
                SE = mean(SE), .groups = 'drop')
    
    CompleteModelFit<-list(Preds=Preds, Preds_NE=Preds_NE)
  
  return(CompleteModelFit)
}





