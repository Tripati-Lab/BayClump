regressionParameters<-function(CompleteModelFit){
  n=nrow(CompleteModelFit$BLM1_fit$BUGSoutput$sims.matrix)
  if(length(CompleteModelFit) == 7){
    
    summ_parameters<-as.data.frame(matrix(NA,5,10))
    colnames(summ_parameters)<-c('Model','Materials','MaterialLevelParameter',"Intercept", "Intercept_95%CI_lower","Intercept_95%CI_upper" ,"Slope", "Slope_95%CI_lower","Slope_95%CI_upper","ESS")
    summ_parameters[1,]<-cbind.data.frame('York','All',NA,CompleteModelFit$Y$a[1], CompleteModelFit$Y$a[1]-1.96*CompleteModelFit$Y$a[2], CompleteModelFit$Y$a[1]+1.96*CompleteModelFit$Y$a[2],CompleteModelFit$Y$b[1], CompleteModelFit$Y$b[1] - CompleteModelFit$Y$b[2], CompleteModelFit$Y$b[1]+1.96*CompleteModelFit$Y$b[2] ,NA)
    summ_parameters[2,]<-cbind.data.frame('Linear model','All',NA,CompleteModelFit$M0$coefficients[1], CompleteModelFit$M0$coefficients[1]-1.96*summary(CompleteModelFit$M0)$coefficients[1,2],CompleteModelFit$M0$coefficients[1]+1.96*summary(CompleteModelFit$M0)$coefficients[1,2], CompleteModelFit$M0$coefficients[2], CompleteModelFit$M0$coefficients[2]-1.96*summary(CompleteModelFit$M0)$coefficients[2,2],CompleteModelFit$M0$coefficients[2]+1.96*summary(CompleteModelFit$M0)$coefficients[2,2],NA)
    summ_parameters[3,]<-cbind.data.frame('Main effects ANCOVA','2',CompleteModelFit$M1$coefficients[3],CompleteModelFit$M1$coefficients[1], CompleteModelFit$M1$coefficients[1]-1.96*summary(CompleteModelFit$M1)$coefficients[1,2], CompleteModelFit$M1$coefficients[1]+1.96*summary(CompleteModelFit$M1)$coefficients[1,2],CompleteModelFit$M1$coefficients[2], CompleteModelFit$M1$coefficients[2]-1.96*summary(CompleteModelFit$M1)$coefficients[2,2],CompleteModelFit$M1$coefficients[2]+1.96*summary(CompleteModelFit$M1)$coefficients[2,2],NA)
    summ_parameters[4,]<-cbind.data.frame('Interaction effects ANCOVA','2',CompleteModelFit$M2$coefficients[3],CompleteModelFit$M2$coefficients[1], CompleteModelFit$M2$coefficients[1]-1.96*summary(CompleteModelFit$M2)$coefficients[1,2], CompleteModelFit$M2$coefficients[1]+1.96*summary(CompleteModelFit$M2)$coefficients[1,2],CompleteModelFit$M2$coefficients[2], CompleteModelFit$M2$coefficients[2]-1.96*summary(CompleteModelFit$M2)$coefficients[1,2],CompleteModelFit$M2$coefficients[2]+1.96*summary(CompleteModelFit$M2)$coefficients[1,2],NA)
    summ_parameters[5,]<-cbind.data.frame('Bayesian Linear Model','All',NA,CompleteModelFit$BLM1_fit$BUGSoutput$summary[1,1], CompleteModelFit$BLM1_fit$BUGSoutput$summary[1,3],CompleteModelFit$BLM1_fit$BUGSoutput$summary[1,7], CompleteModelFit$BLM1_fit$BUGSoutput$summary[2,1], CompleteModelFit$BLM1_fit$BUGSoutput$summary[2,3],CompleteModelFit$BLM1_fit$BUGSoutput$summary[2,7],paste0('alpha=',CompleteModelFit$BLM1_fit$BUGSoutput$summary[1,9], ', beta=',CompleteModelFit$BLM1_fit$BUGSoutput$summary[2,9]))
    
    ##For the main effects ANCOVA
    BLMOu1<-head(CompleteModelFit$BLM2_fit$BUGSoutput$summary,-2)
    Materials<-as.numeric(gsub("[^\\d]+", "", row.names(BLMOu1), perl=TRUE))
    BME_s<-do.call(rbind.data.frame,lapply(unique(na.omit(Materials)), function(z){
      cbind.data.frame('Bayesian Main effects ANCOVA','All',z,BLMOu1[1,1], BLMOu1[1,3],BLMOu1[1,7], BLMOu1[1+z,1], BLMOu1[1+z,3],BLMOu1[1+z,7],paste0('alpha=',BLMOu1[1,9], ', beta=',BLMOu1[z+1,9]))
    }))
    names(BME_s)<-colnames(summ_parameters)
    summ_parameters<-rbind.data.frame(summ_parameters,BME_s)
    
    ##For the interaction effects ANCOVA
    BLMOu2<-head(CompleteModelFit$BLM3_fit$BUGSoutput$summary,-2)
    Materials<-as.numeric(gsub("[^\\d]+", "", row.names(BLMOu1), perl=TRUE))
    BME_s<-do.call(rbind.data.frame,lapply(unique(na.omit(Materials)), function(z){
      cbind.data.frame('Bayesian Interaction effects ANCOVA','All',z,BLMOu2[z,1], BLMOu2[z,3],BLMOu2[z,7], BLMOu2[2+z,1], BLMOu2[2+z,3], BLMOu2[2+z,7],paste0('alpha=',BLMOu2[z,9], ', beta=',BLMOu2[z+2,9]))
    }))
    names(BME_s)<-colnames(summ_parameters)
    
    summ_parameters<-rbind.data.frame(summ_parameters,BME_s)
  }else{
    summ_parameters<-as.data.frame(matrix(NA,3,10))
    colnames(summ_parameters)<-c('Model','Materials','MaterialLevelParameter',"Intercept", "Intercept_95%CI_lower","Intercept_95%CI_upper" ,"Slope", "Slope_95%CI_lower","Slope_95%CI_upper","ESS")
    summ_parameters[1,]<-cbind.data.frame('York','All',NA,CompleteModelFit$Y$a[1], CompleteModelFit$Y$a[1]-1.96*CompleteModelFit$Y$a[2], CompleteModelFit$Y$a[1]+1.96*CompleteModelFit$Y$a[2],CompleteModelFit$Y$b[1], CompleteModelFit$Y$b[1] - CompleteModelFit$Y$b[2], CompleteModelFit$Y$b[1]+1.96*CompleteModelFit$Y$b[2] ,NA)
    summ_parameters[2,]<-cbind.data.frame('Linear model','All',NA,CompleteModelFit$M0$coefficients[1], CompleteModelFit$M0$coefficients[1]-1.96*summary(CompleteModelFit$M0)$coefficients[1,2],CompleteModelFit$M0$coefficients[1]+1.96*summary(CompleteModelFit$M0)$coefficients[1,2], CompleteModelFit$M0$coefficients[2], CompleteModelFit$M0$coefficients[2]-1.96*summary(CompleteModelFit$M0)$coefficients[2,2],CompleteModelFit$M0$coefficients[2]+1.96*summary(CompleteModelFit$M0)$coefficients[2,2],NA)
    summ_parameters[3,]<-cbind.data.frame('Bayesian Linear Model','All',NA,CompleteModelFit$BLM1_fit$BUGSoutput$summary[1,1], CompleteModelFit$BLM1_fit$BUGSoutput$summary[1,3],CompleteModelFit$BLM1_fit$BUGSoutput$summary[1,7], CompleteModelFit$BLM1_fit$BUGSoutput$summary[2,1], CompleteModelFit$BLM1_fit$BUGSoutput$summary[2,3],CompleteModelFit$BLM1_fit$BUGSoutput$summary[2,7],paste0('alpha=',CompleteModelFit$BLM1_fit$BUGSoutput$summary[1,9], ', beta=',CompleteModelFit$BLM1_fit$BUGSoutput$summary[2,9]))
    summ_parameters<-summ_parameters[,-c(2,3)]
     }
  attr(summ_parameters, 'data') <- attr(CompleteModelFit, "data") 
  return(summ_parameters)
  
}


se <- function(x,n) sd(x)/sqrt(n)
