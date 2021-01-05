
#T2--------

predictYork<-function(yorkObject,x, x_error, n_measurements=2, nrep=100){
  
  la<-yorkObject$a[1]-qt(p=.05/2, df=yorkObject$df, lower.tail=F)*yorkObject$a[2]
  ua<-yorkObject$a[1]+qt(p=.05/2, df=yorkObject$df, lower.tail=F)*yorkObject$a[2]
  
  lb<-yorkObject$b[1]-qt(p=.05/2, df=yorkObject$df, lower.tail=F)*yorkObject$b[2]
  ub<-yorkObject$b[1]+qt(p=.05/2, df=yorkObject$df, lower.tail=F)*yorkObject$b[2]
  
  sda<- yorkObject$a[2] * sqrt(yorkObject$df+2)
  sdb<- yorkObject$b[2] * sqrt(yorkObject$df+2)
  
 
  resdat<-as.data.frame(t(sapply(seq_along(x), function(m){
  dist<-sapply(1:nrep, function(w){
     pa<- rnorm(1, mean=yorkObject$a[1], sd=sda)
     pb<-  rnorm(1, mean=yorkObject$b[1], sd=sdb)
     pval<- rnorm(1, mean=x, sd= x_error[m]* sqrt(n_measurements[m]))
     (pval-pa)/pb
    })
   cbind.data.frame(mean=mean(dist), SE=std(dist) , confidence_interval(dist, .95))
  })))
  
  colnames(resdat)<-c('York_mean', 'York_SE', 'York_lower95','York_upper95')
  resdat
}


predictLm<-function(lmObject, x, x_error, n_measurements=2, nrep=100){
  
  a<-do.call(rbind.data.frame,lapply(seq_along(x), function(y){
    dist<-as.data.frame(t(sapply(1:nrep, function(q){
    pval<- rnorm(1, mean=x[y], sd= x_error[y]* sqrt(n_measurements[y]))
   as.data.frame(t(as.data.frame(unlist(inverse.predict(lmObject, pval )))))[,c(1,2,4,5)]
    })))
    dist<-unlist(dist$Prediction)
    cbind.data.frame(mean=mean(dist), SE=std(dist) , confidence_interval(dist, .95))
  }))
 row.names(a)<-NULL
 colnames(a)<-c("LM_mean",'LM_SE', 'LM_lower95', 'LM_upper95')
 a
}

predictLmCovariate<-function(calibrationData, lmObject, x, x_error, n_measurements=2, materials, nrep=100){
  
  a<-do.call(rbind.data.frame,lapply(seq_along(x), function(y){
      dist<-lapply(1:nrep, function(q){
        tryCatch({
        pval<- rnorm(1, mean=x[y], sd= x_error[y]* sqrt(n_measurements[y]))
    invest(lmObject, y0 = pval, x0.name = "T2", newdata = data.frame(Material = materials[y]))
      }, error=function(e){NA})
      })
      dist<-do.call(rbind.data.frame,dist)[,1]
      cbind.data.frame(mean=mean(dist, na.rm=T), SE=std(dist) , confidence_interval(dist, .95))
  }))
  
  row.names(a)<-NULL
  colnames(a)<-c("ANCOVA_mean","ANCOVA_SE", 'ANCOVA_lower95', 'ANCOVA_upper95')
  a
}


predictLmJags<-function(jagsModel, x, x_error, n_measurements=2, nrep=100 ){
  
   mcmc = as.data.frame(jagsModel$BUGSoutput$sims.matrix)
   do.call(rbind.data.frame,lapply(seq_along(x), function(y){
     dist<-as.data.frame(t(sapply(1:nrep, function(q){
       pval<- rnorm(1, mean=x[y], sd= x_error[y]* sqrt(n_measurements[y]))
       
       posterior_target<-as.mcmc((pval-mcmc$alpha)/mcmc$beta)
       a<-cbind.data.frame(mean=mean(posterior_target),
                           se(posterior_target, n=length(posterior_target)) , 
                           HPDinterval(posterior_target))
       row.names(a)<-NULL
       colnames(a)<-c("BLM_mean",'BLM_SE', 'BLM_lower95', 'BLM_upper95')
       a
     })))
     dist<-unlist(dist[,1])
     dist<-cbind.data.frame(mean=mean(dist, na.rm=T), SE=std(dist) , confidence_interval(dist, .95))
     colnames(dist)<-c("BLM_mean",'BLM_SE', 'BLM_lower95', 'BLM_upper95')
     dist
   }))
   
}

predictANCOVA1Jags<-function(jagsModel, x, x_error, n_measurements=2, materials, nrep=100){
  mcmc = as.data.frame(jagsModel$BUGSoutput$sims.matrix)
  NumbCols<-as.numeric(gsub("[^\\d]+", "", colnames(mcmc), perl=TRUE))
  
  do.call(rbind.data.frame,lapply(seq_along(x), function(y){
    dist<-as.data.frame(t(sapply(1:nrep, function(q){
      pval<- rnorm(1, mean=x[y], sd= x_error[y]* sqrt(n_measurements[y]))
    posterior_target<-as.mcmc((pval-mcmc$alpha)/mcmc[which(NumbCols ==materials[y])  ])
    a<-cbind.data.frame(mean=mean(posterior_target),
                        se(posterior_target, n=length(posterior_target)), 
                        HPDinterval(posterior_target))
    row.names(a)<-NULL
    colnames(a)<-c("BmainANCOVA_mean", 'BmainANCOVA_SE','BmainANCOVA_lower95', 'BmainANCOVA_upper95')
    a
    })))
    dist<-unlist(dist[,1])
    dist<-cbind.data.frame(mean=mean(dist, na.rm=T), SE=std(dist) , confidence_interval(dist, .95))
    colnames(dist)<-c("BmainANCOVA_mean", 'BmainANCOVA_SE','BmainANCOVA_lower95', 'BmainANCOVA_upper95')
    dist
  }))
}

predictANCOVA2Jags<-function(jagsModel, x, x_error, n_measurements=2, materials, nrep=100){
  mcmc = as.data.frame(jagsModel$BUGSoutput$sims.matrix)
  
  inte<-mcmc[,grep('alpha', colnames(mcmc))]
  NumbCols_int<-as.numeric(gsub("[^\\d]+", "", colnames(inte), perl=TRUE))
  
  slopes<-mcmc[,grep('beta', colnames(mcmc))]
  NumbCols_slopes<-as.numeric(gsub("[^\\d]+", "", colnames(slopes), perl=TRUE))
  
  do.call(rbind.data.frame,lapply(seq_along(x), function(y){
    dist<-as.data.frame(t(sapply(1:nrep, function(q){
      pval<- rnorm(1, mean=x[y], sd= x_error[y]* sqrt(n_measurements[y]))
    posterior_target<-as.mcmc((pval-inte[,which(NumbCols_int == materials[y])])/slopes[,which(NumbCols_slopes == materials[y]) ])
    a<-cbind.data.frame(mean=mean(posterior_target),
                        se(posterior_target, n=length(posterior_target)), 
                        HPDinterval(posterior_target))
    row.names(a)<-NULL
    colnames(a)<-c("BInteANCOVA_mean",'BInteANCOVA_SE', 'BInteANCOVA_lower95', 'BInteANCOVA_upper95')
    a
    })))
    dist<-unlist(dist[,1])
    dist<-cbind.data.frame(mean=mean(dist, na.rm=T), SE=std(dist) , confidence_interval(dist, .95))
    colnames(dist)<-c("BmainANCOVA_mean", 'BmainANCOVA_SE','BmainANCOVA_lower95', 'BmainANCOVA_upper95')
    dist
  }))
}


##D47 (these doesn't account for uncertainty in input values)------


CalibratePredictYork<-function(yorkObject,x){
  la<-yorkObject$a[1]-qt(p=.05/2, df=yorkObject$df, lower.tail=F)*yorkObject$a[2]
  ua<-yorkObject$a[1]+qt(p=.05/2, df=yorkObject$df, lower.tail=F)*yorkObject$a[2]
  
  lb<-yorkObject$b[1]-qt(p=.05/2, df=yorkObject$df, lower.tail=F)*yorkObject$b[2]
  ub<-yorkObject$b[1]+qt(p=.05/2, df=yorkObject$df, lower.tail=F)*yorkObject$b[2]
  
  
  cbind.data.frame('York_mean'=(x$T2 *yorkObject$b[1])+yorkObject$a[1],
                   'York_SE'=(((x$T2-(la))/(lb) - (x$T2-(ua))/(ub))/2)/sqrt(length(x$T2)),
                   'York_lower95'=(x$T2*ub+(ua)),
                   'York_upper95'=(x$T2*lb+(la))
  )
}

CalibratePredictLm<-function(lmObject, x){
  prvs<-predict(lmObject,  newdata = data.frame(T2=x), se.fit = T)
  lb<-prvs$fit-qt(p=.05/2, df=prvs$df, lower.tail=F)*prvs$se.fit
  up<-prvs$fit+qt(p=.05/2, df=prvs$df, lower.tail=F)*prvs$se.fit
  data.frame("LM_mean"=prvs$fit, "LM_SE"=prvs$se.fit ,'LM_lower95'=lb, 'LM_upper95'=up)
  
}

CalibratePredictLmCovariate<-function(calibrationData, lmObject, x, materials){
      prvs<-predict(lmObject,  newdata = data.frame(T2=x,Material = materials), se.fit = T)
      lb<-prvs$fit-qt(p=.05/2, df=prvs$df, lower.tail=F)*prvs$se.fit
      up<-prvs$fit+qt(p=.05/2, df=prvs$df, lower.tail=F)*prvs$se.fit
      data.frame("ANCOVA_mean"=prvs$fit, 'ANCOVA_SE'=prvs$se.fit , 'ANCOVA_lower95'=lb, 'ANCOVA_upper95'=up)
}


CalibratePredictLmJags<-function(jagsModel, x){
  mcmc = as.data.frame(jagsModel$BUGSoutput$sims.matrix)
  do.call(rbind.data.frame,lapply(seq_along(x), function(y){
    posterior_target<-as.mcmc((x[y]*mcmc$beta+mcmc$alpha))
    a<-cbind.data.frame(mean=mean(posterior_target),
                        se(posterior_target, n=length(posterior_target)), 
                        HPDinterval(posterior_target))
    row.names(a)<-NULL
    colnames(a)<-c("BLM_mean", 'BLM_SE' ,'BLM_lower95', 'BLM_upper95')
    a
  }))
}

CalibratePredictANCOVA1Jags<-function(jagsModel, x, materials){
  mcmc = as.data.frame(jagsModel$BUGSoutput$sims.matrix)
  NumbCols<-as.numeric(gsub("[^\\d]+", "", colnames(mcmc), perl=TRUE))
  
  do.call(rbind.data.frame,lapply(seq_along(x), function(y){
    posterior_target<-as.mcmc((x[y]*mcmc[which(NumbCols ==materials[y])  ]+mcmc$alpha))
    a<-cbind.data.frame(mean=mean(posterior_target),
                        se(posterior_target, n=length(posterior_target)), 
                        HPDinterval(posterior_target))
    row.names(a)<-NULL
    colnames(a)<-c("BmainANCOVA_mean", 'BmainANCOVA_SE' ,'BmainANCOVA_lower95', 'BmainANCOVA_upper95')
    a
  }))
}

CalibratePredictANCOVA2Jags<-function(jagsModel, x, materials){
  mcmc = as.data.frame(jagsModel$BUGSoutput$sims.matrix)
  
  inte<-mcmc[,grep('alpha', colnames(mcmc))]
  NumbCols_int<-as.numeric(gsub("[^\\d]+", "", colnames(inte), perl=TRUE))
  
  slopes<-mcmc[,grep('beta', colnames(mcmc))]
  NumbCols_slopes<-as.numeric(gsub("[^\\d]+", "", colnames(slopes), perl=TRUE))
  
  do.call(rbind.data.frame,lapply(seq_along(x), function(y){
    posterior_target<-as.mcmc((x[y]*slopes[,which(NumbCols_slopes == materials[y]) ]+inte[,which(NumbCols_int == materials[y])]))
    a<-cbind.data.frame(mean=mean(posterior_target),
                        se(posterior_target, n=length(posterior_target)), 
                        HPDinterval(posterior_target))
    row.names(a)<-NULL
    colnames(a)<-c("BInteANCOVA_mean", 'BInteANCOVA_SE' , 'BInteANCOVA_lower95', 'BInteANCOVA_upper95')
    a
  }))
}



