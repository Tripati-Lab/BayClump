RegressionMultipleCI<-function(data, trueAlpha=0.268, trueBeta=0.0369, from=5, to=15){
  
  sampleDataReplicates<- as.data.frame(data)
  colnames(sampleDataReplicates)[1] <- 'model'
  
  uncertaintyModels<-do.call(rbind,lapply(unique(sampleDataReplicates[,1]), function(x){
    
    Theta_mat <- sampleDataReplicates[sampleDataReplicates[,1] == x,]
    
    # Points where to evaluate the model
    x_eval <- seq(from, to, length.out = 100)
    
    # Matrix with the predictions
    fun <- function(x, theta) as.numeric(theta["intercept"]) + (x * as.numeric(theta["slope"]))
    Pred_mat <- apply(Theta_mat, 1, function(theta) fun(x_eval, theta))
    
    
    # Pack the estimates for plotting
    Estims_plot <- cbind(
      model=x,
      x = x_eval, 
      as.data.frame(t(apply(Pred_mat, 1, function(y_est) c(
        median_est = median(y_est), 
        ci_lower_est = quantile(y_est, probs = 0.025, names = FALSE), 
        ci_upper_est = quantile(y_est, probs = 0.975, names = FALSE)
      ))))
    )
    
    
  }))
  uncertaintyModels<-cbind.data.frame(uncertaintyModels, trueD47=(trueAlpha + (trueBeta * x_eval)))
  
  a<-do.call(rbind,lapply(unique(uncertaintyModels$model), function(x){
    sbd<-uncertaintyModels[uncertaintyModels$model ==x, ]
    cbind.data.frame(model=x,maxUpper=median(sbd$ci_upper_est), minLower= median(sbd$ci_lower_est),
                     median=median(sbd$median_est),
                     distance95=median(sbd$ci_upper_est)-median(sbd$ci_lower_est)
    )
  }))
  
  summarystats<-do.call(rbind,lapply(unique(sampleDataReplicates$model), function(x){
    sub<-sampleDataReplicates[sampleDataReplicates$model == x,]
    cbind.data.frame(model=sub$model[1], intercept=median(sub$intercept),
                     lowInt=quantile(sub$intercept, 0.025),
                     upInt=quantile(sub$intercept, 0.975),
                     slope=median(sub$slope),
                     lowSl=quantile(sub$slope, 0.025),
                     upSl=quantile(sub$slope, 0.975)
    )
  }))
  
  
  p<-ggplot(sampleDataReplicates) + 
    #geom_point(data = data,aes(x = x_measured, y = y_measured), fill='black',color='white', shape = 21) + 
    #geom_errorbar(data = data, aes(x = x_measured, ymin = y_measured-y_SD,ymax = y_measured+y_SD), col='black' ) + 
    #geom_errorbarh(data = data, aes(y = y_measured, xmin = x_measured-x_SD,xmax = x_measured+x_SD),col='black')+
    geom_abline(data=sampleDataReplicates,aes(intercept = intercept, slope = slope), alpha=0.5, col='#4e888a')+
    geom_ribbon(data=uncertaintyModels,aes(x = x, y = median_est, ymin = ci_lower_est, ymax = ci_upper_est), fill = "#ffd166")+
    geom_abline(intercept = 0.268, slope = 0.0369, col='#ef476f', size=0.8)+
    geom_line(data=uncertaintyModels,aes(x = x, y = median_est), col='black', size=0.8, lty='dashed')+
    coord_cartesian(ylim = c(0.5, 0.9),
                    xlim=c(10,14)) +
    theme(axis.line =element_blank(),
          axis.ticks = element_line(colour = "black", size = .2),
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
    facet_wrap(~model, ncol=3,scales = "free")
  
  return(list(plot=p, dataUncertainties=uncertaintyModels))
  
}
RegressionSingleCI<-function(data, trueAlpha=0.268, trueBeta=0.0369, from=5, to=15){
  
  sampleDataReplicates<- as.data.frame(data)
  #colnames(sampleDataReplicates)[1] <- 'model'
  
  Theta_mat <- sampleDataReplicates
  
  # Points where to evaluate the model
  x_eval <- seq(from, to, length.out = 100)
  
  # Matrix with the predictions
  fun <- function(x, theta) as.numeric(theta["intercept"]) + (x * as.numeric(theta["slope"]))
  Pred_mat <- apply(Theta_mat, 1, function(theta) fun(x_eval, theta))
  
  
  # Pack the estimates for plotting
  Estims_plot <- cbind(
    model=x,
    x = x_eval, 
    as.data.frame(t(apply(Pred_mat, 1, function(y_est) c(
      median_est = median(y_est), 
      ci_lower_est = quantile(y_est, probs = 0.025, names = FALSE), 
      ci_upper_est = quantile(y_est, probs = 0.975, names = FALSE)
    ))))
  )
  
  
  uncertaintyModels<-cbind.data.frame(Estims_plot, trueD47=(trueAlpha + (trueBeta * x_eval)))
  
  a<-do.call(rbind,lapply(unique(uncertaintyModels$model), function(x){
    sbd<-uncertaintyModels[uncertaintyModels$model ==x, ]
    cbind.data.frame(model=x,maxUpper=median(sbd$ci_upper_est), minLower= median(sbd$ci_lower_est),
                     median=median(sbd$median_est),
                     distance95=median(sbd$ci_upper_est)-median(sbd$ci_lower_est)
    )
  }))
  
  summarystats<-do.call(rbind,lapply(unique(sampleDataReplicates$model), function(x){
    sub<-sampleDataReplicates[sampleDataReplicates$model == x,]
    cbind.data.frame(model=sub$model[1], intercept=median(sub$intercept),
                     lowInt=quantile(sub$intercept, 0.025),
                     upInt=quantile(sub$intercept, 0.975),
                     slope=median(sub$slope),
                     lowSl=quantile(sub$slope, 0.025),
                     upSl=quantile(sub$slope, 0.975)
    )
  }))
  
  
  p<-ggplot(sampleDataReplicates) + 
    #geom_point(data = data,aes(x = x_measured, y = y_measured), fill='black',color='white', shape = 21) + 
    #geom_errorbar(data = data, aes(x = x_measured, ymin = y_measured-y_SD,ymax = y_measured+y_SD), col='black' ) + 
    #geom_errorbarh(data = data, aes(y = y_measured, xmin = x_measured-x_SD,xmax = x_measured+x_SD),col='black')+
    geom_abline(data=sampleDataReplicates,aes(intercept = intercept, slope = slope), alpha=0.5, col='#4e888a')+
    geom_ribbon(data=uncertaintyModels,aes(x = x, y = median_est, ymin = ci_lower_est, ymax = ci_upper_est), fill = "#ffd166")+
    geom_abline(intercept = 0.268, slope = 0.0369, col='#ef476f', size=0.8)+
    geom_line(data=uncertaintyModels,aes(x = x, y = median_est), col='black', size=0.8, lty='dashed')+
    coord_cartesian(ylim = c(0.5, 0.9),
                    xlim=c(10,14)) +
    theme(axis.line =element_blank(),
          axis.ticks = element_line(colour = "black", size = .2),
          axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  
  return(list(plot=p, dataUncertainties=uncertaintyModels))
  
}
