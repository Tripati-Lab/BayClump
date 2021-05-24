RegressionSingleCI<-function(data, from=5, to=15, length.out=100){
  
  sampleDataReplicates<- as.data.frame(data)
  #colnames(sampleDataReplicates)[1] <- 'model'
  
  Theta_mat <- sampleDataReplicates
  
  # Points where to evaluate the model
  x_eval <- seq(from, to, length.out = length.out)
  
  # Matrix with the predictions
  fun <- function(x, theta) as.numeric(theta["intercept"]) + (x * as.numeric(theta["slope"]))
  Pred_mat <- apply(Theta_mat, 1, function(theta) fun(x_eval, theta))
  
  
  # Pack the estimates for plotting
  uncertaintyModels <- cbind(
    x = x_eval, 
    as.data.frame(t(apply(Pred_mat, 1, function(y_est) c(
      median_est = median(y_est), 
      ci_lower_est = quantile(y_est, probs = 0.025, names = FALSE), 
      ci_upper_est = quantile(y_est, probs = 0.975, names = FALSE)
    ))))
  )
  
  
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
