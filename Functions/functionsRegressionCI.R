RegressionSingleCI<-function(data, from, to, length.out=100){
  
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
  
  return(list(uncertaintyModels))
  
}
