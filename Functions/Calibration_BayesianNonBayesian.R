#' Bayesian regressions to calibrate the clumped isotopes paleotermometer using
#' stan.
#' 
#' @param calibrationData The calibration dataset.
#' @param numSavedSteps Number of MCMC iterations to save.
#' @param priors Either \code{Informative}, \code{Weak}, or 
#'               \code{Uninformative} on the slope and intercept.
#' @param samples Number of samples in the \code{calibrationData} to analyze.


fitClumpedRegressions <<- function(calibrationData, 
                                   numSavedSteps = 3000, 
                                   priors = "Informative",
                                   samples=NULL){
  
  if(! priors %in% c("Informative", "Weak", "Uninformative") ){ 
    stop("Priors must be in `Informative`, `Difusse` or `NonInformative`")
  }
  
  if(is.null(samples)){
    warning("Using the full dataset in the calibration step.")
  }else{
    warning("Sampling ", samples, " observations from the dataset.")
    calibrationData <- calibrationData[sample(1:nrow(calibrationData), samples, replace = TRUE), ]
  }
  
  
  if(priors == "Informative"){
    beta_mu =  0.039
    beta_sd = 0.004
    alpha_mu = 0.231
    alpha_sd = 0.065
  }
  
  if(priors == "Uninformative"){
    beta_mu =  0.01
    beta_sd = 0.01
    alpha_mu = 0.01
    alpha_sd = 0.01
  }
  
  if(priors == "Weak"){
    beta_mu =  0.039
    beta_sd = 0.004 * 2
    alpha_mu = 0.231
    alpha_sd = 0.065 * 2
  }
  
  
  ##Models
  fwMod_Errors = "
  data {
    int<lower=0> N; 
    vector[N] y; 
    vector[N] x_meas; 
    real<lower=0> tau;
    real beta_mu;
    real beta_sd;
    real alpha_mu;
    real alpha_sd;
    real mu_x;
    real sigma_x;
  }
  
  parameters {
    vector[N] x;
    real alpha; 
    real beta; 
    real<lower=0> sigma;
  }
  
  model {
    beta ~ normal(beta_mu, beta_sd);
    alpha ~ normal(alpha_mu, alpha_sd);
    
    x ~ normal(mu_x, sigma_x);
    x_meas ~ normal(x, tau); 
    y ~ normal(alpha + beta * x, sigma);
  
    sigma ~ cauchy(0, 5);
  }
  
  generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | 0, sigma);
    }
  }
"
  
  
  fwMod_NE = "
  data {
    int<lower = 0> N;
    vector[N] x;
    vector[N] y;
    real beta_mu;
    real beta_sd;
    real alpha_mu;
    real alpha_sd;
  }

  parameters {
    real alpha;
    real beta;
    real<lower=0> sigma;
  }
  
  model {
    alpha ~ normal(alpha_mu, alpha_sd);
    beta ~ normal(beta_mu, beta_sd);
    sigma ~ cauchy(0, 5);
    y ~ normal(alpha + beta * x, sigma);
  }
  
  generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | 0, sigma);
    }
  }
  
"
  
  fwMod_mixed = "
  data {
    int<lower=0> N;
    int<lower=0> J;
    vector[N] y;
    vector[N] x;
    int Material[N];
    real beta_mu;
    real beta_sd;
    real alpha_mu;
    real alpha_sd;
  }
  
  parameters {
    real<lower=0> sigma;
    vector[J] alpha;
    vector[J] beta;
  }
  
  model {
    alpha ~ normal(alpha_mu, alpha_sd);
    beta ~ normal(beta_mu, beta_sd);
    y ~ normal(alpha[Material] + beta[Material].*x, sigma);
    sigma ~ cauchy(0, 5);
  }
  
  generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | 0, sigma);
    }
  }
"
  
  ## Non-Informative
  
  fwMod_NE2 = "
  data {
    int<lower = 0> N;
    vector[N] x;
    vector[N] y;
  }

  parameters {
    real alpha;
    real beta;
    real<lower=0> sigma;
  }
  
  model {
    sigma ~ cauchy(0, 5);
    y ~ normal(alpha + beta * x, sigma);
  }
  
  generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | 0, sigma);
    }
  }
"
  
  fwMod_Errors2 = "
  data {
    int<lower=0> N; 
    vector[N] y; 
    vector[N] x_meas; 
  }
  
  parameters {
    vector[N] x;
    real tau; 
    real alpha; 
    real beta; 
    real mu_x;
    real sigma_x;
    real<lower=0> sigma;
  }
  
  model {
    x ~ normal(mu_x, sigma_x);
    x_meas ~ normal(x, tau); 
    y ~ normal(alpha + beta * x, sigma);
  
    sigma ~ cauchy(0, 5);
  }
  
  generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | 0, sigma);
    }
  }
"
  
  
  fwMod_mixed2 = "
  data {
    int<lower=0> N;
    int<lower=0> J;
    vector[N] y;
    vector[N] x;
    int Material[N];
  }
  
  parameters {
    real<lower=0> sigma;
    vector[J] alpha;
    vector[J] beta;
  }
  
  model {
    y ~ normal(alpha[Material] + beta[Material].*x, sigma);
    sigma ~ cauchy(0, 5);
  }
  
  generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | 0, sigma);
    }
  }
"
  
  
  #Data
  
  stan_data_NE <- list(N = nrow(calibrationData), 
                       x = calibrationData$Temperature,
                       y = calibrationData$D47, 
                       beta_mu =  beta_mu,
                       beta_sd = beta_sd,
                       alpha_mu = alpha_mu,
                       alpha_sd = alpha_sd)
  
  stan_data_Err <- list(N = nrow(calibrationData), 
                        x_meas = calibrationData$Temperature,
                        y = calibrationData$D47, 
                        tau = mean(calibrationData$TempError), 
                        beta_mu =  beta_mu,
                        beta_sd = beta_sd,
                        alpha_mu = alpha_mu,
                        alpha_sd = alpha_sd,
                        mu_x = mean(calibrationData$Temperature),
                        sigma_x = sd(calibrationData$Temperature))
  
  stan_data_mixed <- list(N = nrow(calibrationData), 
                          x = calibrationData$Temperature,
                          y = calibrationData$D47, 
                          J = length(unique(calibrationData$Material)),
                          Material = as.numeric(calibrationData$Material), 
                          beta_mu =  beta_mu,
                          beta_sd = beta_sd,
                          alpha_mu = alpha_mu,
                          alpha_sd = alpha_sd)
  
  #Parameters for the run
  nChains = 2
  burnInSteps = 1000
  thinSteps = 1
  nIter = ceiling(burnInSteps + (numSavedSteps * thinSteps)/nChains)
  
  #Fit models
  options(mc.cores = parallel::detectCores())
  BLM1_E <- stan(data = stan_data_Err, model_code = if( priors == "Uninformative") {fwMod_Errors2}else{fwMod_Errors}, 
                 chains = nChains, iter = nIter, warmup = burnInSteps,
                 thin = thinSteps, pars = c('alpha', 'beta', 'sigma', 'log_lik'))
  
  BLM1_NE <- stan(data = stan_data_NE, model_code =  if( priors == "Uninformative") {fwMod_NE2}else{fwMod_NE}, 
                  chains = nChains, iter = nIter, warmup = burnInSteps,
                  thin = thinSteps, pars = c('alpha', 'beta', 'sigma', 'log_lik'))
  
  BLM3 <- stan(data = stan_data_mixed, model_code =  if( priors == "Uninformative") {fwMod_mixed2}else{fwMod_mixed}, 
               chains = 2, iter = nIter, warmup = burnInSteps,
               thin = thinSteps, pars = c('alpha', 'beta', 'sigma', 'log_lik'))
  
  ##
  log_lik_BLM1_E = extract_log_lik(BLM1_E, merge_chains = F)
  r_eff_1_BLM1_E = relative_eff(log_lik_BLM1_E)
  log_lik_BLM1_NE = extract_log_lik(BLM1_NE, merge_chains = F)
  r_eff_BLM1_NE = relative_eff(log_lik_BLM1_NE)
  log_lik_BLM3 = extract_log_lik(BLM3, merge_chains = F)
  r_eff_BLM3 = relative_eff(log_lik_BLM3)
  
  loo_BLM1_E <- loo(log_lik_BLM1_E, r_eff = r_eff_1_BLM1_E)
  loo_BLM1_NE <- loo(log_lik_BLM1_NE, r_eff = r_eff_BLM1_NE)
  loo_BLM3 <- loo(log_lik_BLM3, r_eff = r_eff_BLM3)
  
  looComp <- loo_compare(list('BLM1_E' = loo_BLM1_E, 'BLM1_NE' = loo_BLM1_NE, 'BLM3' = loo_BLM3))

  
  CompleteModelFit <- list("BLM1_fit" = BLM1_E
                           ,"BLM1_fit_NoErrors" = BLM1_NE
                           , "BLM3_fit" = BLM3
  )
  
  attr(CompleteModelFit, "loo") <- looComp 
  #attr(CompleteModelFit, "R2s") <- R2sComplete 
  #attr(CompleteModelFit, "DICs") <- DICs 
  
  return(CompleteModelFit)
}


#' Bootstrap York regression models from a calibration dataset
#' 
#' @param data The calibration dataset
#' @param replicates Number of bootstrap replicates
#' @param samples Number of samples per replicate
#' @param D47error The column in data containing the errors in D47

simulateYork_measured <<- function(data, 
                                   replicates, 
                                   samples = NULL, 
                                   D47error = "D47error"){
  do.call(rbind,lapply(1:replicates, function(x){
    dataSub <- data[sample(seq_along(data[,1]), if(is.null(samples)){nrow(data)}else{samples}, replace = T),]
    dataSub$y_SE <- dataSub[,D47error]
    dataSub$x_SE <- abs(dataSub$TempError)
    Reg <- york(cbind.data.frame(dataSub$Temperature, dataSub$x_SE, dataSub$D47, dataSub$y_SE))
    cbind.data.frame("alpha" = Reg$a[1],"beta" = Reg$b[1])
  }))
}


#' Bootstrap an OLS regression models from a calibration dataset
#' 
#' @param data The calibration dataset
#' @param replicates Number of bootstrap replicates
#' @param samples Number of samples per replicate
#' @param D47error The column in data containing the errors in D47

simulateLM_measured <<- function(data, 
                                 replicates, 
                                 samples = NULL, 
                                 D47error="D47error"){
  
  a <- lapply(1:replicates, function(x){
    dataSub <- data[sample(seq_along(data[,1]), if(is.null(samples)){nrow(data)}else{samples}, replace = T),]
    Reg <- summary(lm(D47 ~ Temperature,  dataSub))
    res <- cbind.data.frame("alpha" = Reg$coefficients[1,1],"beta" = Reg$coefficients[2,1])
    attr(res, "R2") <- Reg$r.squared
    res
  })
  
  R2s <- unlist(lapply(a, function(x) attributes(x)$R2))
  R2s <- data.frame(median = median(R2s), lwr = quantile(R2s, 0.025), upr = quantile(R2s, 0.975))
  a <- do.call(rbind,a)
  attr(a, "R2") <- R2s
  return(a)
}


#' Bootstrap a weighted OLS regression models from a calibration dataset
#' 
#' @param data The calibration dataset
#' @param replicates Number of bootstrap replicates
#' @param samples Number of samples per replicate
#' @param D47error The column in data containing the errors in D47


simulateLM_inverseweights <<- function(data, 
                                       replicates, 
                                       samples = NULL, 
                                       D47error="D47error"){
  a <- lapply(1:replicates, function(x){
    dataSub <- data[sample(seq_along(data[,1]), if(is.null(samples)){nrow(data)}else{samples}, replace = T),]
    Reg0 <- lm(D47 ~ Temperature,  dataSub)
    wt <- 1 / lm(abs(Reg0$residuals) ~ Reg0$fitted.values)$fitted.values^2
    Reg <- summary(lm(D47 ~ Temperature,  dataSub, weights = wt))
    res <- cbind.data.frame("alpha" = Reg$coefficients[1,1], "beta" = Reg$coefficients[2,1])
    attr(res, "R2") <- Reg$r.squared
    res
  })
  
  R2s <- unlist(lapply(a, function(x) attributes(x)$R2))
  R2s <- data.frame(median = median(R2s), lwr = quantile(R2s, 0.025), upr = quantile(R2s, 0.975))
  a <- do.call(rbind,a)
  attr(a, "R2") <- R2s
  return(a)
}


#' Bootstrap Deming regression models from a calibration dataset
#' 
#' @param data The calibration dataset
#' @param replicates Number of bootstrap replicates
#' @param samples Number of samples per replicate
#' @param D47error The column in data containing the errors in D47


simulateDeming <<- function(data, 
                            replicates, 
                            samples = NULL, 
                            D47error="D47error"){
  
  do.call(rbind,lapply(1:replicates, function(x){
    dataSub <- data[sample(seq_along(data[,1]), if(is.null(samples)){nrow(data)}else{samples}, replace = T),]
    dataSub$y_SE <- abs(dataSub[,D47error])
    dataSub$x_SE <- abs(dataSub$TempError)
    Reg <- deming(D47 ~ Temperature, dataSub, xstd= x_SE, ystd= y_SE)
    cbind.data.frame("alpha" = Reg$coefficients[1],"beta" = Reg$coefficients[2])
  }))
  
}


#' This function is used to generate CI estimates at given intervals. It is currently
#' used for plotting in BayClump.
#' 
#' @param data A data.frame with two columns labeled beta and alpha. 
#' This should be the result of bootstrapping a given calibration set.
#' @param from the lower limit in x
#' @param to the upper limit in x
#' @param length.out the number of breaks


RegressionSingleCI <-function(data, from, to, length.out=100){
  
  sampleDataReplicates <- as.data.frame(data)
  
  Theta_mat <- sampleDataReplicates
  
  # Points where to evaluate the model
  x_eval <- seq(from, to, length.out = length.out)
  
  # Matrix with the predictions
  fun <- function(x, theta) as.numeric(theta["alpha"]) + (x * as.numeric(theta["beta"]))
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

#' Generate a dataset reflecting the priors used to run the analyses
#' 
#' @param prior Informative or not
#' @param n number of observations to simulate

generatePriorDistCalibration <- function(prior, n=1000){
  if (prior == "Informative") {
    params <- cbind.data.frame(parameter = c("alpha", "beta"),
                               mean = c(0.231,0.039), 
                               sd = c(0.065,0.004))
    params
  } else {
    params <- cbind.data.frame(parameter=c("alpha", "beta"),
                               mean=c(0,0.01), 
                               sd=c(0,0.01))
    params
  }
  
  data <- cbind.data.frame(alpha = rnorm(n, params[1,2], params[1,3]), 
                           beta = rnorm(n, params[2,2], params[2,3]))
  attr(data, "priors") <-  prior
  attr(data, "params") <-  params
  data
}



