#' This function generate temperature predictions (in 10^6/T2) based on a 
#' calibration dataset and target D47. 
#' 
#' @param calData The calibration dataset
#' @param calModel The stan model to be analyzed
#' @param recData The reconstruction dataset
#' @param iter Number of replicates to retain
#' @param priors Whether priors should be \code{Noninformative} or not
#' @param prior_mu Prior on mean temperature (scale = 10^6/T2, T in K)
#' @param prior_sig Prior on sd temperature (scale = 10^6/T2, T in K)

BayesianPredictions <- function(calModel,
                                calData,
                                recData,
                                iter = 1000,
                                priors = "Uninformative",
                                prior_mu = 11,
                                prior_sig = 5){
  
  vects.params <- extract(calModel)
  
  ununfpredMod = "
data {
  int<lower=0> n;
  vector[n] y;

  int<lower=0> posts;
  vector[posts] alpha;              
  vector[posts] beta;
  vector[posts] sigma;
}

parameters {
  matrix[n, posts] x_new;
}

model {
  vector[posts] y_new_hat; 
  for(i in 1:n){
    y_new_hat = alpha + beta .* x_new[i,]';
    y[i] ~ normal(y_new_hat, sigma);
}
}
"

predMod = "
data {
  int<lower=0> n;
  vector[n] y;

  int<lower=0> posts;
  vector[posts] alpha;              
  vector[posts] beta;
  vector[posts] sigma;  
  real prior_mu;
  real prior_sig;
}

parameters {
  matrix[n, posts] x_new;
}

model {
  vector[posts] y_new_hat; 
  for(i in 1:n){
    x_new[i,] ~ normal(prior_mu, prior_sig);
    y_new_hat = alpha + beta .* x_new[i,]';
    y[i] ~ normal(y_new_hat, sigma);
}
}
"

stan_date <- list(n = nrow(recData), 
                  y = recData$D47, 
                  posts = length(vects.params$beta), 
                  alpha = vects.params$alpha,
                  beta = vects.params$beta,
                  sigma = vects.params$sigma,
                  prior_mu = prior_mu,
                  prior_sig = prior_sig)

options(mc.cores = parallel::detectCores())
data.rstan <- stan(data = stan_date, model_code = if(priors == "Uninformative"){ununfpredMod}else{predMod}, 
                   chains = 2, iter = iter, warmup = floor(iter/2),
                   , control = list(adapt_delta = 0.90, max_treedepth = 10)
)

params2 <- extract(data.rstan)
Xouts2 <- params2$x_new
Xdims2 <- dim(Xouts2)
xis <- list()
recs <- lapply(1:Xdims2[2], function(x){
  cbind.data.frame(mean(Xouts2[,x,]), quantile(Xouts2[,x,], c(0.025)),
                   quantile(Xouts2[,x,], c(0.975)))
})

recs <- do.call(rbind, recs)

cbind.data.frame(Sample = recData$Sample, 
                 D47 = recData$D47, 
                 D47error = recData$D47error, 
                 meanTemp = sqrt(10^6/recs[,1]) - 273.15, 
                 Temp_L = sqrt(10^6/recs[,3]) - 273.15, 
                 Temp_H = sqrt(10^6/recs[,2]) - 273.15)


}


#' Generate a dataset reflecting the priors used to run the analyses
#' 
#' @param prior Informative or not
#' @param n number of observations to simulate


generatePriorReconstructions <- function(prior, n=1000){
  
  if(prior == "Informative"){
    params <- cbind.data.frame(parameter=c("alpha", "beta"),
                               mean=c(0.231,0.039), 
                               sd=c(0.065,0.004))
    params
  } else {
    params <- cbind.data.frame(parameter=c("alpha", "beta"),
                               mean=c(0,0.01), 
                               sd=c(0,0.01))
    params
  }
  x <- rnorm(n, 11,0.0001)
  xC <- sqrt(10^6/x)-273.15
  
  data <- cbind.data.frame(alpha=rnorm(n, params[1,2], params[1,3]), 
                           beta=rnorm(n, params[2,2], params[2,3]),
                           x,
                           xC)
  attr(data, "priors") <-  prior
  attr(data, "params") <-  params
  data
}


