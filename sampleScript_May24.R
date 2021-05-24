##Packages
library(pbapply)
library(deming)
library(pbmcapply)
library(IsoplotR)
library(deming)
library(R2jags)
library(data.table)
library(ggplot2)

#Read the relevant scripts in the functions folder
sapply(list.files('Functions', full.names = T), source)

##The dataset with the default structure
data<-read.csv('Data/SampleData.csv')


#Run each model independently using 5 replicates and 
#100 observations per replicate. The number of generations for 
#the Bayesian models is set to 1,000 but we generally need 20.000.
#isMixed won't fit the mixed model (no Material column needed).

simulateYork_measured(data=data, replicates=5, samples=100)
simulateLM_measured(data=data, replicates=5, samples=100)
simulateDeming(data=data, replicates=5, samples=100)
simulateLM_inverseweights(data=data, replicates=5, samples=100)
simulateBLM_measuredMaterial(data=data, replicates=5, samples=100, generations=1000, isMixed=F)
simulateBLM_measuredMaterial(data=data, replicates=5, samples=100, generations=1000, isMixed=T)

##Now, run the seven regression models using 5 replicates per model and sampling 100
##observations per replicate. 
simulateAll(data, replicates=5, samples=100, generations=1000)

##If the user wants to analyze also the mixed model, use isMixed=T
simulateAll(data, replicates=5, samples=100, generations=1000, isMixed=T)



##Regression confidence intervals

##First, if only the confidence interval for a single model is requested
datReplicatesSingle<-simulateYork_measured(data=data, replicates=5, samples=100)
outCI<-RegressionSingleCI(data=datReplicatesSingle, from=5, to=15, length.out=100)
outCI$dataUncertainties





