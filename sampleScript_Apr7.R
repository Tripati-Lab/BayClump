##Packages
library(pbapply)
library(deming)
library(pbmcapply)
library(IsoplotR)
library(deming)
library(R2jags)
library(data.table)

#Read the relevant scripts in the functions folder
sapply(list.files('Functions', full.names = T), source)

##The dataset with the default structure
data<-read.csv('Data/SampleData.csv')

##Run the seven regression models using 5 replicates per model and sampling 100
##observations per replicate. The number of generations for the Bayesian models
##is set to 1,000 but we generally need 20.000. The following line will ignore
##material identity (i.e. the mixed model is not included)
simulateAll(data, replicates=5, samples=100, generations=1000)

##If the user wants to analyze also the mixed model, use isMixed=T

simulateAll(data, replicates=5, samples=100, generations=1000, isMixed=T)


##Thanks for all your patience, Hannah! Hope that this is much more organized!

