##Packages
library(pbapply)
library(deming)
library(pbmcapply)
library(IsoplotR)
library(deming)
library(R2jags)
library(data.table)
library(ggplot2)
library(plyr)

#Read the relevant scripts in the functions folder
sapply(list.files('Functions', full.names = T), source)

##The dataset with the default structure
data<-read.csv('Data/SampleData.csv')
data$TempError<-ifelse(data$TempError==0, 10-2, data$TempError)
data<-data[complete.cases(data),]


#Run each model independently using 5 replicates---------
 
#100 observations per replicate. The number of generations for 
#the Bayesian models is set to 1,000 but we generally need 20.000.
#isMixed won't fit the mixed model (no Material column needed).

simulateYork_measured(data=data, replicates=5)
simulateLM_measured(data=data, replicates=5)
simulateDeming(data=data, replicates=5)
simulateLM_inverseweights(data=data, replicates=5)
simulateBLM_measuredMaterial(data=data, replicates=5, generations=1000, isMixed=F)
simulateBLM_measuredMaterial(data=data, replicates=5, generations=1000, isMixed=T)

##Now, run the seven regression models using 5 replicates per model and sampling 100
##observations per replicate. 
simulateAll(data, replicates=5, samples=100, generations=1000)

##If the user wants to analyze also the mixed model, use isMixed=T
simulateAll(data, replicates=5, samples=100, generations=1000, isMixed=T)



##Regression confidence intervals ---------

##For a data.frame of intercept and slope per replicate

datReplicatesSingle<-simulateYork_measured(data=data, replicates=5, samples=100)
outCI<-RegressionSingleCI(data=datReplicatesSingle, from=5, to=15, length.out=100)
outCI$dataUncertainties



##Get temperature reconstructions--------

#First, one that only needs a data.frame of intercept, 
#slope per replicate, along with target D47, measurement error in target D47 and,
#number of measurements used to estimate error in the target D47.

datReplicatesYork<-simulateYork_measured(data=data, replicates=5, samples=100)
TempRecYork<-TempReconstructionCoefficients(datReplicatesYork, 
                                     targety=c(0.6,0.7), 
                                     errory=c(0.005, 0.01), 
                                     nmeasurements=c(2,3))

TempRecYork$summary  ##The summary per target D47
TempRecYork$replicates #If people want to see the analyses per replicate

#Second, Bayesian predictions. We only need to provide the dataset that was used
#for an initial calibration step in addition to target D47, measurement error in 
#target D47 and, number of measurements used to estimate error in the target D47.


##If we account for different materials in the dataset, we have to indicate their
##identities (relative to the first dataset). hasMaterial will run a mixed model.
TempRecBayesian<-TempReconstructionBayesian(dataCalibration=data,
                                 targety=c(0.6,0.7), 
                                 errory=c(0.005, 0.01), 
                                 materialPred=c(1,2), 
                                 nmeasurements=c(2,3),
                                 hasMaterial=T)

TempRecBayesian$summary
TempRecBayesian$replicates

##If we only have a single material but still would like to run the mixed, the material
##column in the dataset should be set to a constant number (same as in the materialPred)

data2<-data
data2$Material<-1
TempRecBayesian<-TempReconstructionBayesian(data2,
                                            targety=c(0.6,0.7), 
                                            errory=c(0.005, 0.01), 
                                            nmeasurements=c(2,3),
                                            materialPred=c(1,1), 
                                            hasMaterial=T)

##The mixed model can also be excluded from the analyses

TempRecBayesian<-TempReconstructionBayesian(data,
                                            targety=c(0.6,0.7), 
                                            errory=c(0.005, 0.01), 
                                            nmeasurements=c(2,3),
                                            hasMaterial=F)


      