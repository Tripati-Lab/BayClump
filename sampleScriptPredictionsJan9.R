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
library(dplyr)

#Potentially new packages!
library(clumpedr) ##devtools::install_github("isoverse/clumpedr", ref = "dev")
library(data.table)


#Read the relevant scripts in the functions folder
sapply(list.files('Functions', full.names = T), source)

##Load the calibration dataset
calData<-read.csv('Data/SampleData.csv')
calData$T2 <- calData$Temperature
calData$TempError<-ifelse(calData$TempError==0, 10-2, calData$TempError)
calData<-calData[complete.cases(calData),]


################NEW APPROACH

##Non-Bayesian are based on individual replicates per sample
targetD47<-c(0.71,0.72, 0.71)
predictTclassic(calData, targety=targetD47, model='wlm', replicates=5, bootDataset=T, onlyMedian=T)
predictTclassic(calData, targety=targetD47, model='wlm', replicates=5, bootDataset=T, onlyMedian=F)
predictTclassic(calData, targety=targetD47, model='wlm', replicates=5, bootDataset=F, onlyMedian=F)

predictTclassic(calData, targety=targetD47, model='lm', replicates=5, bootDataset=T, onlyMedian=T)
predictTclassic(calData, targety=targetD47, model='lm', replicates=5, bootDataset=T, onlyMedian=F)
predictTclassic(calData, targety=targetD47, model='lm', replicates=5, bootDataset=F, onlyMedian=F)

predictTclassic(calData, targety=targetD47, model='York', replicates=5, bootDataset=T, onlyMedian=T)
predictTclassic(calData, targety=targetD47, model='York', replicates=5, bootDataset=T, onlyMedian=F)
predictTclassic(calData, targety=targetD47, model='York', replicates=5, bootDataset=F, onlyMedian=F)

predictTclassic(calData, targety=targetD47, model='Deming', replicates=5, bootDataset=T, onlyMedian=T)
predictTclassic(calData, targety=targetD47, model='Deming', replicates=5, bootDataset=T, onlyMedian=F)
predictTclassic(calData, targety=targetD47, model='Deming', replicates=5, bootDataset=F, onlyMedian=F)


##Bayesian use mean and error for a given sample
error_targetD47 <- c(0.01, 0.02,0.02)
material <- c(1,2,3)
predictTcBayes(calibrationData=calData[sample(1:nrow(calData),nrow(calData)),], 
               data=cbind(D47=targetD47,error=error_targetD47), generations=20000, 
               hasMaterial=T)

predictTcBayes(calibrationData=calData[sample(1:nrow(calData),nrow(calData)),], 
               data=cbind(D47=targetD47,error=error_targetD47), generations=20000, 
               hasMaterial=F)

predictTcBayes_replicates(calData=calData, 
                          targetD47=targetD47, 
                          error_targetD47=error_targetD47, 
                          material=material, 
                          nrep=2, 
                          hasMaterial=F, 
                          generations=1000)

predictTcBayes_replicates(calData=calData, 
                          targetD47=targetD47, 
                          error_targetD47=error_targetD47, 
                          material=material, 
                          nrep=2, 
                          hasMaterial=T, 
                          generations=1000)



###########Classic approach

simulateYork_measured(data=calData, replicates=5)
simulateLM_measured(data=calData, replicates=5)
simulateDeming(data=calData, replicates=5)
simulateLM_inverseweights(data=calData, replicates=5)
simulateBLM_measuredMaterial(data=calData, replicates=5, generations=1000, isMixed=F)
simulateBLM_measuredMaterial(data=calData, replicates=5, generations=1000, isMixed=T)

