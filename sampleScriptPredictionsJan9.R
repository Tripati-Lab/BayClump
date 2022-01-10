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
#library(clumpedr) ##devtools::install_github("isoverse/clumpedr", ref = "dev")
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
predictTcBayes(calibrationData=calData, 
               data=cbind(D47=targetD47,error=error_targetD47), generations=1000, 
               hasMaterial=T, bootDataset=T, onlyMedian=T, replicates = 5)
predictTcBayes(calibrationData=calData, 
               data=cbind(D47=targetD47,error=error_targetD47), generations=1000, 
               hasMaterial=T, bootDataset=T, onlyMedian=T, replicates = 5)
predictTcBayes(calibrationData=calData, 
               data=cbind(D47=targetD47,error=error_targetD47), generations=1000, 
               hasMaterial=T, bootDataset=T, onlyMedian=T, replicates = 5)


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

a<- simulateYork_measured(data=calData, replicates=5)
b<- simulateLM_measured(data=calData, replicates=5)
c<- simulateDeming(data=calData, replicates=5)
d<- simulateLM_inverseweights(data=calData, replicates=5)
e<- simulateBLM_measuredMaterial(data=calData, replicates=5, generations=1000, isMixed=F)
f<- simulateBLM_measuredMaterial(data=calData, replicates=5, generations=1000, isMixed=T)

classicCalibration(reps = a, targetD47=targetD47, error_targetD47=error_targetD47)
classicCalibration(reps = b, targetD47=targetD47, error_targetD47=error_targetD47)
classicCalibration(reps = c, targetD47=targetD47, error_targetD47=error_targetD47)
classicCalibration(reps = d, targetD47=targetD47, error_targetD47=error_targetD47)
classicCalibration(reps = e$BLM_Measured_errors, targetD47=targetD47, error_targetD47=error_targetD47)
classicCalibration(reps = e$BLM_Measured_no_errors, targetD47=targetD47, error_targetD47=error_targetD47)
classicCalibration(reps = f$BLMM_Measured_errors, targetD47=targetD47, error_targetD47=error_targetD47)





