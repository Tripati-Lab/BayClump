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


##Non-Bayesian are based on individual replicates per sample
targetD47<-c(0.71,0.72, 0.71)
predictTclassic(calData, targety=targetD47, model='wlm')
predictTclassic(calData, targety=targetD47, model='lm')
predictTclassic(calData, targety=targetD47, model='York')
predictTclassic(calData, targety=targetD47, model='Deming')


##Bayesian use mean and error for a given sample
error_targetD47 <- c(0.01, 0.02,0.02)
material <- c(1,2,3)
predictTcBayes(calibrationData=calData[sample(1:nrow(calData),nrow(calData)),], 
               data=cbind(D47=targetD47,error=error_targetD47), generations=20000, 
               hasMaterial=T)

predictTcBayes(calibrationData=calData[sample(1:nrow(calData),nrow(calData)),], 
               data=cbind(D47=targetD47,error=error_targetD47), generations=20000, 
               hasMaterial=F)






