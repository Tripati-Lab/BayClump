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
library(doBy)

#Read the relevant scripts in the functions folder
sapply(list.files('Functions', full.names = T), source)

##Load the calibration dataset
calData<-read.csv('Data/SampleData.csv')
calData$T2 <- (10^6)/(calData$Temperature + 273.15)^2 

##Pipeline
PipCriteria <- read_excel("~/MEGA/Projects/1_ClumpedIsotopes/PredSummary_Jun18.xlsx", 
                                sheet = "selected")


##Modestou
targetD47<-c(0.714,
             0.715,
             0.723,
             0.714,
             0.713,
             0.708,
             0.712,
             0.706,
             0.705,
             0.709,
             0.705,
             0.706,
             0.695,
             0.71)
error_targetD47<-c(0.0067,
                   0.0067,
                   0.0066,
                   0.0067,
                   0.0067,
                   0.007,
                   0.0064,
                   0.0065,
                   0.007,
                   0.0071,
                   0.0068,
                   0.0066,
                   0.0071,
                   0.0067)



infTempBayesian<-clumpipe(calData=calData,
         PipCriteria=PipCriteria, 
         targetD47=targetD47, 
         error_targetD47=error_targetD47, 
         nrep=1000,
         BayesianOnly=T)

infTempBest<-clumpipe(calData=calData,
                  PipCriteria=PipCriteria, 
                  targetD47=targetD47, 
                  error_targetD47=error_targetD47, 
                  nrep=1000,
                  BayesianOnly=F)



