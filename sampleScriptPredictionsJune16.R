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

calData<-calData[complete.cases(calData),]

##Two replicates of the same sample
targetD47<-c(0.71,0.72, 0.71)
n<-predictTclassic(calData, targetD47, model='wlm')




################################################
#####Fit only a single replicate per model######
################################################

##If the user wants to get predictions given a LM
reg<-summary(lm(calData$D47 ~ calData$T2))


predictTcNonBayes(data=cbind(targetD47,error_targetD47), 
                  slope=reg$coefficients[2,1], 
                  slpcnf=reg$coefficients[2,2], 
                  intercept=reg$coefficients[1,1], 
                  intcnf=reg$coefficients[1,2])


##If the user wants to get predictions given a York model
Reg<-york(cbind.data.frame(calData$T2, 
                           calData$TempError, 
                           calData$D47,
                           calData$D47error))
predictTcNonBayes(cbind(targetD47,error_targetD47), 
                  Reg$b[1], Reg$b[2], Reg$a[1], Reg$a[2] )

##If the user wants to get predictions given some a weighted regression model
Reg<-summary(lm( D47 ~ T2 ,  calData, weights = D47error))
WeightedPred<-predictTcNonBayes(cbind(targetD47,error_targetD47), 
                                Reg$coefficients[2,1],Reg$coefficients[2,2], Reg$coefficients[1,1], 
                                Reg$coefficients[1,2] )

##If the user wants to get predictions given some a Deming regression model
Reg<-deming(D47 ~ T2, calData, xstd= 1/calData$TempError, ystd= 1/calData$D47error)
predictTcNonBayes(cbind(targetD47,error_targetD47), Reg$coefficients[2],
                  as.numeric(strsplit(capture.output(Reg)[8], ' ')[[1]][7]),
                  Reg$coefficients[1], 
                  as.numeric(strsplit(capture.output(Reg)[7], ' ')[[1]][3] ))

##If the user wants to get predictions under a Bayesian framework
##(simple replicate; inside and outside a Bayesian framework)
predictTcBayes(calibrationData=calData, 
               data=cbind(targetD47,error_targetD47),
               generations=5000, 
               hasMaterial=F)


######What if the user wants to get predictions based on especific regression parameters?######
######Imagine the users want to use the output of the calibration step or their own dataset?##########

Slope=0.0449; SlopeSE=0.001; Intercept=0.167; InterceptSE=0.01
predictTcNonBayes(data=cbind(targetD47,error_targetD47), 
                  slope=Slope, 
                  slpcnf=SlopeSE, 
                  intercept=Intercept, 
                  intcnf=InterceptSE)
##If there is no uncertainty in regression parameters, use 0

########Please don't implement this one yet!!###########
########It's too take consuming#########################


##If the user wants to get predictions under a Bayesian framework
##(multiple replicates; inside a Bayesian framework)

single_rep<-function(i){
  tryCatch({
    res<- predictTcBayes(calibrationData=calData[sample(1:nrow(calData),nrow(calData)),], 
                         data=cbind(targetD47,error_targetD47), generations=20000, 
                         hasMaterial=F, onlyWithinBayesian = T)
    cbind.data.frame(replicate=i, res)
  }, error=function(e){})
}

pbmclapply(1:2, single_rep, mc.cores = 4)


