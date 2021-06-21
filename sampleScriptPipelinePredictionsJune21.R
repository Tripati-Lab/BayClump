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

##Closest functions
dist <- function(df2, D47error, TError, bestCase=T){
  dt <- data.table((df2$Dataset_D47Error-D47error)^2+(df2$Dataset_TError-TError)^2)
  return(if(bestCase){which.min(dt$V1)}else{which.max(dt$V1)})
}

distD47 <- function(df2, TargetD47error, Target_D47, bestCase=T){
  dt <- data.table((df2$Target_D47Error-TargetD47error)^2+(df2$Target_D47-Target_D47)^2)
  return(if(bestCase){which.min(dt$V1)}else{which.max(dt$V1)})
}

#Read the relevant scripts in the functions folder
sapply(list.files('Functions', full.names = T), source)

##Load the calibration dataset
calData<-read.csv('Data/SampleData.csv')
calData$T2 <- (10^6)/(calData$Temperature + 273.15)^2 

##Pipeline
PipCriteria <- read_excel("~/MEGA/Projects/1_ClumpedIsotopes/PredSummary_Jun18.xlsx", 
                                sheet = "selected")

##The user need to provide more than one target D47 and an error
targetD47<-c(0.71,0.72, 0.89,0.61); error_targetD47<-c(0.02,0.03,0.02,0.03)


##Find overall error scenario for the dataset

targetScenario<-unlist(PipCriteria[dist(df2=PipCriteria, 
     D47error=mean(calData$D47error), 
     TError=mean(calData$TempError), bestCase = T),1])


subPipCriteria<-PipCriteria[PipCriteria$Error_scenario ==targetScenario, ]

##Find model for each target D47

compilationModels<-cbind.data.frame(targetD47, error_targetD47, do.call(rbind,lapply(1:length(targetD47), function(x){
  #Don't change the bestCase argument
  selRow<-distD47(subPipCriteria, TargetD47error=error_targetD47[x], Target_D47=targetD47[x], bestCase = T)
  subPipCriteria[selRow,]
})))




##Bayesian predictions

BP<-compilationModels[compilationModels$PredictionType == 'Bayesian predictions',]
NBP<-compilationModels[compilationModels$PredictionType != 'Bayesian predictions',]


##Perform non-Bayesian predictions
prediction_NBP<-do.call(rbind,lapply(1:nrow(NBP), function(x){
  SM<-NBP[x,'ModelSelected']
  
  ##Run the model
  if(SM == 'LM'){
   valsReps<- simulateLM_measured(data=calData, replicates=1000, samples=50)
    
  }
  if(SM=='York'){
    valsReps<- simulateYork_measured(data=calData, replicates=1000, samples=50)
    
  }
  if(SM=='Deming'){
    valsReps<-  simulateDeming(data=calData, replicates=1000, samples=50)
    
  }
  
  if(SM=='invweLM'){
    valsReps<- simulateLM_inverseweights(data=calData, replicates=1000, samples=50)
    
  }
  
  if(length(grep('BLM',SM))>0){
    ##Just do the one with errors for now
    valsReps<-simulateBLM_measuredMaterial(data=calData, replicates=1000, samples=50, generations=10000, isMixed=F)[[1]]
  }
  
  ##Get the predictions
  
  predictTcNonBayes(data=cbind(NBP$targetD47[x],NBP$error_targetD47[x]), 
                    slope=median(valsReps$slope), 
                    slpcnf=CItoSE(quantile(valsReps$slope, 0.975), quantile(valsReps$slope, 0.025)), 
                    intercept=median(valsReps$intercept), 
                    intcnf=CItoSE(quantile(valsReps$intercept, 0.975), quantile(valsReps$intercept, 0.025)))[1,]
  
}))


##Perform Bayesian predictions

singleRep<-function(i) {predictTcBayes(calibrationData=calData, 
               data=cbind(BP$targetD47,BP$error_targetD47),
               generations=5000, 
               hasMaterial=F, onlyWithinBayesian=T)
}

totalRep<-pbmclapply(1:5,singleRep, mc.cores = 4)







##Fit the best model
targetScenario$PredictionType
targetScenario$ModelSelected






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


