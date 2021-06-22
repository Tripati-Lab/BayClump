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
library(readxl)
library(ggthemr)
ggthemr('pale', layout = 'scientific', spacing = 2, type = 'inner')

#Read the relevant scripts in the functions folder
sapply(list.files('Functions', full.names = T), source)

##Load the calibration dataset
calData<-read.csv('Data/SampleData.csv')
calData$T2 <- (10^6)/(calData$Temperature + 273.15)^2 

##Pipeline
PipCriteria <- read_excel("~/MEGA/Projects/1_ClumpedIsotopes/PredSummary_Jun18.xlsx", 
                                sheet = "selected")


##Dataset
Miocene_Datasets <- read_excel("~/Downloads/Miocene Datasets.xlsx", 
                               sheet = "Modestou2020 (Averages)")

##Modestou------
targetD47<-Miocene_Datasets$`D47 (ICDES90)`
error_targetD47<-Miocene_Datasets$seD47

#Updated frame
PredUpdatedFrame<-predictTcNonBayes(data=cbind(targetD47,error_targetD47), 
                  slope=0.0449, 
                  slpcnf=0.001, 
                  intercept=0.167, 
                  intcnf=0.01)

PredOriginalFrame<-predictTcNonBayes(data=cbind(Miocene_Datasets$`D47 (Bernasconi Ref. Frame)`,error_targetD47), 
                                    slope=0.0449, 
                                    slpcnf=0.001, 
                                    intercept=0.167, 
                                    intcnf=0.01)

infTempBayesian<-clumpipe(calData=calData,
         PipCriteria=PipCriteria, 
         targetD47=targetD47, 
         error_targetD47=error_targetD47, 
         nrep=100,
         BayesianOnly=T)


infTempBest<-clumpipe(calData=calData,
                  PipCriteria=PipCriteria, 
                  targetD47=targetD47, 
                  error_targetD47=error_targetD47, 
                  nrep=1000,
                  BayesianOnly=F)


Data_complete <-rbindlist(list(cbind(group="OriginalFrame",PredOriginalFrame),
               cbind(group="NewFrame",PredUpdatedFrame),
               cbind(group="clumpipe",   infTempBest[,c(19,20,21:23)]),
               cbind(group="Bayesian",  infTempBayesian[,c(1,2,21:23)])), fill = T)

Data_complete$Age<- rep(Miocene_Datasets$Age, nrow(Data_complete)/nrow(Miocene_Datasets))
Data_complete$type<-ifelse(is.na(Data_complete$type), "Parameter Uncertainty", Data_complete$type)

write.csv(Data_complete, 'Modestou.predictions.csv')


ModestouPlot<-ggplot()+
  geom_ribbon(data=Data_complete[Data_complete$type== 'Parameter Uncertainty' & Data_complete$group != 'OriginalFrame' ,], aes(x=Age, y=Tc,ymin=lwr, ymax=upr,fill=group, color = group), alpha=0.1, linetype='dotted')+
  geom_ribbon(data=NULL,aes(ymin=infTempBest[,'Tc'], ymax=infTempBayesian[,'Tc'], x=
                              Miocene_Datasets$Age), fill='#FFDB58', alpha=0.9)+
  geom_line(data=Data_complete[Data_complete$type== 'Parameter Uncertainty' & Data_complete$group != 'OriginalFrame' ,], 
            aes(x=Age, y=Tc, color=group), size=1)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
  
pdf('Modestou.pdf', 7,5)
print(ModestouPlot)
dev.off()


ModestouPlotOriginal<-ggplot(data=Data_complete[Data_complete$group %in% c('OriginalFrame','NewFrame') ,])+
  geom_ribbon(aes(x=Age, y=Tc,ymin=lwr, ymax=upr,fill=type, color = type), alpha=0.1, linetype='dotted')+
  geom_line(aes(x=Age, y=Tc, color=type), size=1)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  facet_wrap(~group)


pdf('Modestou_IgnoreParameters.pdf', 10,5)
print(ModestouPlotOriginal)
dev.off()


###XXXXXX--------
Miocene_Datasets <- read_excel("~/Downloads/Miocene Datasets.xlsx", 
                               sheet = "Leutert2020 (Averages)")

##Leutert2020------
targetD47<-Miocene_Datasets$`D47 (ICDES90)`
error_targetD47<-Miocene_Datasets$seD47

#Updated frame
PredUpdatedFrame<-predictTcNonBayes(data=cbind(targetD47,error_targetD47), 
                                    slope=0.0449, 
                                    slpcnf=0.001, 
                                    intercept=0.167, 
                                    intcnf=0.01)

PredOriginalFrame<-predictTcNonBayes(data=cbind(Miocene_Datasets$`D47 (Bernasconi Ref. Frame)`,error_targetD47), 
                                     slope=0.0449, 
                                     slpcnf=0.001, 
                                     intercept=0.167, 
                                     intcnf=0.01)

infTempBayesian<-clumpipe(calData=calData,
                          PipCriteria=PipCriteria, 
                          targetD47=targetD47, 
                          error_targetD47=error_targetD47, 
                          nrep=100,
                          BayesianOnly=T)


infTempBest<-clumpipe(calData=calData,
                      PipCriteria=PipCriteria, 
                      targetD47=targetD47, 
                      error_targetD47=error_targetD47, 
                      nrep=1000,
                      BayesianOnly=F)


Data_complete <-rbindlist(list(cbind(group="OriginalFrame",PredOriginalFrame),
                               cbind(group="NewFrame",PredUpdatedFrame),
                               cbind(group="clumpipe",   infTempBest[,c(19,20,21:23)]),
                               cbind(group="Bayesian",  infTempBayesian[,c(1,2,21:23)])), fill = T)

Data_complete$Age<- rep(Miocene_Datasets$Age, nrow(Data_complete)/nrow(Miocene_Datasets))
Data_complete$type<-ifelse(is.na(Data_complete$type), "Parameter Uncertainty", Data_complete$type)


write.csv(Data_complete, 'Leutert2020.predictions.csv')



ModestouPlot<-ggplot()+
  geom_ribbon(data=Data_complete[Data_complete$type== 'Parameter Uncertainty' & Data_complete$group != 'OriginalFrame' ,], aes(x=Age, y=Tc,ymin=lwr, ymax=upr,fill=group, color = group), alpha=0.1, linetype='dotted')+
  geom_ribbon(data=NULL,aes(ymin=infTempBest[,'Tc'], ymax=infTempBayesian[,'Tc'], x=
                              Miocene_Datasets$Age), fill='#FFDB58', alpha=0.9)+
  geom_line(data=Data_complete[Data_complete$type== 'Parameter Uncertainty' & Data_complete$group != 'OriginalFrame' ,], 
            aes(x=Age, y=Tc, color=group), size=1)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))

pdf('Leutert2020.pdf', 7,5)
print(ModestouPlot)
dev.off()


ModestouPlotOriginal<-ggplot(data=Data_complete[Data_complete$group %in% c('OriginalFrame','NewFrame') ,])+
  geom_ribbon(aes(x=Age, y=Tc,ymin=lwr, ymax=upr,fill=type, color = type), alpha=0.1, linetype='dotted')+
  geom_line(aes(x=Age, y=Tc, color=type), size=1)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  facet_wrap(~group)


pdf('Leutert2020_IgnoreParameters.pdf', 10,5)
print(ModestouPlotOriginal)
dev.off()


