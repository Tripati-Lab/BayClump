library(pbapply)
library(deming)
library(pbmcapply)
library(IsoplotR)
library(deming)
library(R2jags)
library(data.table)

sapply(list.files('Functions', full.names = T), source)

data<-read.csv('Data/SampleData.csv')
simulateAll(data, replicates=5, samples=100, generations=1000, isMixed=F)
simulateAll(data, replicates=5, samples=100, generations=1000, isMixed=T)
