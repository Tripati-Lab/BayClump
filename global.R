if("tidyverse" %in% installed.packages() == TRUE) {remove.packages("tidyverse")} # Tidyverse breaks the app by calling sample.int instead of sample if it's even installed

# To work in shinyapps.io, each package has be loaded individually instead of with lapply
library(shiny)
library(shinydashboard)
library(plyr)
library(dplyr)
library(xtable)
library(viridis)
library(data.table)
library(plotly)
library(bibtex)
library(readxl)
library(lme4)
library(R2jags)
library(IsoplotR)
library(investr)
library(DT)
library(knitcitations)
library(deming)
library(shinyBS)
library(ggridges)
library(ggpubr)
library(parallel)
library(openxlsx)
library(bib2df)
#library(loo)
#library(coda)

# Create automatic bibliography containing loaded packages and their dependencies
bibtex::write.bib(loadedNamespaces(), file = "Rpackages.bib", append = FALSE, verbose = TRUE)

######################################################################################################

# Load necessary data

BayClump_calibration_template <<- read.csv("Data/BayClump_calibration_template.csv")
BayClump_reconstruction_template <<- read.csv("Data/BayClump_reconstruction_template.csv")
Petersen <<- read.csv("Data/SampleData.csv") # PETERSEN IN 10^6/T^2
Anderson <<- read.csv("Data/SampleData2.csv") # ANDERSON IN 10^6/T^2
PetersenAnderson <<- rbind(Petersen, Anderson) # BOTH

# Create empty objects to receive data
Material <- NULL
Complete_Calibration_List <- NULL

# Improve reproducibility
set.seed(4)

# Load necessary functions
sapply(list.files('Functions', full.names = T), source)


