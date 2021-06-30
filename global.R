if("tidyverse" %in% installed.packages() == TRUE) {remove.packages("tidyverse")} # Tidyverse breaks the app by calling sample.int instead of sample if it's even installed

# To work in shinyapps.io, each package has be loaded individually instead of with lapply
library(shiny)
library(shinydashboard)
library(isoreader)
library(dplyr)
library(xtable)
library(viridis)
library(data.table)
library(plotly)
library(bibtex)
library(readxl)
library(lme4)
library(rjags)
library(R2jags)
library(IsoplotR)
library(chemCal)
library(investr)
library(MCMCvis)
library(DT)
library(knitcitations)
library(pbapply)
library(deming)
library(shinyBS)
library(ggridges)
library(ggpubr)
library(pbmcapply)
library(openxlsx)
library(clumpedr)
library(bib2df)

# Create automatic bibliography containing loaded packages and their dependencies
bibtex::write.bib(loadedNamespaces(), file = "Rpackages.bib", append = FALSE, verbose = TRUE)

######################################################################################################

# Load necessary data

BayClump_calibration_template <- read.csv("BayClump_calibration_template.csv")
BayClump_reconstruction_template <- read.csv("BayClump_reconstruction_template.csv")
Petersen <- read.csv("Data/SampleData.csv") # PETERSEN PLACEHOLDER
Anderson <- read.csv("Data/SampleData2.csv") # ANDERSON PLACEHOLDER
PetersenAnderson <- rbind(Petersen, Anderson)

# Create empty objects to receive data
Material <- NULL
Complete_Calibration_List <- NULL

# Improve reproducibility
set.seed(4)

# Load necessary functions
sapply(list.files('Functions', full.names = T), source)


