if("tidyverse" %in% installed.packages() == TRUE) {remove.packages("tidyverse")} # Tidyverse breaks the app by calling sample.int instead of sample if it's even installed

# To work in shinyapps.io, each package has be loaded individually instead of with lapply
library(shiny)
library(shinydashboard)
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

packages <- c("shiny", "shinydashboard", "dplyr", "xtable", "viridis", "data.table", "shinyBS", "plotly", "bibtex",
              "readxl", "lme4", "ggplot2", "rjags", "R2jags", "IsoplotR", "chemCal", "investr", "MCMCvis",
              "ggridges", "ggpubr", "DT", "knitcitations", "openxlsx", "pbapply", "deming", "pbmcapply", "clumpedr")


#new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)
#lapply(packages, library, character.only = TRUE)

bibtex::write.bib(packages, file = "Rpackages.bib", append = FALSE, verbose = TRUE)

######################################################################################################

# Load necessary data

BayClump_calibration_template <- read.csv("BayClump_calibration_template.csv")
BayClump_reconstruction_template <- read.csv("BayClump_reconstruction_template.csv")
Petersen <- read.csv("Data/SampleData.csv") # PETERSEN PLACEHOLDER
Anderson <- read.csv("Data/SampleData2.csv") # ANDERSON PLACEHOLDER
PetersenAnderson <- rbind(Petersen, Anderson)
Material <- NULL
Complete_Calibration_List <- NULL
# Load necessary scripts
sapply(list.files('Functions', full.names = T), source)

