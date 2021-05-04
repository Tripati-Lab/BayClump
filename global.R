
packages <- c("shiny", "shinydashboard", "dplyr", "xtable", "viridis", "data.table", "shinyBS", "plotly", "bibtex",
              "readxl", "lme4", "ggplot2", "R2jags", "IsoplotR", "chemCal", "investr", "MCMCvis", "ggplot2", "ggthemr",
              "ggridges", "data.table", "ggpubr", "DT", "knitcitations", "tictoc", "openxlsx", "pbapply", "deming", "pbmcapply")

new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, require, character.only = TRUE)

bibtex::write.bib(packages, file = "Rpackages.bib", append = FALSE, verbose = TRUE)

######################################################################################################

# Load necessary data

BayClump_calibration_template <- read.csv("BayClump_calibration_template.csv")
BayClump_reconstruction_template <- read.csv("BayClump_reconstruction_template.csv")
Petersen <- read.csv("Data/SampleData.csv")
Material <- NULL
Complete_Calibration_List <- NULL
# Load necessary scripts
sapply(list.files('Functions', full.names = T), source)
