
packages <- c("shiny", "shinydashboard", "leaflet", "dplyr", "xtable", "viridis", "data.table", "shinyBS", "plotly",
              "readxl", "lme4", "ggplot2", "R2jags", "IsoplotR", "chemCal", "investr", "MCMCvis", "ggplot2", "ggthemr",
              "ggridges", "data.table", "ggpubr", "DT")

lapply(packages, require, character.only = TRUE)

######################################################################################################

# Load necessary data

BayClump_template <- read.csv("BayClump_template.csv")
Petersen <- read.csv("Data/SampleData.csv")

# Load necessary scripts
sapply(list.files('Functions', full.names = T), source)
