
packages <- c("shiny", "shinydashboard", "leaflet", "dplyr", "xtable", "viridis", "data.table", "shinyBS", "plotly")

lapply(packages, require, character.only = TRUE)

######################################################################################################

# Load necessary data

BayClump_template <- read.csv("BayClump_template.csv")
