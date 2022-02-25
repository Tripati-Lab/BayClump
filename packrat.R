## Install the package from CRAN
install.packages("packrat")

## Initialize packrat
packrat::init("~/Documents/Research/UCLA/Tripati lab/Bayesian calibrations/BayClump-main")

## Later on, install some package
## It will be installed in ~/path/to/your/project/packrat/lib
install.packages("shiny")
install.packages("shinydashboard")
install.packages("plyr")
install.packages("dplyr")
install.packages("xtable")
install.packages("parallel")
install.packages("viridis")
install.packages("data.table")
install.packages("plotly")
install.packages("bibtex")
install.packages("readxl")
install.packages("lme4")
install.packages("R2jags")
install.packages("IsoplotR")
install.packages("investr")
install.packages("DT")
install.packages("knitcitations")
install.packages("pbapply")
install.packages("deming")
install.packages("shinyBS")
install.packages("ggridges")
install.packages("ggpubr")
install.packages("pbmcapply")
install.packages("openxlsx")
install.packages("devtools")
install.packages("bib2df")

## Take a snapshot of installed packages
packrat::snapshot()
