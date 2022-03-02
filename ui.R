# Define UI for application

header <-  dashboardHeader(title = "BayClump: Bayesian methods for clumped isotope paleothermometry",
                           titleWidth = 700
)



body <- dashboardBody(
  tags$head(tags$style(HTML('
      .form-group, .selectize-control {
           margin-bottom: 3px;
      }
      .form-group, .selectize-control {
           margin-top: 3px;
      }     
      .form-group, .selectize-control {
           margin-left: 3px;
      }
      .form-group, .selectize-control {
           margin-right: 3px;
      }
      .box-body {
          padding-bottom: 3px;
      }
      .box-body {
          padding-top: 2px;
      }
      .box {margin: 2px;}'
  )
  )
  ),
  
  ################# Tabs #################
  tabItems(
    tabItem(tabName = "calibration",
            fluidRow(
              box(width = 4, 
                  title = "Step 1: Calibration Options", solidHeader = TRUE,
                  column(12,
                         radioButtons("calset", "Select calibration set",
                                      choices = c('Román-Palacios et al. - Model 1' = 'model1',
                                                  'Román-Palacios et al. - Model 2' = 'model2',
                                                  "Use Román-Palacios et al. - Models 1 & 2" = 'model1and2',
                                                  "Upload my own calibration data" = 'mycal',
                                                  "Use all" = "all")
                         ),
                         
                         # Download template
                         downloadButton("BayClump_cal_temp", label = "Download calibration data template"),
                         
                         # Upload data
                         fileInput("calibrationdata", "Select calibration data file", accept = ".csv" 
                         ),
                         
                         # Uncertainties
                      #   radioButtons("uncertainties", " ",
                      #                c("My data contain uncertainties" = "myuncertainties",
                      #                  "Use uncertainties from Daëron et al., 2020" = "usedaeron")
                      #   ),
                         
                         # Reference frame
                         div(style="display:inline-block; vertical-align:top;", 
                             radioButtons("refframe", "Reference frame",
                                          c("Use the Intercarb carbon dioxide equilibrium scale at 90°C (I-CDES 90)" = "icdes90",
                                            "I am using my own calibration data in a different reference frame" = "myown")
                             ),                    
                             bsTooltip('refframe', "I-CDES 90 must be chosen if using calibration data from Román-Palacios et al.",
                                       placement = "bottom", trigger = "hover",
                                       options = NULL)
                         ),
                      
                         # Misc options
                         tags$b("Miscellaneous options"),
                         checkboxInput("scale" ,'Scale data')
                  )
                  
              ),
              
              # Model selection
              box(width = 4,
                  title = "Step 2: Select Models", solidHeader = TRUE,
                  column(12, "For help choosing an appropriate number of bootstrap replicates or the temperature range for CI estimation, see the User Manual",
                         numericInput("replication", label = "Number of bootstrap replicates for every model", 
                                      1000, min = 2, max = 10000),
                        sliderInput("range", label = HTML(paste0("Range for temperature to use for CI estimation (10",tags$sup("6"),"/T",tags$sup("2"),")")),
                                    min = 0, max = 30, value = c(1, 14)),
                        numericInput("generations", label = "Number of iterations for Bayesian models", 
                                     20000, min = 1000, max = 1000000),
                        checkboxInput("multicore", "Multicore", TRUE),
                         uiOutput("myList"),
                         selectInput("priors", label = "Bayesian priors", 
                                    choices = c("Informative", "Diffuse"), selected = "Informative"),
                         checkboxInput("simulateLM_measured", "Linear model", FALSE),
                         checkboxInput("simulateLM_inverseweights", "Inverse weighted linear model", FALSE),
                         checkboxInput("simulateYork_measured", "York regression", FALSE),
                         checkboxInput("simulateDeming", "Deming regression", FALSE),
                         checkboxInput("simulateBLM_measuredMaterial", "Bayesian simple linear model", FALSE),
                         checkboxInput("simulateBLMM_measuredMaterial", "Bayesian mixed model", FALSE),
                         bsTooltip('simulateBLM_measuredMaterial', "Running Bayesian models can take a few minutes. Please be patient.",
                                   placement = "bottom", trigger = "hover",
                                   options = NULL),
                         bsTooltip('simulateBLMM_measuredMaterial', "Running Bayesian models can take a few minutes. Please be patient.",
                                   placement = "bottom", trigger = "hover",
                                   options = NULL),
                         
                         # Summary stats panel
                         tableOutput("contents"),
                         
                         # Run models
                         div(style="display:inline-block; vertical-align:top;", 
                             actionButton('runmods', "Run selected models", 
                                          icon = icon("cogs", lib = "font-awesome")
                             )
                         ),
                         div(style="display:inline-block; vertical-align:top;", 
                             actionButton('reset', "Reset ALL", 
                                          icon = icon("trash-alt", lib = "font-awesome"))),
                         bsTooltip('reset', "Warning: This resets EVERYTHING, including data inputs",
                                   placement = "bottom", trigger = "hover",
                                   options = NULL),
                         verbatimTextOutput("modresults")
                         
                  )
              ),
              box(width = 4,
                  title = "Step 3: Calibration Output", solidHeader = TRUE,
                  column(12, "Truncated output from each selected model",
                         tags$h4("Linear model"),
                         verbatimTextOutput("lmcal", placeholder = TRUE),
                         tags$h4("Inverse weighted linear model"),
                         verbatimTextOutput("lminversecal", placeholder = TRUE),
                         tags$h4("York regression"),
                         verbatimTextOutput("york", placeholder = TRUE),
                         tags$h4("Deming regression"),
                         verbatimTextOutput("deming", placeholder = TRUE),
                         tags$h4("Bayesian simple linear model - no errors"),
                         verbatimTextOutput("blinnoerr", placeholder = TRUE),
                         tags$h4("Bayesian simple linear model - with errors"),
                         verbatimTextOutput("blinwerr", placeholder = TRUE),
                         tags$h4("Bayesian mixed model - with errors"),
                         verbatimTextOutput("blinmwerr", placeholder = TRUE)
                  ),
                  
                  # Download all calibration data
                  downloadButton("downloadcalibrations", label = "Download full calibration output"),
                  downloadButton("downloadBayesian", label = "Download raw results for Bayesian models"),
                  downloadButton("downloadPosteriorCalibration", label = "Download posterior one Bayesian replicate"),
                  downloadButton("downloadPriorsCalibration", label = "Download priors")
              )
            )
    ),
    
    #Calibration plot tab
    tabItem(tabName = "plots",
            fluidRow(
              box(width = 12,
                  column(12,
                         title = "Calibration plots", solidHeader = TRUE,
                         plotlyOutput("rawcaldata")
                         
                  ))),
            fluidRow(
              box(width = 12,
                  column(12,
                         plotlyOutput("lmcalibration")
                  ))),
            
            fluidRow(
              box(width = 12,
                  column(12,
                         plotlyOutput("lminversecalibration")
                  ))),
            
            fluidRow(
              box(width = 12,
                  column(12,
                         plotlyOutput("yorkcalibration")
                  ))),
            
            fluidRow(
              box(width = 12,
                  column(12,
                         plotlyOutput("demingcalibration")
                  ))),
            
            fluidRow(
              box(width = 12,
                  column(12,
                         plotlyOutput("bayeslincalibration")
                  ))),
            
            fluidRow(
              box(width = 12,
                  column(12,
                         plotlyOutput("bayesmixedcalibration")
                  )))
            
            
    ),
    
    
    #Reconstruction tab
    tabItem(tabName = "reconstruction",
            fluidRow(
              box(width = 5, 
                  title = "Step 1: Reconstruction setup", solidHeader = TRUE,
                  column(12, tags$b("Parameters for the selected models are automatically transferred from the Calibration tab to this tab. If you are interested in running classic reconstructions, please calibrate the relevant models calibration before."),
                         tags$br(),
                         # Download templates
                         downloadButton("BayClump_reconstruction_template.csv", label = "Download reconstruction data template"),
                         
                         # Upload data
                         fileInput("reconstructiondata", "Select reconstruction data file", accept = ".csv" 
                         ),
                         
                         # Summary stats panel
                         tableOutput("contents2"),
                         numericInput("replicationRec", label = "Number of bootstrap replicates for every model", 
                                      100, min = 2, max = 10000),
                         tags$b("Models to run:"),
                         checkboxInput("simulateLM_measuredRec", "Linear model", FALSE),
                         checkboxInput("simulateLM_inverseweightsRec", "Inverse weighted linear model", FALSE),
                         checkboxInput("simulateYork_measuredRec", "York regression", FALSE),
                         checkboxInput("simulateDemingRec", "Deming regression", FALSE),
                         checkboxInput("simulateBLM_measuredMaterialRec", "Bayesian simple linear model", FALSE),
                         checkboxInput("simulateBLMM_measuredMaterialRec", "Bayesian mixed model", FALSE),
                         tags$b("Factor in parameter uncertainty?"),
                         checkboxInput("AccountErrorDataset", "Yes"),
                         bsTooltip('bayesianPredictions', "This can take several minutes for large datasets",
                                   placement = "bottom", trigger = "hover",
                                   options = NULL),
                         
                         div(style="display:inline-block; vertical-align:top;", 
                             actionButton('runrec', "Run reconstructions", 
                                          icon = icon("cogs", lib = "font-awesome")
                             )
                         ),
                         verbatimTextOutput("recresults"),
                         # Download all reconstruction data
                         downloadButton("downloadreconstructions", label = "Download reconstruction output"),
                         downloadButton("downloadreconstructionsPosterior", label = "Download posterior reconstruction output"),
                         downloadButton("downloadPriorsReconstruction", label = "Download priors")
                         
                  )
              ),
              box(width = 7,
                  title = "Step 2: Temperature reconstructions", solidHeader = TRUE,
                  column(12, "Truncated output from each selected model",
                         tableOutput("lmrecswun"),
                         tableOutput("lminverserecswun"),
                         tableOutput("yorkrecswun"),
                         tableOutput("demingrecswun"),
                         tableOutput("Bpredictions"),
                         tableOutput("BpredictionsErrors"),
                         tableOutput("BpredictionsBLMM")
                  )
              )
            ),  
            
    ),
    
    #User manual
    tabItem(tabName = "usermanual",
            fluidRow(
              column(12,
                     withMathJax(includeMarkdown("userguide.Rmd"))
              ))),
    
    #BayWatch
    tabItem(tabName = "demo",
            fluidRow(
              column(12,
                     "Recorded demo here"))),
    
    #Citations tab
    tabItem(tabName = "citations",
            fluidRow(
              column(12,
                     DT::dataTableOutput("bibTable")
              ))),
    
    #Manuscript
    tabItem(tabName = "manuscript",
            fluidRow(
              column(12, "A static link to the publisher's website will be embedded here upon acceptance for publication"
                     #includeHTML("Our manuscript here")
              ))),
    
    #Download
    tabItem(tabName = "download",
            fluidRow(
              column(12,
                     "BayClump will be made available as a standalone Electron app upon acceptance for publication, and a download link will be provided here")))
  )
)





############ Dashboard Sidebar ############


sidebar <- dashboardSidebar(width = 200,
                            sidebarMenu(
                              menuItem("Calibrations", tabName = "calibration", 
                                       icon = icon("drafting-compass", lib = "font-awesome")
                              ),
                              
                              menuItem("Calibration Plots", tabName = "plots", 
                                       icon =icon("chart-bar", lib = "font-awesome")
                              ),
                              
                              menuItem("Reconstructions", tabName = "reconstruction", 
                                       icon = icon("chart-area", lib = "font-awesome")
                              ),
                              
                              menuItem("User Manual", tabName = "usermanual", 
                                       icon = icon("question-circle", lib = "font-awesome")
                              ),
                              
                              menuItem("BayWatch", tabName = "demo", 
                                       icon = icon("glasses", lib = "font-awesome")
                              ),
                              
                              menuItem("Citations", tabName = "citations", 
                                       icon = icon("align-right", lib = "font-awesome")
                              ),
                              
                              menuItem("Manuscript", tabName = "manuscript", 
                                       icon = icon("scroll", lib = "font-awesome")
                              ),
                              
                              menuItem("Download BayClump", tabName = "download", 
                                       icon = icon("laptop-code", lib = "font-awesome")
                              ),
                              
                              menuItem("Contacts", icon = icon("address-card", lib = "font-awesome"),
                                       a(actionButton(inputId = "email1", label = "Email", 
                                                      icon = icon("envelope", lib = "font-awesome")),
                                         href="mailto:hannah.carroll@paleoeco.org"),
                                       a(actionButton(inputId = "github", label = "GitHub", 
                                                      icon = icon("github", lib = "font-awesome"),
                                                      onclick ="window.open('https://github.com/Tripati-Lab/BayClump/', '_blank')")),
                                       a(actionButton(inputId = "bugs", label = "Submit a Bug Report", 
                                                      icon = icon("bug", lib = "font-awesome"),
                                                      onclick ="window.open('https://github.com/Tripati-Lab/BayClump/issues', '_blank')")
                                       )
                              )
                            )
)


########### Dashboard page ####################
dashboardPage(
  header,
  sidebar,
  body
)


