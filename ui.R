# Define UI for application

header <-  dashboardHeader(title = "BayClump: Bayesian methods for clumped isotope paleothermometry",
                           titleWidth = 700)

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
                    radioButtons("uncertainties", " ",
                                 c("My data contain uncertainties" = "myuncertainties",
                                   "Use uncertainties from Daëron et al., 2020" = "usedaeron")
                                 ),
                    
                    # Reference frame
                    div(style="display:inline-block; vertical-align:top;", 
                        radioButtons("refframe", "Reference frame",
                                c("Use the carbon dioxide equilibrium scale (CDES)" = "cdes",
                                  "I am using my own calibration data in a different reference frame" = "myown")
                                ),                    
                        bsTooltip('refframe', "CDES must be chosen if using calibration data from Román-Palacios et al. (Petersen et al., 2019)",
                                  placement = "bottom", trigger = "hover",
                                  options = NULL)
                      ),


                    # Misc options
                    checkboxGroupInput('misc', "Miscellaneous options",
                    c("Scale data" = 'scale',
                      "Calculate uncertainties" = 'calcuncertainties')
                    )
                    
                  )
                        ),
                
                # Model selection
                box(width = 4,
                    title = "Step 2: Select Models", solidHeader = TRUE,
                  column(12, "The Bayesian model can take up to three minutes to run",





                          checkboxInput("simulateLM_measured", "Linear model", FALSE),
                          checkboxInput("simulateLM_inverseweights", "Inverse weighted linear model", FALSE),
                          checkboxInput("simulateYork_measured", "York regression", FALSE),
                          checkboxInput("simulateDeming", "Deming regression", FALSE),
                          checkboxInput("simulateBLM_measuredMaterial", "Bayesian simple linear models", FALSE),
                        #  checkboxInput("simulateBLMM_measuredMaterial", "Bayesian simple and mixed linear models", FALSE),
                         
                        #  checkboxInput("fitBayesianMainANCOVASimple", "Bayesian main effects ANCOVA", FALSE),
                        #  checkboxInput("fitBayesianInteractionANCOVASimple", "Bayesian interaction effects ANCOVA", FALSE),

                    bsTooltip('simulateBLM_measuredMaterial', "Running the Bayesian model can take a few minutes. Please be patient.",
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
                    title = "Step 3: Output Options", solidHeader = TRUE,
                  column(12, verbatimTextOutput("lmcal"),
                             verbatimTextOutput("lminversecal"),
                             verbatimTextOutput("york"),
                             verbatimTextOutput("deming"),
                             verbatimTextOutput("blin")
                         ),
                  
                  # Download all calibration data
                  downloadButton("downloadcalibrations", label = "Download calibration output")
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
                    )))
              ),


       #Reconstruction tab
       tabItem(tabName = "reconstruction",
               fluidRow(
                 box(width = 6, 
                     title = "Step 1: Reconstruction setup", solidHeader = TRUE,
                     column(12,
                            radioButtons("calset2", "Pick your calibration",
                                         choices = c("Use calibration data and options from Calibration tab" = 'usecaltab',
                                                     'Use Román-Palacios et al. defaults (Petersen et al., 2019)' = 'petersen2')
                            ),
                            
                            # Download templates
                            downloadButton("BayClump_reconstruction_template.csv", label = "Download reconstruction data template"),
                            
                            # Upload data
                            textInput("path2", "File:"),
                            actionButton("browse2", "Browse", 
                                         icon = icon("search", lib = "font-awesome")
                            ),
                            
                            actionButton("upload3", "Upload data for reconstructions", 
                                         icon = icon("upload", lib = "font-awesome")
                            ),
                            
                            checkboxInput("confirm", "My calibration data and reconstruction data are in the same reference frame")
                     )
                 ),
                 box(width = 6,
                     "More reconstruction options here!!!")
               ),  
              fluidRow(
               box(width = 12,
                    column(12,
                  title = "Reconstruction dummy data", solidHeader = TRUE,
                    plotlyOutput("examplefig2")
                    )
                  )),
              fluidRow(
                box(width = 12,
                    column(12,
                           title = "More dummy data", solidHeader = TRUE,
                           dataTableOutput("table1")
                           ))
              )),

      #User manual
      tabItem(tabName = "usermanual",
              fluidRow(
                column(12,
                  includeMarkdown("userguide.Rmd")))),
    
    #BayWatch
    tabItem(tabName = "demo",
            fluidRow(
              column(12,
                     "Recorded demo here"))),
    
      #Petersen et al cal set
      tabItem(tabName = "petersen",
                column(12,
              fluidRow("Petersen et al. calibration data",
                  dataTableOutput("Petersendat")
                  ))), 
    
      #Citations tab
      tabItem(tabName = "citations",
              fluidRow(
                column(12#,
                  #includeHTML("citations.html")
                  ))),
    
      #Manuscript
      tabItem(tabName = "manuscript",
              fluidRow(
                column(12, "Our manuscript here (if open access)"
                  #includeHTML("Our manuscript here")
                  ))),
    
    #Download
    tabItem(tabName = "download",
            fluidRow(
              column(12,
                     "Package BayClump as an Electron app for desktop and put download links here!")))
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
  
  menuItem("Petersen et al., 2019", tabName = "petersen", 
           icon = icon("table", lib = "font-awesome")
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


