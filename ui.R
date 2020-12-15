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
                                       choices = c('Román-Palacios et al. (Petersen et al., 2019)' = 'petersen',
                                                   "Use my calibration data" = 'mycal',
                                                   "Use both" = "both")
                                       ),
                
                    # Download template
                    downloadButton("BayClump_template.csv", label = "Download calibration data template"),
                    
                    # Upload data
                    textInput("path", "File:"),
                    actionButton("browse", "Browse", 
                                 icon = icon("search", lib = "font-awesome")
                                 ),
                    actionButton("upload", "Upload calibration data", 
                                 icon = icon("upload", lib = "font-awesome")
                                 ),
                    
                    # Uncertainties
                    radioButtons("uncertainties", " ",
                                 c("My data contain uncertainties" = "myuncertainties",
                                   "Use uncertainties from Daëron et al." = "usedaeron")
                                 ),
                    
                    # Digestion temperature
                    div(style="display:inline-block; vertical-align:top;", 
                        selectInput(width = "175px",
                                    "digesttemp", "Digestion temperature",
                                c(" " = "NULL",
                                  "0°C" = "zero",
                                  "25°C" = "twentyfive",
                                  "70°C" = "seventy",
                                  "90°C" = "ninety")
                                )
                      ),
                    div(style="display:inline-block; vertical-align:top; ", 
                        textInput(width = "175px",
                                  "othertemp", "OR input temperature (°C)"
                              )
                        ),
                    
                    # Misc options
                    checkboxInput('scale', 'Scale data'),
                    checkboxInput('calcuncertainties', "Calculate uncertanties")
                    
                  )
                        ),
                
                # Model selection
                box(width = 4,
                    title = "Step 2: Select models", solidHeader = TRUE,
                  column(12,
                         checkboxInput('linear', "Linear model"),
                         checkboxInput('bayes', "Bayesian model"),
                         checkboxInput('york', "York regression"),
                  
                    # Summary stats panel
                    plotOutput(
                    "dummyplot"
                    ),
                    
                    # Run models
                    div(style="display:inline-block; vertical-align:top;", 
                      actionButton('runmodels', "Run selected models", 
                                 icon = icon("cogs", lib = "font-awesome")
                  )
                    ),
                  div(style="display:inline-block; vertical-align:top;", 
                      actionButton('reset', "Reset ALL", 
                                   icon = icon("trash-alt", lib = "font-awesome"))),
                                               bsTooltip('reset', "Warning: This resets EVERYTHING, including data inputs",
                                                         placement = "bottom", trigger = "hover",
                                                         options = NULL)
                      
                  
                )
                ),
                box(width = 4,
                    title = "Step 3: Output options", solidHeader = TRUE,
                  column(12,
                        htmlOutput("dummytext")
                         )
        )
        )
    ),
                
      #Calibration plot tab
      tabItem(tabName = "plots",
              fluidRow(
               box(width = 12,
                    column(12,
                  title = "Calibration plots", solidHeader = TRUE,
                    plotlyOutput("examplefig")
                 # dataTableOutput("caldata")
                  )))),
    
       #Reconstruction tab
       tabItem(tabName = "reconstruction",
               fluidRow(
                 box(width = 6, 
                     title = "Step 1: Calibration Options", solidHeader = TRUE,
                     column(12,
                            radioButtons("calset2", "Select calibration set",
                                         choices = c("Use calibration data and options from Calibration tab" = 'usecaltab',
                                                     'Use Román-Palacios et al. (Petersen et al., 2019)' = 'petersen2',
                                                     "Use my calibration data" = 'mycal2',
                                                     "Use both Román-Palacios et al. (Petersen et al., 2019) and my calibration data" = "both2")
                            ),
                            
                            # Download templates
                            downloadButton("BayClump_cal_template.csv", label = "Download calibration data template"),
                            downloadButton("BayClump_reconstruction_template.csv", label = "Download reconstruction data template"),
                            
                            # Upload data
                            textInput("path2", "File:"),
                            actionButton("browse2", "Browse", 
                                         icon = icon("search", lib = "font-awesome")
                            ),
                            actionButton("upload2", "Upload calibration data", 
                                         icon = icon("upload", lib = "font-awesome")
                            ),
                            
                            actionButton("upload3", "Upload data for reconstructions", 
                                         icon = icon("upload", lib = "font-awesome")
                            ),
                            
                            # Uncertainties
                            radioButtons("uncertainties2", " ",
                                         c("My data contain uncertainties" = "myuncertainties",
                                           "Use uncertainties from Daëron et al." = "usedaeron")
                            ),
                            
                            # Digestion temperature
                            div(style="display:inline-block; vertical-align:top;", 
                                selectInput(width = "175px",
                                            "digesttemp2", "Digestion temperature",
                                            c(" " = "NULL",
                                              "0°C" = "zero",
                                              "25°C" = "twentyfive",
                                              "70°C" = "seventy",
                                              "90°C" = "ninety")
                               )
                            ),
                            div(style="display:inline-block; vertical-align:top; ", 
                               textInput(width = "175px",
                                          "othertemp2", "OR input temperature (°C)"
                                )
                            ),
                            
                            # Misc options
                            checkboxInput('scale2', 'Scale data'),
                            checkboxInput('calcuncertainties2', "Calculate uncertanties")
                            
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
              fluidRow(
                column(12, "Petersen data here"
                 # dataTableOutput("petersenetal")
                  ))),
    
      #Citations tab
      tabItem(tabName = "citations",
              fluidRow(
                column(12,"Citations here"
                  #dataTableOutput("citations")
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
  menuItem("Calibration Setup", tabName = "calibration", 
           icon = icon("drafting-compass", lib = "font-awesome")
           ),
    
  menuItem("Calibration Plots", tabName = "plots", 
           icon =icon("chart-bar", lib = "font-awesome")
           ),

  menuItem("Reconstruction", tabName = "reconstruction", 
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


