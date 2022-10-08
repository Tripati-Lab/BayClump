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
    tabItem(tabName = "bayclump",
            fluidRow( align = "center",
              column(12,
                     h2("Welcome to BayClump!"),
                     br(),
                     h4("BayClump is developed to support and facilitate the use of Bayesian models and analyses involving clumped isotope calibration datasets and
                     associated temperature reconstructions. BayClump, a self-contained a Shiny Dashboard application, facilitates the comparisons of Bayesian 
                     and frequentist models for analyzing clumped isotopes datasets for paleoclimatic reconstructions."
                     ),
                     h4("BayClump is divided into two main sections. First, users can fit linear models using published datasets or their own
                        calibration data. Second, based on the calibration analyses, users can subsequently conduct temperature reconstructions
                        based on parameter distributions derived from the initial calibration step."),
                     h4("We provide additional documentation under the User Manual tab, include details on how to run basic analyses
                        in BayClump under the BayWatch tab, and list additional resources including the source code, citations, 
                        among others, in the remaining tabs."),
                     h4("Finally, note that the code implemented in BayClump is part the an R package bayclumpr. Users interested in
                        analyzing their own data with more flexibility should definitively explore bayclumpr as an option."),
                     br(),
                     h4("Enjoy using the app and please do reach out if you have any questions, concerns, or comments!")
                     ))),
    tabItem(tabName = "calibration",
            fluidRow(
              box(width = 4, 
                  title = h4("Step 1: Calibration Options"), solidHeader = FALSE,
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
                         
                         # Summary stats panel
                         tableOutput("contents"),
                         
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
                         #tags$b("Miscellaneous options"),
                         #checkboxInput("scale" ,'Scale data')
                  )
                  
              ),
              
              # Model selection
              box(width = 4,
                  title = h4("Step 2: Select Models"), solidHeader = FALSE,
                  column(12, #h5("For help choosing an appropriate number of bootstrap replicates or the temperature range for CI estimation, see the User Manual"),
                         numericInput("replication", label = "Number of bootstrap replicates", 
                                      100, min = 2, max = 10000),
                        numericInput("generations", label = "Number of posterior samples to analize", 
                                     3000, min = 100, max = 7000),
                         checkboxInput("cal.ols", "Linear model", FALSE),
                         checkboxInput("cal.wols", "Inverse weighted linear model", FALSE),
                         checkboxInput("cal.york", "York regression", FALSE),
                         checkboxInput("cal.deming", "Deming regression", FALSE),
                         checkboxInput("cal.bayesian", "Bayesian linear models", FALSE),
                         uiOutput('priors'),

                         # Run models
                         div(style="display:inline-block; vertical-align:top;", 
                             actionButton('runmods', "Run selected models", 
                                          icon = icon("cogs", lib = "font-awesome", verify_fa = FALSE)
                             )
                         ),
                    br(),     
                    verbatimTextOutput("modresults")
                         
                  )
              ),
              box(width = 4,
                  title = h4("Step 3: Calibration Output"), solidHeader = FALSE,
                  column(12,
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
                         tags$h4("Bayesian mixed model - no errors"),
                         verbatimTextOutput("blinmwerr", placeholder = TRUE)
                  ),
                  
                  # Download all calibration data
                  downloadButton("downloadcalibrations", label = "Raw results"),
                  downloadButton("downloadBayesian", label = "Bayesian summary"),
                  downloadButton("downloadPosteriorCalibration", label = "Bayesian posteriors")
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
                  column(12, tags$b("Parameters for the selected models are automatically transferred from the Calibration tab to this tab."),
                         tags$br(),
                         # Download templates
                         downloadButton("BayClump_reconstruction_template.csv", label = "Download reconstruction data template"),
                         
                         # Upload data
                         fileInput("reconstructiondata", "Select reconstruction data file", accept = ".csv" 
                         ),
                         
                         # Summary stats panel
                         tableOutput("contents2"),
                         uiOutput("myList2"),
                         tags$b("Models to run:"),
                         checkboxInput("cal.olsRec", "Linear model", FALSE),
                         checkboxInput("cal.wolsRec", "Inverse weighted linear model", FALSE),
                         checkboxInput("cal.yorkRec", "York regression", FALSE),
                         checkboxInput("cal.demingRec", "Deming regression", FALSE),
                         checkboxInput("rec.bayesian", "Bayesian linear models", FALSE),
                         #numericInput("TPriorMean", label = "Mean value for the prior distribution on temperature (10^6/T^2)", 
                          #            11, min = 0, max = 20),
                         #numericInput("TPriorSd", label = "Standard deviation value for the prior distribution on temperature (10^6/T^2)", 
                          #            0.01, min = 0.0001, max = 0.09),
                         bsTooltip('BayesianCalibrationsRec', "This can take several minutes to an hour for large datasets",
                                   placement = "bottom", trigger = "hover",
                                   options = NULL),
                         bsTooltip('simpleInversion', "Uses IPI method from McClelland et al. (2021)",
                                   placement = "bottom", trigger = "hover",
                                   options = NULL),
                         
                         div(style="display:inline-block; vertical-align:top;", 
                             actionButton('runrec', "Run reconstructions", 
                                          icon = icon("cogs", lib = "font-awesome", verify_fa = FALSE)
                             )
                         ),
                         verbatimTextOutput("recresults"),
                         # Download all reconstruction data
                         downloadButton("downloadreconstructions", label = "Download results") #,
                         #downloadButton("downloadreconstructionsPosterior", label = "Download posterior reconstruction output"),
                         #downloadButton("downloadPriorsReconstruction", label = "Download priors")
                         
                  )
              ),
              box(width = 7,
                  title = "Step 2: Temperature reconstructions", solidHeader = TRUE,
                  column(12,
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
    tabItem(tabName = "quick",
            fluidRow(
              column(12,
                     withMathJax(includeMarkdown("doc/quick.Rmd"))
              ))),
    
    tabItem(tabName = "calTab",
            fluidRow(
              column(12,
                     withMathJax(includeMarkdown("doc/calTab.Rmd"))
              ))),
    

    tabItem(tabName = "updatesTab",
            fluidRow(
              column(12,
                     withMathJax(includeMarkdown("doc/updatesTab.Rmd"))
              ))),
    
    tabItem(tabName = "recTab",
            fluidRow(
              column(12,
                     withMathJax(includeMarkdown("doc/recTab.Rmd"))
              ))),
    
    tabItem(tabName = "calPlots",
            fluidRow(
              column(12,
                     withMathJax(includeMarkdown("doc/calPlots.Rmd"))
              ))),
    
    tabItem(tabName = "models",
            fluidRow(
              column(12,
                     withMathJax(includeMarkdown("doc/models.Rmd"))
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
              box(width = 12, title = "Note: This is a preprint and has not yet been accepted for publication. Updates will be posted here.",
                  solidHeader = TRUE,
              column(12, 
                     uiOutput("msframe")
                     )
                   )
            )
            ),
    
    #Download
    tabItem(tabName = "download",
            fluidRow(
              column(12,
                     "BayClump will be made available as a standalone Electron app upon acceptance for publication, and a download link will be provided here")))
  )
)





############ Dashboard Sidebar ############


sidebar <- dashboardSidebar(width = 200,
                            tags$script(JS("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';")),
                            sidebarMenu(
                              menuItem("BayClump", tabName = "bayclump", 
                                       icon = icon("thermometer", lib = "font-awesome", verify_fa = FALSE)
                              ),
                              menuItem("Calibrations", tabName = "calibration", 
                                       icon = icon("drafting-compass", lib = "font-awesome", verify_fa = FALSE)
                              ),
                              
                              menuItem("Calibration Plots", tabName = "plots", 
                                       icon =icon("chart-bar", lib = "font-awesome", verify_fa = FALSE)
                              ),
                              

                              menuItem("Reconstructions", tabName = "reconstruction", 
                                       icon = icon("chart-area", lib = "font-awesome", verify_fa = FALSE)
                              ),
                              
                              menuItem("User Manual", tabName = "usermanual", 
                                       icon = icon("question-circle", lib = "font-awesome", verify_fa = FALSE),
                                       menuSubItem('Quick start',
                                                   tabName = 'quick'),
                                       menuSubItem('Calibrations',
                                                   tabName = 'calTab'),
                                       menuSubItem('Calibration plots',
                                                   tabName = 'calPlots'),
                                       menuSubItem('Reconstructions',
                                                   tabName = 'recTab'),
                                       menuSubItem('Models',
                                                   tabName = 'models'),
                                       menuSubItem('Updates!',
                                                   tabName = 'updatesTab')
                              ),
                              
                              menuItem("Citations", tabName = "citations", 
                                       icon = icon("align-right", lib = "font-awesome", verify_fa = FALSE)
                              ),
                              
                              menuItem("Manuscript", tabName = "manuscript", 
                                       icon = icon("scroll", lib = "font-awesome", verify_fa = FALSE)
                              ),
                              
                              menuItem("Download BayClump", tabName = "download", 
                                       icon = icon("laptop-code", lib = "font-awesome", verify_fa = FALSE)
                              ),
                              
                              menuItem("Contacts", icon = icon("address-card", lib = "font-awesome"), verify_fa = FALSE,
                                       a(actionButton(inputId = "email1", label = "Email", 
                                                      icon = icon("envelope", lib = "font-awesome", verify_fa = FALSE)),
                                         href="mailto:hannah.carroll@paleoeco.org"),
                                       a(actionButton(inputId = "github", label = "GitHub", 
                                                      icon = icon("github", lib = "font-awesome", verify_fa = FALSE),
                                                      onclick ="window.open('https://github.com/Tripati-Lab/BayClump/', '_blank')")),
                                       a(actionButton(inputId = "bugs", label = "Submit a Bug Report", 
                                                      icon = icon("bug", lib = "font-awesome", verify_fa = FALSE),
                                                      onclick ="window.open('https://github.com/Tripati-Lab/BayClump/issues', '_blank')")
                                       )
                              )
                            )
)


########### Dashboard page ####################
dashboardPage(
  skin = "black",
  header,
  sidebar,
  body
  ,tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
  )
)


