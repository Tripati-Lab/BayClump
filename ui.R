# Define UI for application

header <-  dashboardHeader(#title = "BayClump: Bayesian methods for clumped isotope paleothermometry",
                           titleWidth = 750,    title = shinyDashboardLogoDIY(
                             boldText = "BayClump"
                             ,mainText = ""
                             ,textSize = 22
                             ,badgeText = "v1.0"
                             ,badgeTextColor = "black"
                             ,badgeTextSize = 2
                             ,badgeBackColor = "white"
                             ,badgeBorderRadius = 3
                           )
)




body <- dashboardBody(
  tags$head(tags$style(HTML('
      .main-sidebar { font-size: 15px; }
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
                     h1("Welcome to BayClump!"),
                     h5("Bayesian methods for clumped isotope paleothermometry"),
                     br(),
                     
                     tags$img(
                       src = "bayclump.png",
                       alt = "Success, you are correct",
                       width = 200
                       #, height = 25
                     ), 
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
            fluidRow( align = "center",
                      column(12,
                             h1("Calibrations"),
                             h5("Please select a combination of calibration dataset(s) and regressions models. 
                                Information in this tab will be passed on to both the calibrations plot and reconstructions tab.
                                Feel free to either select any of the pre-loaded datasets or upload your own. 
                                If you decide to use your local dataset for the analyses, please make sure to follow
                                the structure outlined in the template.")
                      )),
            fluidRow(
              box(width = 4, 
                  title = h4(HTML("<b>Step 1: Calibration Options</b>"), align='center'), solidHeader = FALSE,
                  column(12,
                         radioButtons("calset", "Select calibration set",
                                      choices = c('Román-Palacios et al. - Model 1' = 'model1',
                                                  'Román-Palacios et al. - Model 2' = 'model2',
                                                  "Use Román-Palacios et al. - Models 1 & 2" = 'model1and2',
                                                  "Upload my own calibration data" = 'mycal',
                                                  "Use all" = "all")
                         ),
                         
                         # Summary stats panel
                         column(12, align="center", tableOutput("contents"),
                         
                         # Download template
                         downloadButton("BayClump_cal_temp", label = "Download template")
),

                         # Reference frame
                         # div(style="display:inline-block; vertical-align:top;", 
                         #     radioButtons("refframe", "Reference frame",
                         #                  c("Use the Intercarb carbon dioxide equilibrium scale at 90°C (I-CDES 90)" = "icdes90",
                         #                    "I am using my own calibration data in a different reference frame" = "myown")
                         #     ),                    
                         #     bsTooltip('refframe', "I-CDES 90 must be chosen if using calibration data from Román-Palacios et al.",
                         #               placement = "bottom", trigger = "hover",
                         #               options = NULL)
                         # ),
                         # 
                         # Misc options
                         #tags$b("Miscellaneous options"),
                         #checkboxInput("scale" ,'Scale data')
                  )
                  
              ),
              
              # Model selection
              box(width = 4,
                  title = h4(HTML("<b>Step 2: Select Models</b>"), align='center'), solidHeader = FALSE,
                  column(12, #h5("For help choosing an appropriate number of bootstrap replicates or the temperature range for CI estimation, see the User Manual"),
                         numericInput("replication", label = "Number of bootstrap replicates", 
                                      100, min = 2, max = 10000),
                        numericInput("generations", label = "Number of posterior samples", 
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
                  title = h4(HTML("<b>Step 3: Calibration Output</b>"), align='center'), solidHeader = FALSE,
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
                  column(12, align="center",
                         downloadButton("downloadcalibrations", label = "Raw results"),
                         downloadButton("downloadBayesian", label = "Bayesian summary"),
                         downloadButton("downloadPosteriorCalibration", label = "Bayesian posteriors")
                  )
              )
            )
    ),
    
    #Calibration plot tab
    tabItem(tabName = "plots",
            fluidRow( align = "center",
                      column(12,
                             h1("Calibration plots"),
                             h5("BayClump allows for a quick visualization of the calibrations performed
                                within the app. Interactive visualizations, conducted using plotly, can be downloaded to
                                still images.")
                      )),
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
            fluidRow( align = "center",
                      column(12,
                             h1("Reconstructions tab"),
                             h5("This section of BayClump focuses on deriving teperature reconstructions
                                from existing calibrations. Note that  parameters for the selected models are automatically transferred from the Calibration tab.
                                Please follow the structure of the template file provided in BayClump to upload your own data.")
                      )),
            fluidRow(
              box(width = 5, 
                  title = h4(HTML("<b>Step 1: Reconstruction setup</b>") , align='center') ,solidHeader = TRUE,
                  column(12, #("Parameters for the selected models are automatically transferred from the Calibration tab."),
                         #tags$br(),
                         column(12, align="center",
                         # Download templates
                         downloadButton("BayClump_reconstruction_template.csv", label = "Download template"),
                         
                         # Upload data
                         fileInput("reconstructiondata", "Select reconstruction data file", accept = ".csv" 
                         )),
                         
                         # Summary stats panel
                         tableOutput("contents2"),
                         uiOutput("myList2"),
                         tags$b("Models to run:"),
                         checkboxInput("cal.olsRec", "Linear model", FALSE),
                         checkboxInput("cal.wolsRec", "Inverse weighted linear model", FALSE),
                         checkboxInput("cal.yorkRec", "York regression", FALSE),
                         checkboxInput("cal.demingRec", "Deming regression", FALSE),
                         checkboxInput("rec.bayesian", "Bayesian linear models", FALSE),
                         
                         column(12, align="center",
                         div(style="display:inline-block; vertical-align:top;", 
                             actionButton('runrec', "Run reconstructions", 
                                          icon = icon("cogs", lib = "font-awesome", verify_fa = FALSE)
                             )
                         ),
                         verbatimTextOutput("recresults")) #,
                         #downloadButton("downloadreconstructionsPosterior", label = "Download posterior reconstruction output"),
                         #downloadButton("downloadPriorsReconstruction", label = "Download priors")
                         
                  )
              ),
              box(width = 7,
                  title = h4(HTML("<b>Step 2: Temperature reconstructions</b>"), align = 'center'), solidHeader = TRUE,
                  column(12,
                         tableOutput("lmrecswun"),
                         tableOutput("lminverserecswun"),
                         tableOutput("yorkrecswun"),
                         tableOutput("demingrecswun"),
                         tableOutput("Bpredictions"),
                         tableOutput("BpredictionsErrors"),
                         tableOutput("BpredictionsBLMM"),
                         # Download all reconstruction data
                         column(12, align="center",
                         downloadButton("downloadreconstructions", label = "Download results")
                         )
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
            
            fluidRow( align = "center",
                      column(12,
                             h1("Baywatch"),
                             h5("A recorded tutorial of BayClump.")
                      ))),
    
    #Citations tab
    tabItem(tabName = "citations",
            fluidRow(
              column(12,
                     DT::dataTableOutput("bibTable")
              ))),
    
    #Manuscript
    tabItem(tabName = "manuscript", 
            fluidRow( align = "center",
                      column(12,
                             h1("Manuscript tab"),
                             h5("Below we provide a copy of the latest version of the manuscript associated with BayClump. 
                                Note: This is a preprint and has not yet been accepted for publication. Updates will be posted here.")
                      )),
            fluidRow(
              uiOutput("msframe")
            )
            ),
    
    #Download
    tabItem(tabName = "download",
            fluidRow( align = "center",
                      column(12,
                             h1("Download BayClump"),
                             h5("BayClump will be made available as a standalone 
                                Electron app upon acceptance for publication, and a download link 
                                will be provided here")
                      )))
  )
)





############ Dashboard Sidebar ############


sidebar <- dashboardSidebar(width = 200,
                            tags$script(JS("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';")),
                            sidebarMenu(
                              menuItem("BayClump", tabName = "bayclump", 
                                       icon = icon("thermometer", lib = "font-awesome", verify_fa = FALSE)
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
                              menuItem("Baywatch", tabName = "demo", 
                                       icon = icon("glasses", lib = "font-awesome", verify_fa = FALSE)
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
                              ),
                              menuItem("Citations", tabName = "citations", 
                                       icon = icon("align-right", lib = "font-awesome", verify_fa = FALSE)
                              )
                            )
)


theme <- create_theme(
  theme = "journal",
  adminlte_color(
    green = "#3fff2d",
    blue = "#2635ff",
    red = " #ff2b2b",
    yellow = "#feff6e",
    fuchsia = "#ff5bf8",
    navy = "#374c92",
    purple = "#615cbf",
    maroon = "#b659c9",
    light_blue = "#F78A54"
  ),
  adminlte_sidebar(
    dark_bg = "#FAE6C2",
    dark_hover_bg = "#F78A54",
    dark_color = "#2E3440",
    dark_submenu_color = "#2E3440",
    light_submenu_hover_color = "#81A1C1",
    light_submenu_bg = "#81A1C1"
  ),
  adminlte_global(
    content_bg = "white",
    box_bg = "#b2d4d9", 
    info_box_bg = "#93B8B1"
  )
)


########### Dashboard page ####################
dashboardPage(
  freshTheme = theme,
  #skin = "blue-light",
  header,
  sidebar,
  body
)


