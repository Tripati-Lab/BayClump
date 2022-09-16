# Define server logic
server <- function(input, output, session) { 
  options(shiny.maxRequestSize=800*1024^2) 
  
  # Show package citations
  get_path <- reactive({
    path <- file.path(paste0(getwd()), paste("Rpackages", ".bib", sep=""))
    return(path)
  })
  
  get_bib <- reactive({
    ## insert your get bib logic here
    pkgbib <- bibtex::read.bib("Rpackages.bib")
    
    df <- bib2df(get_path()) %>% dplyr::select(BIBTEXKEY, NOTE, AUTHOR, TITLE, YEAR, JOURNAL, VOLUME, PAGES, URL)
    df$AUTHOR <- unlist(lapply(df$AUTHOR, paste, collapse = ", "))
    df <- df %>% arrange(AUTHOR)
    return(df)
  })
  
  output$bibTable <- DT::renderDataTable({
    
    bib <- get_bib()
    bib
    
  }, caption = "This is an automatically generated list of R packages used to render and run BayClump. Citation information is provided by package authors.", options = list(pageLength = 20, info = FALSE))
  
  # Calibration tab
  
  output$BayClump_cal_temp <- downloadHandler(
    filename = "BayClump_calibration_template.csv",
    content = function(file) {
      write.csv(BayClump_calibration_template, file, row.names = FALSE)
    }
  )
  
  
  
  calibrationData = reactive({
    switch(input$calset,
           "model1" = return(Petersen),
           "model2" = return(Anderson),
           "model1and2" = return(PetersenAnderson),
           "mycal" = reactiveValues({
             req(input$calibrationdata)
             n_rows = length(count.fields(input$calibrationdata$datapath))
             df_out = read.csv(input$calibrationdata$datapath)
             return(df_out)
           }),
           "all" = reactiveValues({
             req(input$calibrationdata)
             n_rows = length(count.fields(input$calibrationdata$datapath))
             df_out = read.csv(input$calibrationdata$datapath)
             alldat <- rbind(df_out, PetersenAnderson)
             return(alldat)
           })
    )
  })

  #For parameter estimates
  if(exists("wb")) rm(wb) # Delete any existing workbook in preparation for new results
  wb <- createWorkbook("calibration output") # Prepare a workbook for calibration outputs
  
  #For convergence
  if(exists("wb3")) rm(wb3) # Delete any existing workbook in preparation for new results
  wb3 <- createWorkbook("Bayesian output") # Prepare a workbook for calibration outputs
  
  
  if(exists("wb4")) rm(wb4) # Delete any existing workbook in preparation for new results
  wb4 <- createWorkbook("Bayesian posterior output") # Prepare a workbook for calibration outputs
  
  #For convergence
  if(exists("wb5")) rm(wb5) # Delete any existing workbook in preparation for new results
  wb5 <- createWorkbook("Bayesian reconstruction posterior output") # Prepare a workbook for calibration outputs
  
  
  observeEvent(calibrationData(),{
  output$myList <-  renderUI({
    numericInput("samples", min = 3, max = nrow(calibrationData()), 
                label = paste0("Number of observations per bootstrap sample (max. recommended: ", nrow(calibrationData()) ,")" ), 
                value =  nrow(calibrationData()))
  })
  })
  

  observe({
    output$contents <- renderTable({
      calsummary <- calibrationData() %>%
        summarize(
          "Total samples" = length(calibrationData()$Sample.Name),
          "Unique samples" = length(unique(calibrationData()$Sample.Name)),
          "Total replicates" = sum(calibrationData()$N),
          "Materials" = length(unique(calibrationData()$Material))
        )
      return(calsummary)
    }, 
    rownames=FALSE, options = list(pageLength = 1, info = FALSE)
    )
  }) 
  

  toListen <- reactive({
    list(input$runmods)
  })
  
  modresult <- eventReactive(toListen() , {
    
    if ("all" %in% input$calibrationdata) {
      print(noquote("Please upload calibration data first"))
    } 
    

    priors <<- input$priors
    replicates <<- input$replication
    ngenerationsBayes <<- input$generations
    ngenerationsBayesianPredictions <<- input$generations

    ##Download priors (calibration; not reactive...)
    if(exists("wb6")) rm(wb6)
    wb6 <- createWorkbook("Priors - Calibration")
    if("Settings" %in% names(wb6) == TRUE) 
    {removeWorksheet(wb6, "Settings")}
    if("Distributions" %in% names(wb6) == TRUE) 
    {removeWorksheet(wb6, "Distributions") }
    pd <- generatePriorDistCalibration(prior = priors)
    addWorksheet(wb6, "Settings") # Add a blank sheet
    writeData(wb6, sheet = "Settings", attr(pd, "params")) # Write regression data
    addWorksheet(wb6, "Distributions") # Add a blank sheet
    writeData(wb6, sheet = "Distributions", pd) # Write regression data
    
    output$downloadPriorsCalibration <- downloadHandler(
      filename = function() { 
        paste("Priors_calibration_", Sys.time(), ".xlsx", sep="")
      },
      content = function(file) {
        saveWorkbook(wb6, file, overwrite = TRUE)
      }
    )
    
    ##Download priors (reconstructions)
    if(exists("wb7")) rm(wb7)
    wb7 <- createWorkbook("Priors - Reconstruction")
    if("Settings" %in% names(wb7) == TRUE) 
    {removeWorksheet(wb7, "Settings")}
    if("Distributions" %in% names(wb7) == TRUE) 
    {removeWorksheet(wb7, "Distributions") }
    pd <- generatePriorReconstructions(prior = priors)
    addWorksheet(wb7, "Settings") # Add a blank sheet
    writeData(wb7, sheet = "Settings", attr(pd, "params")) # Write regression data
    addWorksheet(wb7, "Distributions") # Add a blank sheet
    writeData(wb7, sheet = "Distributions", pd) # Write regression data
    
    output$downloadPriorsReconstruction <- downloadHandler(
      filename = function() { 
        paste("Priors_reconstruction_", Sys.time(), ".xlsx", sep="")
      },
      content = function(file) {
        saveWorkbook(wb7, file, overwrite = TRUE)
      }
    )
    
    
    
    
    # Remove existing worksheets from wb on "run" click, if any
    if("Linear regression" %in% names(wb) == TRUE) 
    {removeWorksheet(wb, "Linear regression") & removeWorksheet(wb, "Linear regression CI")}
    if("Inverse linear regression" %in% names(wb) == TRUE) 
    {removeWorksheet(wb, "Inverse linear regression") & removeWorksheet(wb, "Inverse linear regression CI")}
    if("York regression" %in% names(wb) == TRUE) 
    {removeWorksheet(wb, "York regression") & removeWorksheet(wb, "York regression CI")}
    if("Deming regression" %in% names(wb) == TRUE) 
    {removeWorksheet(wb, "Deming regression") & removeWorksheet(wb, "Deming regression CI")}
    if("Bayesian model no errors" %in% names(wb) == TRUE) 
    {removeWorksheet(wb, "Bayesian model no errors") & removeWorksheet(wb, "Bayesian model no errors CI")}
    if("Bayesian model with errors" %in% names(wb) == TRUE) 
    {removeWorksheet(wb, "Bayesian model with errors") & removeWorksheet(wb, "Bayesian model with errors CI")}
    if("Bayesian mixed w errors" %in% names(wb) == TRUE)
    {removeWorksheet(wb, "Bayesian mixed w errors") & removeWorksheet(wb, "Bayesian mixed w errors CI")}
    
    
    ##Also for the Bayesian sheet
    if("Bayesian model no errors" %in% names(wb3) == TRUE) 
    {removeWorksheet(wb3, "Bayesian model no errors")}
    if("Bayesian model with errors" %in% names(wb3) == TRUE) 
    {removeWorksheet(wb3, "Bayesian model with errors") }
    if("Bayesian mixed w errors" %in% names(wb3) == TRUE)
    {removeWorksheet(wb3, "Bayesian mixed w errors")}
    
    ##Also for the Bayesian posterior sheet
    if("Bayesian model no errors" %in% names(wb4) == TRUE) 
    {removeWorksheet(wb4, "Bayesian model no errors")}
    if("Bayesian model with errors" %in% names(wb4) == TRUE) 
    {removeWorksheet(wb4, "Bayesian model with errors") }
    if("Bayesian mixed w errors" %in% names(wb4) == TRUE)
    {removeWorksheet(wb4, "Bayesian mixed w errors")}
    

    
    lmcals <<- NULL
    lminversecals <<- NULL
    yorkcals <<- NULL
    demingcals <<- NULL
    bayeslincals <<- NULL

    
    # Calibration data
    
    calData <<- NULL
    calData <<- calibrationData()

    samples <<- if(is.null(input$samples)){nrow(calData) }else{input$samples}
    samples <<- ifelse(samples==2, 3, samples)

    # Recode NA or 0 error values to dummy value
    calData$D47error[calData$D47error == 0] <<- 0.000001
    calData$TempError[calData$TempError == 0] <<- 0.000001
    calData$D47error[is.na(calData$D47error)] <<- 0.000001
    calData$TempError[is.na(calData$TempError)] <<- 0.000001
    calData$Material[is.na(calData$Material)] <<- 1
    
    
    
    calData$Material <<- factor(calData$Material, labels = seq(1:length(unique(calData$Material))))
    if(min(as.numeric(calData$Material)) != 1){
    print(noquote("The sequence of Materials now start from 1"))
    }
    
    
    ##Limits of the CI
    
    minLim <- ifelse(input$range[1]==0, min(calData$Temperature),input$range[1])
    maxLim <- ifelse(input$range[2]==0, max(calData$Temperature),input$range[2])
    

    
    if(input$scale == TRUE) {
      calData$Temperature <<- scale(calData$Temperature)
      calData$TempError <<- scale(calData$TempError)
      calData$D47 <<- scale(calData$D47)
      calData$D47error <<- scale(calData$D47error)
    }
    
    NegErrors <- any(calData$D47error <= 0) | any(calData$TempError <= 0)
    if(NegErrors) {
      print(noquote("Invalid input: 0 or negative uncertainty values"))
      }

    
    if(input$simulateLM_measured == FALSE &
       input$simulateLM_inverseweights == FALSE &
       input$simulateYork_measured == FALSE &
       input$simulateDeming == FALSE &
       input$BayesianCalibrations == FALSE) {
      print(noquote("Please select at least one model"))
      }
    
   if(NegErrors == FALSE){
    
    if(input$simulateLM_measured != FALSE |
       input$simulateLM_inverseweights != FALSE |
       input$simulateYork_measured != FALSE |
       input$simulateDeming != FALSE |
       input$BayesianCalibrations != FALSE) {
      
      withProgress(message = "Running selected models, please wait", {
        
        if(input$simulateLM_measured == FALSE) {
        }
        if(input$simulateLM_inverseweights == FALSE) {
        }
        if(input$simulateYork_measured == FALSE) {
        }
        if(input$simulateDeming == FALSE) {
        }
        if(input$BayesianCalibrations == FALSE) {
        }
        
        totalModels <- c(input$simulateLM_measured, input$simulateLM_inverseweights, input$simulateYork_measured,
                         input$simulateDeming, input$BayesianCalibrations)
        TotProgress <- length(which(totalModels==T))
        
        if(input$simulateLM_measured != FALSE) {
          sink(file = "out/linmodtext.txt", type = "output")
          lmcals <<- simulateLM_measured(calData, replicates = replicates, samples = samples)
          sink()
          incProgress(1/TotProgress, detail="...Done fitting the OLS...")
          
          lmci <<- RegressionSingleCI(data = lmcals, from = minLim, to = maxLim)
          lmcalci <- as.data.frame(lmci)
          
          output$lmcalibration <- renderPlotly({
            lmfig <- plot_ly(calibrationData()
            )
            lmfig <- lmfig %>%
              add_trace(x = ~calibrationData()$Temperature, 
                        y = ~D47,
                        type = "scatter", 
                        mode = "markers", 
                        marker = list(color = "black"),
                        opacity = 0.5,
                        name = "Raw data",
                        text = as.character(calibrationData()$Sample.Name),
                        hovertemplate = paste(
                          "<b>Sample: %{text}</b><br><br>",
                          "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                          "Δ<sub>47</sub> (‰): %{y}<br>",
                          "Mineralogy: ", as.character(calibrationData()$Mineralogy),"<br>",
                          "Type: ", as.character(calibrationData()$Material),
                          "<extra></extra>"))
            lmfig <- lmfig %>% 
              add_ribbons(data = lmcalci,
                          x = ~x,
                          y = ~median_est,
                          ymin = ~ci_lower_est,
                          ymax = ~ci_upper_est,
                          line = list(color = "#ffd166"),
                          fillcolor = "#ffd166",
                          opacity = 0.5,
                          name = "95% CI",
                          hovertemplate = paste(
                            "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                            "Δ<sub>47</sub> (‰): %{y}<br>")) %>%
              add_lines(data = lmcalci,
                        x = ~x,
                        y = ~median_est,
                        name = "Mean estimate",
                        line = list(color = "black", dash = "dash"),
                        hovertemplate = paste(
                          "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                          "Δ<sub>47</sub> (‰): %{y}<br>"))
            lmfig <- lmfig %>% layout(title = "<b> Linear calibration model </b>",
                                      legend=list(title=list(text="Legend")),
                                      xaxis = list(title = "Temperature (10<sup>6</sup>/T<sup>2</sup>)", hoverformat = ".1f"), 
                                      yaxis = list(title = "Δ<sub>47</sub> (‰)", hoverformat = ".3f"))
            
            return(lmfig)
          })
          
          addWorksheet(wb, "Linear regression") # Add a blank sheet
          addWorksheet(wb, "Linear regression CI") # Add a blank sheet 
          
          lmcalci2 <- lmcalci
          names(lmcalci2) <- c("10^6/T^2", "D47_median_est",	"D47_ci_lower_est",	"D47_ci_upper_est")
          
          writeData(wb, sheet = "Linear regression", lmcals) # Write regression data
          writeData(wb, sheet = "Linear regression CI", lmcalci2)
          
          
          
          # cat("Linear regression complete \n *R^2=", round(unlist(attributes(lmcals)$R2[1],4)),
          #     " (95% CI, ",round(unlist(attributes(lmcals)$R2[2],4)),"-",round(unlist(attributes(lmcals)$R2[3],4)),")"
          # )
          

          output$lmcal <- renderPrint({
            do.call(rbind.data.frame,apply(lmcals, 2, function(x){
              cbind.data.frame(Mean= round(mean(x), 4), `SE`=round(sd(x)/sqrt(length(x)),7))
            }))
          })
          
        }
        
        if(input$simulateLM_inverseweights != FALSE) {
          sink(file = "out/inverselinmodtext.txt", type = "output")
          lminversecals <<- simulateLM_inverseweights(calData, replicates = replicates, samples = samples)
          sink()
          incProgress(1/TotProgress, detail="...Done fitting weighted OLS...")
          
          lminverseci <- RegressionSingleCI(data = lminversecals, from = minLim, to = maxLim)
          lminversecalci <- as.data.frame(lminverseci)
          
          output$lminversecalibration <- renderPlotly({
            lminversefig <- plot_ly(data = calibrationData()
            )
            lminversefig <- lminversefig %>% 
              add_trace(x = ~calibrationData()$Temperature, 
                        y = ~D47,
                        type = "scatter", 
                        mode = "markers", 
                        marker = list(color = "black"),
                        opacity = 0.5,
                        name = "Raw data",
                        text = as.character(calibrationData()$Sample.Name),
                        hovertemplate = paste(
                          "<b>Sample: %{text}</b><br><br>",
                          "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                          "Δ<sub>47</sub> (‰): %{y}<br>",
                          "Mineralogy: ", as.character(calibrationData()$Mineralogy),"<br>",
                          "Type: ", as.character(calibrationData()$Material),
                          "<extra></extra>")) %>%
              add_ribbons(data = lminversecalci,
                          x = ~x,
                          y = ~median_est,
                          ymin = ~ci_lower_est,
                          ymax = ~ci_upper_est,
                          line = list(color = "#ffd166"),
                          fillcolor = "#ffd166",
                          opacity = 0.5,
                          name = "95% CI",
                          hovertemplate = paste(
                            "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                            "Δ<sub>47</sub> (‰): %{y}<br>")) %>%
              add_lines(data = lminversecalci,
                        x = ~x,
                        y = ~median_est,
                        name = "Mean estimate",
                        line = list(color = "black", dash = "dash"),
                        hovertemplate = paste(
                          "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                          "Δ<sub>47</sub> (‰): %{y}<br>"))
            lminversefig <- lminversefig %>% layout(title = "<b> Inverse linear calibration model </b>",
                                                    legend=list(title=list(text="Legend")),
                                                    xaxis = list(title = "Temperature (10<sup>6</sup>/T<sup>2</sup>)", hoverformat = ".1f"), 
                                                    yaxis = list(title = "Δ<sub>47</sub> (‰)", hoverformat = ".3f"))
            
            return(lminversefig)
          })
          
          addWorksheet(wb, "Inverse linear regression") # Add a blank sheet
          addWorksheet(wb, "Inverse linear regression CI") # Add a blank sheet 
          
          lminversecalci2 <- lminversecalci
          names(lminversecalci2) <- c("10^6/T^2", "D47_median_est",	"D47_ci_lower_est",	"D47_ci_upper_est")
          
          writeData(wb, sheet = "Inverse linear regression", lminversecals) # Write regression data
          writeData(wb, sheet = "Inverse linear regression CI", lminversecalci2)
          
          # cat("\nInverse regression complete \n *R^2=", round(unlist(attributes(lminversecals)$R2[1],4)),
          #     " (95% CI, ",round(unlist(attributes(lminversecals)$R2[2],4)),"-",round(unlist(attributes(lminversecals)$R2[3],4)),")"
          # )
          
          output$lminversecal <- renderPrint({
            do.call(rbind.data.frame,apply(lminversecals, 2, function(x){
              cbind.data.frame(Mean= round(mean(x), 4),`SE`=round(sd(x)/sqrt(length(x)),7))
            }))
          })
          
        }
        
        if(input$simulateYork_measured != FALSE) {
          sink(file = "out/yorkmodtext.txt", type = "output")
          yorkcals <<- simulateYork_measured(calData, replicates = replicates, samples = samples)
          sink()
          incProgress(1/TotProgress, detail="...Done fitting York regression...")
          
          yorkci <- RegressionSingleCI(data = yorkcals, from = minLim, to = maxLim)
          yorkcalci <- as.data.frame(yorkci)
          
          output$yorkcalibration <- renderPlotly({
            yorkfig <- plot_ly(data = calibrationData()
            )
            yorkfig <- yorkfig %>% 
              add_trace(x = ~calibrationData()$Temperature, 
                        y = ~D47,
                        type = "scatter", 
                        mode = "markers", 
                        marker = list(color = "black"),
                        opacity = 0.5,
                        name = "Raw data",
                        text = as.character(calibrationData()$Sample.Name),
                        hovertemplate = paste(
                          "<b>Sample: %{text}</b><br><br>",
                          "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                          "Δ<sub>47</sub> (‰): %{y}<br>",
                          "Mineralogy: ", as.character(calibrationData()$Mineralogy),"<br>",
                          "Type: ", as.character(calibrationData()$Material),
                          "<extra></extra>")) %>%
              add_ribbons(data = yorkcalci,
                          x = ~x,
                          y = ~median_est,
                          ymin = ~ci_lower_est,
                          ymax = ~ci_upper_est,
                          line = list(color = "#ffd166"),
                          fillcolor = "#ffd166",
                          opacity = 0.5,
                          name = "95% CI",
                          hovertemplate = paste(
                            "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                            "Δ<sub>47</sub> (‰): %{y}<br>")) %>%
              add_lines(data = yorkcalci,
                        x = ~x,
                        y = ~median_est,
                        name = "Mean estimate",
                        line = list(color = "black", dash = "dash"),
                        hovertemplate = paste(
                          "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                          "Δ<sub>47</sub> (‰): %{y}<br>"))
            yorkfig <- yorkfig %>% layout(title = "<b> York calibration model </b>",
                                          legend=list(title=list(text="Legend")),
                                          xaxis = list(title = "Temperature (10<sup>6</sup>/T<sup>2</sup>)", hoverformat = ".1f"), 
                                          yaxis = list(title = "Δ<sub>47</sub> (‰)", hoverformat = ".3f"))
            
            return(yorkfig)
          })
          
          addWorksheet(wb, "York regression") # Add a blank sheet
          addWorksheet(wb, "York regression CI") # Add a blank sheet 
          
          yorkcalci2 <- yorkcalci
          names(yorkcalci2) <- c("10^6/T^2", "D47_median_est",	"D47_ci_lower_est",	"D47_ci_upper_est")
          
          writeData(wb, sheet = "York regression", yorkcals) # Write regression data
          writeData(wb, sheet = "York regression CI", yorkcalci2)
          
          # cat("\nYork regression complete")
          
          output$york <- renderPrint({
            do.call(rbind.data.frame,apply(yorkcals, 2, function(x){
              cbind.data.frame(Mean= round(mean(x), 4),`SE`=round(sd(x)/sqrt(length(x)),7))
            }))
          })
          
        }
        
        if(input$simulateDeming != FALSE) {
          sink(file = "out/demingmodtext.txt", type = "output")
          demingcals <<- simulateDeming(calData, replicates = replicates, samples = samples)
          sink()
          incProgress(1/TotProgress, detail="...Done fitting Deming regression model...")
          
          demingci <- RegressionSingleCI(data = demingcals, from = minLim, to = maxLim)
          demingcalci <- as.data.frame(demingci)
          
          output$demingcalibration <- renderPlotly({
            demingfig <- plot_ly(data = calibrationData()
            )
            demingfig <- demingfig %>% 
              add_trace(x = ~calibrationData()$Temperature, 
                        y = ~D47,
                        type = "scatter", 
                        mode = "markers", 
                        marker = list(color = "black"),
                        opacity = 0.5,
                        name = "Raw data",
                        text = as.character(calibrationData()$Sample.Name),
                        hovertemplate = paste(
                          "<b>Sample: %{text}</b><br><br>",
                          "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                          "Δ<sub>47</sub> (‰): %{y}<br>",
                          "Mineralogy: ", as.character(calibrationData()$Mineralogy),"<br>",
                          "Type: ", as.character(calibrationData()$Material),
                          "<extra></extra>")) %>%
              add_ribbons(data = demingcalci,
                          x = ~x,
                          y = ~median_est,
                          ymin = ~ci_lower_est,
                          ymax = ~ci_upper_est,
                          line = list(color = "#ffd166"),
                          fillcolor = "#ffd166",
                          opacity = 0.5,
                          name = "95% CI",
                          hovertemplate = paste(
                            "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                            "Δ<sub>47</sub> (‰): %{y}<br>")) %>%
              add_lines(data = demingcalci,
                        x = ~x,
                        y = ~median_est,
                        name = "Mean estimate",
                        line = list(color = "black", dash = "dash"),
                        hovertemplate = paste(
                          "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                          "Δ<sub>47</sub> (‰): %{y}<br>"))
            demingfig <- demingfig %>% layout(title = "<b> Deming calibration model </b>",
                                              legend=list(title=list(text="Legend")),
                                              xaxis = list(title = "Temperature (10<sup>6</sup>/T<sup>2</sup>)", hoverformat = ".1f"), 
                                              yaxis = list(title = "Δ<sub>47</sub> (‰)", hoverformat = ".3f"))
            
            return(demingfig)
          })
          
          addWorksheet(wb, "Deming regression") # Add a blank sheet
          addWorksheet(wb, "Deming regression CI") # Add a blank sheet 
          
          demingcalci2 <- demingcalci
          names(demingcalci2) <- c("10^6/T^2", "D47_median_est",	"D47_ci_lower_est",	"D47_ci_upper_est")
          
          writeData(wb, sheet = "Deming regression", demingcals) # Write regression data
          writeData(wb, sheet = "Deming regression CI", demingcalci2)
          
          #cat("\nDeming regression complete")
          output$deming <- renderPrint({
            do.call(rbind.data.frame,apply(demingcals, 2, function(x){
              cbind.data.frame(Mean= round(mean(x), 4),`SE`=round(sd(x)/sqrt(length(x)),7))
            }))
          })
          
        }
        
        if(input$BayesianCalibrations != FALSE) {
          
          sink(file = "out/Bayeslinmodtext.txt", type = "output")
          bayeslincals <<- fitClumpedRegressions(calibrationData = calData, 
                                                 priors = priors,
                                                 numSavedSteps = ngenerationsBayes,
                                                 samples = samples)

          PostBLM1_fit_NoErrors <-do.call(rbind, mcmc.list(
            lapply(1:ncol(bayeslincals$BLM1_fit_NoErrors), function(x) {
            mcmc(as.array(bayeslincals$BLM1_fit_NoErrors)[,x,])
            })))
          
          PostBLM1_fit <- do.call(rbind, mcmc.list(
            lapply(1:ncol(bayeslincals$BLM1_fit), function(x) {
              mcmc(as.array(bayeslincals$BLM1_fit)[,x,])
            })))
          
          PostBLM3_fit <- do.call(rbind, mcmc.list(
            lapply(1:ncol(bayeslincals$BLM3_fit), function(x) {
              mcmc(as.array(bayeslincals$BLM3_fit)[,x,])
            })))
          

          
          sink()
          incProgress(1/TotProgress, detail="...Done fitting the Bayesian linear regression model...")
          
          bayeslincinoerror <- RegressionSingleCI(data = PostBLM1_fit_NoErrors, from = minLim, to = maxLim)
          bayeslincalcinoerror <- as.data.frame(bayeslincinoerror)
          bayeslinciwitherror <- RegressionSingleCI(data = PostBLM1_fit, from = minLim, to = maxLim)
          bayeslincalciwitherror <- as.data.frame(bayeslinciwitherror)
          
          output$bayeslincalibration <- renderPlotly({
            bayeslinfig <- plot_ly(data = calibrationData()
            )
            bayeslinfig <- bayeslinfig %>% 
              add_trace(x = ~calibrationData()$Temperature, 
                        y = ~D47,
                        type = "scatter", 
                        mode = "markers", 
                        marker = list(color = "black"),
                        opacity = 0.5,
                        name = "Raw data",
                        text = as.character(calibrationData()$Sample.Name),
                        hovertemplate = paste(
                          "<b>Sample: %{text}</b><br><br>",
                          "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                          "Δ<sub>47</sub> (‰): %{y}<br>",
                          "Mineralogy: ", as.character(calibrationData()$Mineralogy),"<br>",
                          "Type: ", as.character(calibrationData()$Material),
                          "<extra></extra>")) %>%
              add_ribbons(data = bayeslincalcinoerror,
                          x = ~x,
                          y = ~median_est,
                          ymin = ~ci_lower_est,
                          ymax = ~ci_upper_est,
                          line = list(color = "#ffd166"),
                          fillcolor = "#ffd166",
                          opacity = 0.5,
                          name = "95% CI - no error",
                          hovertemplate = paste(
                            "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                            "Δ<sub>47</sub> (‰): %{y}<br>")) %>%
              add_lines(data = bayeslincalcinoerror,
                        x = ~x,
                        y = ~median_est,
                        name = "Mean estimate - no error",
                        line = list(color = "black", dash = "dash"),
                        hovertemplate = paste(
                          "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                          "Δ<sub>47</sub> (‰): %{y}<br>")) %>%
              add_ribbons(data = bayeslincalciwitherror,
                          x = ~x,
                          y = ~median_est,
                          ymin = ~ci_lower_est,
                          ymax = ~ci_upper_est,
                          line = list(color = "#446455"),
                          fillcolor = "#446455",
                          opacity = 0.5,
                          name = "95% CI - with error",
                          hovertemplate = paste(
                            "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                            "Δ<sub>47</sub> (‰): %{y}<br>")) %>%
              add_lines(data = bayeslincalciwitherror,
                        x = ~x,
                        y = ~median_est,
                        name = "Mean estimate - with error",
                        line = list(color = "#446455", dash = "dash"),
                        hovertemplate = paste(
                          "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                          "Δ<sub>47</sub> (‰): %{y}<br>"))
            bayeslinfig <- bayeslinfig %>% layout(title = "<b> Bayesian linear calibration model </b>",
                                                  legend=list(title=list(text="Legend")),
                                                  xaxis = list(title = "Temperature (10<sup>6</sup>/T<sup>2</sup>)", hoverformat = ".1f"), 
                                                  yaxis = list(title = "Δ<sub>47</sub> (‰)", hoverformat = ".3f"))
            
            return(bayeslinfig)
          })
          
          addWorksheet(wb, "Bayesian model no errors") # Add a blank sheet
          addWorksheet(wb, "Bayesian model no errors CI") # Add a blank sheet 
          addWorksheet(wb, "Bayesian model with errors") # Add a blank sheet
          addWorksheet(wb, "Bayesian model with errors CI") # Add a blank sheet 
          
          bayeslincalcinoerror2 <- bayeslincalcinoerror
          names(bayeslincalcinoerror2) <- c("10^6/T^2", "D47_median_est",	"D47_ci_lower_est",	"D47_ci_upper_est")
          
          bayeslincalciwitherror2 <- bayeslincalciwitherror
          names(bayeslincalciwitherror2) <- c("10^6/T^2", "D47_median_est",	"D47_ci_lower_est",	"D47_ci_upper_est")
          

          writeData(wb, sheet = "Bayesian model no errors", PostBLM1_fit_NoErrors) # Write regression data
          writeData(wb, sheet = "Bayesian model no errors CI", bayeslincalcinoerror2)
          
          writeData(wb, sheet = "Bayesian model with errors", PostBLM1_fit) # Write regression data
          writeData(wb, sheet = "Bayesian model with errors CI", bayeslincalciwitherror2)
          
          
          ##For the Bayesian sheet
          conv_BLM <- summary(bayeslincals$BLM1_fit_NoErrors)$summary
          conv_BLM_errors <- summary(bayeslincals$BLM1_fit)$summary
          conv_BLM <- cbind.data.frame(parameter=row.names(conv_BLM),conv_BLM)
          conv_BLM_errors <- cbind.data.frame(parameter=row.names(conv_BLM_errors),conv_BLM_errors)
          
          
          addWorksheet(wb3, "Bayesian model no errors") # Add a blank sheet
          addWorksheet(wb3, "Bayesian model with errors") # Add a blank sheet
          writeData(wb3, sheet = "Bayesian model no errors", conv_BLM) # Write regression data
          writeData(wb3, sheet = "Bayesian model with errors", conv_BLM_errors) # Write regression data
          
          ##For the posterior sheet
          
          addWorksheet(wb4, "Bayesian model no errors") # Add a blank sheet
          addWorksheet(wb4, "Bayesian model with errors") # Add a blank sheet
          writeData(wb4, sheet = "Bayesian model no errors", PostBLM1_fit_NoErrors) # Write regression data
          writeData(wb4, sheet = "Bayesian model with errors", PostBLM1_fit) # Write regression data
          
          
          outBLM <- summary(bayeslincals$BLM1_fit_NoErrors)$summary
          
          output$blinnoerr <- renderPrint({
            round(outBLM[c(1:2),c(1,2)],7)
          })
          
          outBLMerrors <- summary(bayeslincals$BLM1_fit)$summary

          output$blinwerr <- renderPrint({
             round(outBLMerrors[c(1:2),c(1,2)],7)
          })
          
          outBLMM <- summary(bayeslincals$BLM3_fit)$summary
          outBLMM <- as.data.frame(outBLMM[grep("alpha|beta", row.names(outBLMM)),] )

          output$blinmwerr <- renderPrint(
            round(outBLMM[,c(1,2)],7)
          ) 

          outBLMM2 <- PostBLM3_fit
          outBLMM2 <- as.data.frame(outBLMM2[,grep("alpha|beta", colnames(outBLMM2))] )

          nmat <- ncol(outBLMM2)/2
          
          bayeslmminciwitherror <- lapply(seq(1, ncol(outBLMM2), by=2), function(x){
            subset <- outBLMM2[,c(x,x+1)]
            colnames(subset) <- c("alpha", "beta")
            RegressionSingleCI(data = subset, from = minLim, to = maxLim)[[1]]
          })
          
          bayeslmminciwitherror <- rbindlist(bayeslmminciwitherror, idcol="Material")
          bayeslmmincalciwitherror <- as.data.frame(bayeslmminciwitherror)

          addWorksheet(wb, "Bayesian mixed w errors") # Add a blank sheet
          addWorksheet(wb, "Bayesian mixed w errors CI") # Add a blank sheet 
          
          bayeslmmincalciwitherror2 <- bayeslmmincalciwitherror
          names(bayeslmmincalciwitherror2) <- c("material", "10^6/T^2", "D47_median_est",	"D47_ci_lower_est",	"D47_ci_upper_est")
          
          writeData(wb, sheet = "Bayesian mixed w errors", PostBLM3_fit) # Write regression data
          writeData(wb, sheet = "Bayesian mixed w errors CI", bayeslmmincalciwitherror2)
          
          ##For the Bayesian sheet
          conv_BLMM <- summary(bayeslincals$BLM3_fit)$summary
          conv_BLMM <- cbind.data.frame(parameter=row.names(conv_BLMM),conv_BLMM)
          
          addWorksheet(wb3, "Bayesian mixed w errors") # Add a blank sheet
          writeData(wb3, sheet = "Bayesian mixed w errors", conv_BLMM) # Write regression data
          
          addWorksheet(wb4, "Bayesian mixed w errors") # Add a blank sheet
          writeData(wb4, sheet = "Bayesian mixed w errors", PostBLM3_fit) # Write regression data
          
          output$bayesmixedcalibration <- renderPlotly({
            bayesmixedfig <- plot_ly(data = calibrationData()
            )
            bayesmixedfig <- bayesmixedfig %>% 
              add_trace(x = ~calibrationData()$Temperature, 
                        y = ~D47,
                        type = "scatter", 
                        mode = "markers", 
                        marker = list(color = "black"),
                        opacity = 0.5,
                        name = "Raw data",
                        text = as.character(calibrationData()$Sample.Name),
                        hovertemplate = paste(
                          "<b>Sample: %{text}</b><br><br>",
                          "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                          "Δ<sub>47</sub> (‰): %{y}<br>",
                          "Mineralogy: ", as.character(calibrationData()$Mineralogy),"<br>",
                          "Type: ", as.character(calibrationData()$Material),
                          "<extra></extra>")) 
            colsTarget <- rainbow(nmat)
            materials <- unique(bayeslmmincalciwitherror$Material)
            
            for(y in 1:nmat){
              
              bayesmixedfig <- bayesmixedfig %>% 
                          add_ribbons(data = bayeslmmincalciwitherror[bayeslmmincalciwitherror$Material==materials[y],],
                          x = ~x,
                          y = ~median_est,
                          ymin = ~ci_lower_est,
                          ymax = ~ci_upper_est,
                          line = list(color = colsTarget[y]),
                          fillcolor = colsTarget[y],
                          opacity = 0.5,
                          name = "95% CI - with error",
                          hovertemplate = paste(
                            "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                            "Δ<sub>47</sub> (‰): %{y}<br>")) %>%
              add_lines(data = bayeslmmincalciwitherror[bayeslmmincalciwitherror$Material==materials[y],],
                        x = ~x,
                        y = ~median_est,
                        name = "Mean estimate - with error",
                        line = list(color = colsTarget[y], dash = "dash"),
                        hovertemplate = paste(
                          "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                          "Δ<sub>47</sub> (‰): %{y}<br>"))
            }
            
            bayesmixedfig <- bayesmixedfig %>% layout(title = "<b> Bayesian mixed model </b>",
                                                  legend=list(title=list(text="Legend")),
                                                  xaxis = list(title = "Temperature (10<sup>6</sup>/T<sup>2</sup>)", hoverformat = ".1f"), 
                                                  yaxis = list(title = "Δ<sub>47</sub> (‰)", hoverformat = ".3f"))
            
            return(bayesmixedfig)
          })
          
          lootab <- attr(bayeslincals,"loo")
          
          cat("Model comparison in loo\n\n")
          cat("")
          print(lootab)
          
          # cat(paste0("\nBayesian linear model complete \n *with errors \n   *R^2=", round(attr(bayeslincals,"R2s")[1,2],4),
          #            " (95% CI, ",round(attr(bayeslincals,"R2s")[1,3],4),"-",round(attr(bayeslincals,"R2s")[1,4],4),")",
          #            "\n   *DIC=", round(attr(bayeslincals,"DICs")[1],4),
          #            "\n *without errors\n   *R^2=", round(attr(bayeslincals,"R2s")[2,2],4),
          #            " (95% CI, ",round(attr(bayeslincals,"R2s")[2,3],4),"-",round(attr(bayeslincals,"R2s")[2,4],4),")",
          #            "\n   *DIC=", round(attr(bayeslincals,"DICs")[2],4),
          #            "\n Bayesian linear mixed model complete\n   *R^2=", round(attr(bayeslincals,"R2s")[3,2],4),
          #            " (95% CI, ",round(attr(bayeslincals,"R2s")[3,3],4),"-",round(attr(bayeslincals,"R2s")[3,4],4),")",
          #            "\n   *DIC=", round(attr(bayeslincals,"DICs")[3],4)
          # )
          #)

      }
        
      })
      
    }
   }
    
  })
  
  output$modresults <- renderPrint({
    modresult()
  })
  
  output$downloadcalibrations <- downloadHandler(
    filename = function() { 
      paste("Calibration_output_", Sys.time(), ".xlsx", sep="")
    },
    
    content = function(file) {
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  
  output$downloadBayesian <- downloadHandler(
    filename = function() { 
      paste("Bayesian_output_", Sys.time(), ".xlsx", sep="")
    },
    
    content = function(file) {
      saveWorkbook(wb3, file, overwrite = TRUE)
    }
  )
  
  output$downloadPosteriorCalibration <- downloadHandler(
    filename = function() { 
      paste("Bayesian_posterior_output_", Sys.time(), ".xlsx", sep="")
    },
    
    content = function(file) {
      saveWorkbook(wb4, file, overwrite = TRUE)
    }
  )
  
  
  ########################### Example code
  #    output$calstorun<-renderText({
  #      return(caltorun)
  #})
  ###################################################
  
  
  # Calibration plots tab
  
  observe({
    minlength <- length(unique(calibrationData()$Mineralogy))
    if( !all(is.na(calibrationData()$Mineralogy)) == TRUE ){
    output$rawcaldata <- renderPlotly({
      rawcalfig <- plot_ly(calibrationData(), 
                           x = ~Temperature, 
                           y = ~D47, 
                           type = "scatter", 
                           mode = "markers", 
                           #linetype = ~as.factor(Material), 
                           color = ~as.factor(Mineralogy),
                           colors = viridis_pal(option = "D", end = 0.9)(minlength),
                           opacity = 0.6,
                           error_y = ~list(array = ~D47error, color = "#000000"),
                           error_x = ~list(array = ~TempError, color = "#000000"),
                           text = as.character(calibrationData()$Sample.Name),
                           hovertemplate = paste(
                             "<b>Sample: %{text}</b><br><br>",
                             "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                             "Δ<sub>47</sub> (‰): %{y}<br>",
                             "Mineralogy: ", as.character(calibrationData()$Mineralogy),"<br>",
                             "Type: ", as.character(calibrationData()$Material),
                             "<extra></extra>"))
      rawcalfig <- rawcalfig %>% layout(title = "<b> Raw calibration data from user input </b>",
                                        legend=list(title=list(text="Material and mineralogy")),
                                        xaxis = list(title = "Temperature (10<sup>6</sup>/T<sup>2</sup>)"), 
                                        yaxis = list(title = "Δ<sub>47</sub> (‰)", hoverformat = ".3f"))
      
      return(rawcalfig)
    })
    }else{
      output$rawcaldata <- renderPlotly({
        rawcalfig <- plot_ly(calibrationData(), 
                             x = ~Temperature, 
                             y = ~D47, 
                             type = "scatter", 
                             mode = "lines+markers", 
                             linetype = ~as.factor(Material), 
                             color = ~as.factor(Material),
                             colors = viridis_pal(option = "D", end = 0.9)(length(unique(calibrationData()$Material))),
                             opacity = 0.6,
                             error_y = ~list(array = ~D47error, color = "#000000"),
                             error_x = ~list(array = ~TempError, color = "#000000"),
                             text = as.character(calibrationData()$Sample.Name),
                             hovertemplate = paste(
                               "<b>Sample: %{text}</b><br><br>",
                               "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                               "Δ<sub>47</sub> (‰): %{y}<br>",
                               "Type: ", as.character(calibrationData()$Material),
                               "<extra></extra>"))
        rawcalfig <- rawcalfig %>% layout(title = "<b> Raw calibration data from user input </b>",
                                          legend=list(title=list(text="Material")),
                                          xaxis = list(title = "Temperature (10<sup>6</sup>/T<sup>2</sup>)"), 
                                          yaxis = list(title = "Δ<sub>47</sub> (‰)", hoverformat = ".3f"))
        
        return(rawcalfig)
      })
    }
  })
  
  # Reconstruction tab
  
  output$BayClump_reconstruction_template.csv <- downloadHandler(
    filename = "BayClump_reconstruction_template.csv",
    content = function(file) {
      write.csv(BayClump_reconstruction_template, file, row.names = FALSE)
    }
  )

  reconstructionData = reactive({
    req(input$reconstructiondata)
    n_rows = length(count.fields(input$reconstructiondata$datapath))
    df_out = read.csv(input$reconstructiondata$datapath)
    return(df_out)
  })
  
  if(exists("wb2")) rm(wb2) # Delete any existing workbook in preparation for new results
  wb2 <- createWorkbook("reconstruction output") # Prepare a workbook for reconstruction outputs
  
  observe({
    output$contents2 <- renderTable({
      recsummary <- reconstructionData() %>%
        summarize(
          "Unique samples" = length(unique(reconstructionData()$Sample)),
          "Total replicates" = sum(reconstructionData()$N),
          "Materials" = length(unique(reconstructionData()$Material))
        )
      return(recsummary)
    }, 
    rownames=FALSE, options = list(pageLength = 1, info = FALSE)
    )
  })
  
  recresult <- eventReactive(input$runrec, {
    
    if(is.null(input$reconstructiondata)) {print(noquote("Please upload reconstruction data first"))}
    if(!is.null(input$reconstructiondata)) {
      
      recData <- NULL
      recData <- reconstructionData()
      #recData$Material <- factor(recData$Material, labels = seq(1:length(unique(recData$Material))))
      #if(min(as.numeric(recData$Material)) != 1){
      #  print(noquote("The sequence of Materials now start from 1"))
      #}
  
      
      # Remove existing worksheets from wb2 on "run" click, if any
      if("Linear" %in% names(wb2) == TRUE) 
      {removeWorksheet(wb2, "Linear")}
      if("Linear w no uncertainty" %in% names(wb2) == TRUE) 
      {removeWorksheet(wb2, "Linear w no uncertainty")}
      
      
      if("Inverse linear" %in% names(wb2) == TRUE) 
      {removeWorksheet(wb2, "Inverse linear")}
      if("Inverse linear w no uncertainty" %in% names(wb2) == TRUE) 
      {removeWorksheet(wb2, "Inverse linear w no uncertainty")}
      
      if("York" %in% names(wb2) == TRUE) 
      {removeWorksheet(wb2, "York")}
      if("York w no uncertainty" %in% names(wb2) == TRUE) 
      {removeWorksheet(wb2, "York w no uncertainty")}
      
      if("Deming" %in% names(wb2) == TRUE) 
      {removeWorksheet(wb2, "Deming")}
      if("Deming w no uncertainty" %in% names(wb2) == TRUE) 
      {removeWorksheet(wb2, "Deming w no uncertainty")}
      
      if("Bayes and error" %in% names(wb2) == TRUE) 
      {removeWorksheet(wb2, "Bayes and error") }
      if("Bayes w error no uncertainty" %in% names(wb2) == TRUE) 
      {removeWorksheet(wb2, "Bayes w error no uncertainty") }
      
      if("Bayesian predictions" %in% names(wb2) == TRUE) 
      {removeWorksheet(wb2, "Bayesian predictions")}
      if("Bayesian mixed model" %in% names(wb2) == TRUE)
      {removeWorksheet(wb2, "Bayesian mixed model")}
      if("Bayesian linear model, errors" %in% names(wb2) == TRUE)
      {removeWorksheet(wb2, "Bayesian linear model, errors")}
      if("Bayesian linear model" %in% names(wb2) == TRUE)
      {removeWorksheet(wb2, "Bayesian linear model")}
      if("Bayesian linear mixed model" %in% names(wb2) == TRUE)
      {removeWorksheet(wb2, "Bayesian linear mixed model")}
      
      
      ##Also for the Bayesian posterior sheet
      if("Bayesian model no errors" %in% names(wb5) == TRUE) 
      {removeWorksheet(wb5, "Bayesian model no errors")}
      if("Bayesian model with errors" %in% names(wb5) == TRUE) 
      {removeWorksheet(wb5, "Bayesian model with errors") }
      if("Bayesian mixed w errors" %in% names(wb5) == TRUE)
      {removeWorksheet(wb5, "Bayesian mixed w errors")}
      if("Bayesian linear mixed model" %in% names(wb5) == TRUE)
      {removeWorksheet(wb5, "Bayesian linear mixed model")}
      
      if("Bayesian model no errors (C)" %in% names(wb5) == TRUE) 
      {removeWorksheet(wb5, "Bayesian model no errors (C)")}
      if("Bayesian model with errors (C)" %in% names(wb5) == TRUE) 
      {removeWorksheet(wb5, "Bayesian model with errors (C)") }
      if("Bayesian mixed w errors, Convergence" %in% names(wb5) == TRUE)
      {removeWorksheet(wb5, "Bayesian mixed w errors, Convergence")}
      if("Bayesian linear mixed model (C)" %in% names(wb5) == TRUE)
      {removeWorksheet(wb5, "Bayesian linear mixed model (C)")}
      
        totalModelsRecs <- c(input$simulateLM_measuredRec, input$simulateLM_inverseweightsRec, input$simulateYork_measuredRec,
                         input$simulateDemingRec, input$BayesianCalibrationsRec)
        totalModelsRecs <- length(which(totalModelsRecs==T))

        withProgress(message = "Running selected reconstructions, please wait", {
          
          #OLS
          if( input$simulateLM_measuredRec) {
            
            if(is.null(lmcals)) { print(noquote("Please run the calibration step for linear models first")) }else{
              
              sink("out/LMpredictions.txt", type = "output")
            
              lmrec <<-  predictTc(calData = calData,
                                 recData = recData,
                                 obCal = lmcals,
                                 clumpedClassic = input$simpleInversion)
            
            sink()
            df1 <- lmrec
            
            names(df1) <- c("Sample", "Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "LW, 95% CI", "Up, 95% CI")
            rownames(df1) <- NULL
            
           
            output$lmrecswun <- renderTable({
              df1$`Δ47 (‰)` <- formatC(df1$`Δ47 (‰)`, digits = 3, format = "f")
              df1$`Δ47 (‰) error` <- formatC(df1$`Δ47 (‰) error`, digits = 4, format = "f")
              df1$`Temperature (°C)` <- formatC(df1$`Temperature (°C)`, digits = 1, format = "f")
              df1$`LW, 95% CI` <- formatC(df1$`LW, 95% CI`, digits = 3, format = "f")
              df1$`Up, 95% CI` <- formatC(df1$`Up, 95% CI`, digits = 3, format = "f")
              head(df1)
            },
              caption = "Linear model",
              caption.placement = getOption("xtable.caption.placement", "top"),
              rownames = FALSE,
              spacing = "m",
              align = "c"
            )
            
            addWorksheet(wb2, "Linear") # Add a blank sheet
            writeData(wb2, sheet = "Linear", df1) # Write reconstruction data
            print(noquote("Linear reconstruction complete"))
            
          
            }
            }
          
          #Inverse weighted linear model 
          if( input$simulateLM_inverseweightsRec ) {
            if(is.null(lminversecals)) { print(noquote("Please run the calibration step for weighted OLS models first")) }else{
              
            sink("out/wLMpredictions.txt", type = "output")
              
              lminverserec <<- predictTc(calData = calData,
                        recData = recData,
                        obCal = lminversecals,
                        clumpedClassic = input$simpleInversion)

              sink()
            incProgress(1/totalModelsRecs, detail="...Done fitting the weighted OLS...")
            
            
            df3<-lminverserec
            
            names(df3) <- c("Sample", "Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "LW, 95% CI", "Up, 95% CI")
            rownames(df3) <- NULL
            
            output$lminverserecswun <- renderTable({
              df3$`Δ47 (‰)` <- formatC(df3$`Δ47 (‰)`, digits = 3, format = "f")
              df3$`Δ47 (‰) error` <- formatC(df3$`Δ47 (‰) error`, digits = 4, format = "f")
              df3$`Temperature (°C)` <- formatC(df3$`Temperature (°C)`, digits = 1, format = "f")
              df3$`LW, 95% CI` <- formatC(df3$`LW, 95% CI`, digits = 3, format = "f")
              df3$`Up, 95% CI` <- formatC(df3$`Up, 95% CI`, digits = 3, format = "f")
              head(df3)
            },
              caption = "Inverse weighted linear model",
              caption.placement = getOption("xtable.caption.placement", "top"),
            rownames = FALSE,
            spacing = "m",
            align = "c"
              
            )
            
            addWorksheet(wb2, "Inverse linear") # Add a blank sheet
            writeData(wb2, sheet = "Inverse linear", df3) # Write reconstruction data
            print(noquote("Inverse weighted linear reconstruction complete"))
            
            
          
            }
            }
          
          # York regression
          if( input$simulateYork_measuredRec) {
            if(is.null(yorkcals)) { print(noquote("Please run the calibration step for York models first")) }else{
              
              sink("out/Yorkpredictions.txt", type = "output")
              


            yorkrec <<-  predictTc(calData = calData,
                                   recData = recData,
                                   obCal = yorkcals,
                                   clumpedClassic = input$simpleInversion)
            
            
              sink()
            incProgress(1/totalModelsRecs, detail="...Done fitting the York regression...")
            

            df5 <- yorkrec
            names(df5) <- c("Sample", "Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "LW, 95% CI", "Up, 95% CI")
            rownames(df5) <- NULL
            
            
            output$yorkrecswun <- renderTable({
              
              df5$`Δ47 (‰)` <- formatC(df5$`Δ47 (‰)`, digits = 3, format = "f")
              df5$`Δ47 (‰) error` <- formatC(df5$`Δ47 (‰) error`, digits = 4, format = "f")
              df5$`Temperature (°C)` <- formatC(df5$`Temperature (°C)`, digits = 1, format = "f")
              df5$`LW, 95% CI` <- formatC(df5$`LW, 95% CI`, digits = 3, format = "f")
              df5$`Up, 95% CI` <- formatC(df5$`Up, 95% CI`, digits = 3, format = "f")
              head(df5)
            },
              caption = "York regression",
              caption.placement = getOption("xtable.caption.placement", "top"),
            rownames = FALSE,
            spacing = "m",
            align = "c"
              
            )
            
            addWorksheet(wb2, "York") # Add a blank sheet
            writeData(wb2, sheet = "York", df5) # Write reconstruction data
            print(noquote("York reconstruction complete"))

          
            }
            }
          
          # Deming regression
          if( input$simulateDemingRec) {
            if(is.null(demingcals) ) { print(noquote("Please run the calibration step for Deming models first")) }else{
              
              sink("out/Demingpredictions.txt", type = "output")
              
              
              
            demingrec <<- predictTc(calData = calData,
                                    recData = recData,
                                    obCal = demingcals,
                                    clumpedClassic = input$simpleInversion)
            
            
            sink()
            incProgress(1/totalModelsRecs, detail="...Done fitting the Deming regression...")
            
            
            df7 <- demingrec
            names(df7) <- c("Sample", "Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "LW, 95% CI", "Up, 95% CI")
            rownames(df7) <- NULL

            output$demingrecswun <- renderTable({
              
              df7$`Δ47 (‰)` <- formatC(df7$`Δ47 (‰)`, digits = 3, format = "f")
              df7$`Δ47 (‰) error` <- formatC(df7$`Δ47 (‰) error`, digits = 4, format = "f")
              df7$`Temperature (°C)` <- formatC(df7$`Temperature (°C)`, digits = 1, format = "f")
              df7$`LW, 95% CI` <- formatC(df7$`LW, 95% CI`, digits = 3, format = "f")
              df7$`Up, 95% CI` <- formatC(df7$`Up, 95% CI`, digits = 3, format = "f")
              head(df7)
            },
              caption = "Deming regression",
              caption.placement = getOption("xtable.caption.placement", "top"),
            rownames = FALSE,
            spacing = "m",
            align = "c"
              
            )
            
            
            addWorksheet(wb2, "Deming") # Add a blank sheet
            writeData(wb2, sheet = "Deming", df7) # Write reconstruction data
            print(noquote("Deming reconstruction complete"))
 
            }
            }

          #BML
          if(input$BayesianCalibrationsRec){
            
            if(is.null(bayeslincals)) { print(noquote("Please run the calibration step for Bayesian linear models first")) }else{
              
              ##This function runs only Bayesian predictions
              sink("out/bayespredictions.txt", type = "output")
              
              #TPriorMean <- as.numeric(input$TPriorMean) #Not being used rn
              #TPriorSd <- as.numeric(input$TPriorSd) #Not being used rn
              
              infTempBayesian <- BayesianPredictions(calModel = bayeslincals$BLM1_fit_NoErrors,
                                                     recData = recData
                                                    )
              
              infTempBayesian_NE <- BayesianPredictions(calModel = bayeslincals$BLM1_fit_NoErrors,
                                                     recData = recData
              )
              
              ##Need to implement the mixed one...
              
              infTempBayesian_Mixed <- infTempBayesian_NE
              
              sink()

              df0 <- infTempBayesian
              names(df0) <- c("Sample", "Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "LW, 95% CI", "Up, 95% CI")
              rownames(df0) <- NULL
              
              output$BpredictionsErrors <- renderTable({
                df0$`Δ47 (‰)` <- formatC(df0$`Δ47 (‰)`, digits = 3, format = "f")
                df0$`Δ47 (‰) error` <- formatC(df0$`Δ47 (‰) error`, digits = 4, format = "f")
                df0$`Temperature (°C)` <- formatC(df0$`Temperature (°C)`, digits = 1, format = "f")
                df0$`LW, 95% CI` <- formatC(df0$`LW, 95% CI`, digits = 3, format = "f")
                df0$`Up, 95% CI` <- formatC(df0$`Up, 95% CI`, digits = 3, format = "f")
                head(df0)
              },
              caption = "Bayesian predictions (BLM_errors)",
              caption.placement = getOption("xtable.caption.placement", "top"),
              rownames = FALSE,
              spacing = "m",
              align = "c"
              
              )

              
              ##Without errors
              df0.1 <- infTempBayesian_NE
              names(df0.1) <- c("Sample", "Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "LW, 95% CI", "Up, 95% CI")
              rownames(df0.1) <- NULL
              
              output$Bpredictions <- renderTable({
                
                df0.1$`Δ47 (‰)` <- formatC(df0.1$`Δ47 (‰)`, digits = 3, format = "f")
                df0.1$`Δ47 (‰) error` <- formatC(df0.1$`Δ47 (‰) error`, digits = 4, format = "f")
                df0.1$`Temperature (°C)` <- formatC(df0.1$`Temperature (°C)`, digits = 1, format = "f")
                df0.1$`LW, 95% CI` <- formatC(df0.1$`LW, 95% CI`, digits = 3, format = "f")
                df0.1$`Up, 95% CI` <- formatC(df0.1$`Up, 95% CI`, digits = 3, format = "f")
                head(df0.1)
              },
              caption = "Bayesian predictions (BLM without errors)",
              caption.placement = getOption("xtable.caption.placement", "top"),
              rownames = FALSE,
              spacing = "m",
              align = "c"
              
              )
              
              
              df0.2 <- infTempBayesian_Mixed
              names(df0.2) <- c("Sample", "Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "LW, 95% CI", "Up, 95% CI")
              rownames(df0.2) <- NULL
              
              output$BpredictionsBLMM <- renderTable({
                df0.2$`Δ47 (‰)` <- formatC(df0.2$`Δ47 (‰)`, digits = 3, format = "f")
                df0.2$`Δ47 (‰) error` <- formatC(df0.2$`Δ47 (‰) error`, digits = 4, format = "f")
                #df0.2$Material <- formatC(df0.2$Material, digits = 1, format = "f")
                df0.2$`Temperature (°C)` <- formatC(df0.2$`Temperature (°C)`, digits = 1, format = "f")
                df0.2$`LW, 95% CI` <- formatC(df0.2$`LW, 95% CI`, digits = 3, format = "f")
                df0.2$`Up, 95% CI` <- formatC(df0.2$`Up, 95% CI`, digits = 3, format = "f")
                head(df0.2)
              },
              caption = "Bayesian predictions under a Bayesian linear mixed model",
              caption.placement = getOption("xtable.caption.placement", "top"),
              rownames = FALSE,
              spacing = "m",
              align = "c"
              )
              
              
              addWorksheet(wb2, "Bayesian linear model") # Add a blank sheet
              addWorksheet(wb2, "Bayesian linear mixed model") # Add a blank sheet
              addWorksheet(wb2, "Bayesian linear model, errors") # Add a blank sheet
              writeData(wb2, sheet = "Bayesian linear model", df0.1)
              writeData(wb2, sheet = "Bayesian linear model, errors", df0)
              writeData(wb2, sheet = "Bayesian linear mixed model", df0.2)
              

              addWorksheet(wb5, "Bayesian model no errors") # Add a blank sheet
              addWorksheet(wb5, "Bayesian model with errors") # Add a blank sheet
              addWorksheet(wb5, "Bayesian linear mixed model") # Add a blank sheet
              writeData(wb5, sheet = "Bayesian model no errors", infTempBayesian_NE) # Write regression data
              writeData(wb5, sheet = "Bayesian model with errors", infTempBayesian ) # Write regression data
              writeData(wb5, sheet = "Bayesian linear mixed model", infTempBayesian_Mixed ) # Write regression data
              

              print(noquote("Bayesian linear reconstructions complete"))
              
              
            
              
            }
          }
          
          
          
        })
      
    }
  })
  
  
  output$recresults <- renderPrint({
    recresult()
  })
  
  output$downloadreconstructions <- downloadHandler(
    filename = function() { 
      paste("Reconstruction_output_", Sys.time(), ".xlsx", sep="")
    },
    
    content = function(file) {
      saveWorkbook(wb2, file, overwrite = TRUE)
    }
  )
  
  output$downloadreconstructionsPosterior <- downloadHandler(
    filename = function() { 
      paste("Reconstruction_posterior_output_", Sys.time(), ".xlsx", sep="")
    },
    
    content = function(file) {
      saveWorkbook(wb5, file, overwrite = TRUE)
    }
  )
  
  # Manuscript tab
  
  output$msframe <- renderUI({
    tags$iframe(style="height:600px; width:100%", src="essoar.10507995.2.pdf")
  })
  
}
  