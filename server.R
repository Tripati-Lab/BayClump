# Define server logic
server <- function(input, output, session) { 
  options(shiny.maxRequestSize=800*1024^2) 
  
  #Number of generations for Bayesian predictions
  ngenerationsBayesianPredictions = 20000
  
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
           'model1' = return(Petersen),
           'model2' = return(Anderson),
           'model1and2' = return(PetersenAnderson),
           'mycal' = reactiveValues({
             req(input$calibrationdata)
             n_rows = length(count.fields(input$calibrationdata$datapath))
             df_out = read.csv(input$calibrationdata$datapath)
             return(df_out)
           }),
           'all' = reactiveValues({
             req(input$calibrationdata)
             n_rows = length(count.fields(input$calibrationdata$datapath))
             df_out = read.csv(input$calibrationdata$datapath)
             alldat <- rbind(df_out, PetersenAnderson)
             return(alldat)
           })
    )
  })

  
  if(exists("wb")) rm(wb) # Delete any existing workbook in preparation for new results
  wb <- createWorkbook("calibration output") # Prepare a workbook for calibration outputs
  
  observe({
    output$contents <- renderTable({
      calsummary <- calibrationData() %>%
        summarize(
          "Unique samples" = length(unique(calibrationData()$Sample.Name)),
          "Total replicates" = sum(calibrationData()$N),
          "Mineralogies" = length(unique(calibrationData()$Mineralogy)),
          "Materials" = length(unique(calibrationData()$Material))
        )
      return(calsummary)
    }, 
    rownames=FALSE, options = list(pageLength = 1, info = FALSE)
    )
  }) 
  
  modresult <- eventReactive(input$runmods, {
    
    if ('all' %in% input$calibrationdata) {
      print(noquote("Please upload calibration data first"))
    } 
    
    hasMaterial <<- ifelse( is.na(calibrationData()$Material), FALSE, TRUE )
    
    # Update the number of bootstrap replicates to run based on user selection
    #replicates <- ifelse(input$replication == "50", 50,
    #                     ifelse(input$replication == "100", 100,
    #                            ifelse(input$replication == "500", 500,
    #                                   ifelse(input$replication == "1000", 1000, NA))))
    
    replicates <- input$replication
    # Bayesian n generations
    ngenerationsBayes <- 20000#ifelse(input$ngenerationsBayesian == "1000", 1000,
                         #ifelse(input$ngenerationsBayesian == "5000", 5000,
                        #        ifelse(input$ngenerationsBayesian == "10000", 10000,
                         #              ifelse(input$ngenerationsBayesian == "20000", 20000, NA))))
    
    
    # Remove existing worksheets from wb on 'run' click, if any
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
    
    calData <<- NULL
    calData <<- calibrationData()

    # Recode NA or 0 error values to dummy value
    calData$D47error[calData$D47error == 0] <<- 0.000001
    calData$TempError[calData$TempError == 0] <<- 0.000001
    calData$D47error[is.na(calData$D47error)] <<- 0.000001
    calData$TempError[is.na(calData$TempError)] <<- 0.000001
    calData$Material[is.na(calData$Material)] <<- 1
    
    # For future implementation:
   # if(input$uncertainties == "usedaeron") { # Placeholder for Daeron et al. uncertainties
  #    calData$TempError <<- 1
  #  }
    
    if(input$scale == TRUE) {
      calData$Temperature <<- scale(calData$Temperature)
      calData$TempError <<- scale(calData$TempError)
      calData$D47 <<- scale(calData$D47)
      calData$D47error <<- scale(calData$D47error)
    }
    
    lmcals <<- NULL
    lminversecals <<- NULL
    yorkcals <<- NULL
    demingcals <<- NULL
    bayeslincals <<- NULL
    bayesmixedcals <<- NULL
    
    if(input$simulateLM_measured == FALSE &
       input$simulateLM_inverseweights == FALSE &
       input$simulateYork_measured == FALSE &
       input$simulateDeming == FALSE &
       input$simulateBLM_measuredMaterial == FALSE &
       input$simulateBLMM_measuredMaterial == FALSE) {print(noquote("Please select at least one model"))}
    
    if(input$simulateLM_measured != FALSE |
       input$simulateLM_inverseweights != FALSE |
       input$simulateYork_measured != FALSE |
       input$simulateDeming != FALSE |
       input$simulateBLM_measuredMaterial != FALSE |
       input$simulateBLMM_measuredMaterial != FALSE) {
      
      withProgress(message = 'Running selected models, please wait', {
        
        if(input$simulateLM_measured == FALSE) {
        }
        if(input$simulateLM_inverseweights == FALSE) {
        }
        if(input$simulateYork_measured == FALSE) {
        }
        if(input$simulateDeming == FALSE) {
        }
        if(input$simulateBLM_measuredMaterial == FALSE) {
        }
        if(input$simulateBLMM_measuredMaterial == FALSE) {
        }
        
        if(input$simulateLM_measured != FALSE) {
          sink(file = "linmodtext.txt", type = "output")
          lmcals <<- simulateLM_measured(calData, replicates = replicates)
          sink()
          
          lmci <- RegressionSingleCI(data = lmcals, from = min(calData$Temperature), to = max(calData$Temperature))
          lmcalci <- as.data.frame(lmci)
          
          output$lmcalibration <- renderPlotly({
            lmfig <- plot_ly(calibrationData()
            )
            lmfig <- lmfig %>%
              add_trace(x = ~calibrationData()$Temperature, 
                        y = ~D47,
                        type = 'scatter', 
                        mode = 'markers', 
                        marker = list(color = 'black'),
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
                          line = list(color = '#ffd166'),
                          fillcolor = '#ffd166',
                          opacity = 0.5,
                          name = '95% CI',
                          hovertemplate = paste(
                            "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                            "Δ<sub>47</sub> (‰): %{y}<br>")) %>%
              add_lines(data = lmcalci,
                        x = ~x,
                        y = ~median_est,
                        name = 'Median estimate',
                        line = list(color = "black", dash = 'dash'),
                        hovertemplate = paste(
                          "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                          "Δ<sub>47</sub> (‰): %{y}<br>"))
            lmfig <- lmfig %>% layout(title = '<b> Linear calibration model </b>',
                                      legend=list(title=list(text='Legend')),
                                      xaxis = list(title = 'Temperature (10<sup>6</sup>/T<sup>2</sup>)', hoverformat = '.1f'), 
                                      yaxis = list(title = 'Δ<sub>47</sub> (‰)', hoverformat = '.3f'))
            
            return(lmfig)
          })
          
          addWorksheet(wb, "Linear regression") # Add a blank sheet
          addWorksheet(wb, "Linear regression CI") # Add a blank sheet 
          
          lmcalci2 <- lmcalci
          names(lmcalci2) <- c("10^6/T^2", "D47_median_est",	"D47_ci_lower_est",	"D47_ci_upper_est")
          
          writeData(wb, sheet = "Linear regression", lmcals) # Write regression data
          writeData(wb, sheet = "Linear regression CI", lmcalci2)
          
          print(noquote("Linear regression complete"))
          output$lmcal <- renderPrint({
            do.call(rbind.data.frame,apply(lmcals, 2, function(x){
              cbind.data.frame(Median= round(median(x), 4), `Lower 95% CI`=round(quantile(x, 0.025), 4),  `Upper 95% CI`=round(quantile(x, 0.975), 4))
            }))
          })
          
        }
        
        if(input$simulateLM_inverseweights != FALSE) {
          sink(file = "inverselinmodtext.txt", type = "output")
          lminversecals <<- simulateLM_inverseweights(calData, replicates = replicates)
          sink()
          
          lminverseci <- RegressionSingleCI(data = lminversecals, from = min(calData$Temperature), to = max(calData$Temperature))
          lminversecalci <- as.data.frame(lminverseci)
          
          output$lminversecalibration <- renderPlotly({
            lminversefig <- plot_ly(data = calibrationData()
            )
            lminversefig <- lminversefig %>% 
              add_trace(x = ~calibrationData()$Temperature, 
                        y = ~D47,
                        type = 'scatter', 
                        mode = 'markers', 
                        marker = list(color = 'black'),
                        opacity = 0.5,
                        name = 'Raw data',
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
                          line = list(color = '#ffd166'),
                          fillcolor = '#ffd166',
                          opacity = 0.5,
                          name = '95% CI',
                          hovertemplate = paste(
                            "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                            "Δ<sub>47</sub> (‰): %{y}<br>")) %>%
              add_lines(data = lminversecalci,
                        x = ~x,
                        y = ~median_est,
                        name = 'Median estimate',
                        line = list(color = "black", dash = 'dash'),
                        hovertemplate = paste(
                          "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                          "Δ<sub>47</sub> (‰): %{y}<br>"))
            lminversefig <- lminversefig %>% layout(title = '<b> Inverse linear calibration model </b>',
                                                    legend=list(title=list(text='Legend')),
                                                    xaxis = list(title = 'Temperature (10<sup>6</sup>/T<sup>2</sup>)', hoverformat = '.1f'), 
                                                    yaxis = list(title = 'Δ<sub>47</sub> (‰)', hoverformat = '.3f'))
            
            return(lminversefig)
          })
          
          addWorksheet(wb, "Inverse linear regression") # Add a blank sheet
          addWorksheet(wb, "Inverse linear regression CI") # Add a blank sheet 
          
          lminversecalci2 <- lminversecalci
          names(lminversecalci2) <- c("10^6/T^2", "D47_median_est",	"D47_ci_lower_est",	"D47_ci_upper_est")
          
          writeData(wb, sheet = "Inverse linear regression", lminversecals) # Write regression data
          writeData(wb, sheet = "Inverse linear regression CI", lminversecalci2)
          
          print(noquote("Inverse linear regression complete"))
          output$lminversecal <- renderPrint({
            do.call(rbind.data.frame,apply(lminversecals, 2, function(x){
              cbind.data.frame(Median= round(median(x), 4), `Lower 95% CI`=round(quantile(x, 0.025), 4),  `Upper 95% CI`=round(quantile(x, 0.975), 4))
            }))
          })
          
        }
        
        if(input$simulateYork_measured != FALSE) {
          sink(file = "yorkmodtext.txt", type = "output")
          yorkcals <<- simulateYork_measured(calData, replicates = replicates)
          sink()
          
          yorkci <- RegressionSingleCI(data = yorkcals, from = min(calData$Temperature), to = max(calData$Temperature))
          yorkcalci <- as.data.frame(yorkci)
          
          output$yorkcalibration <- renderPlotly({
            yorkfig <- plot_ly(data = calibrationData()
            )
            yorkfig <- yorkfig %>% 
              add_trace(x = ~calibrationData()$Temperature, 
                        y = ~D47,
                        type = 'scatter', 
                        mode = 'markers', 
                        marker = list(color = 'black'),
                        opacity = 0.5,
                        name = 'Raw data',
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
                          line = list(color = '#ffd166'),
                          fillcolor = '#ffd166',
                          opacity = 0.5,
                          name = '95% CI',
                          hovertemplate = paste(
                            "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                            "Δ<sub>47</sub> (‰): %{y}<br>")) %>%
              add_lines(data = yorkcalci,
                        x = ~x,
                        y = ~median_est,
                        name = 'Median estimate',
                        line = list(color = "black", dash = 'dash'),
                        hovertemplate = paste(
                          "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                          "Δ<sub>47</sub> (‰): %{y}<br>"))
            yorkfig <- yorkfig %>% layout(title = '<b> York calibration model </b>',
                                          legend=list(title=list(text='Legend')),
                                          xaxis = list(title = 'Temperature (10<sup>6</sup>/T<sup>2</sup>)', hoverformat = '.1f'), 
                                          yaxis = list(title = 'Δ<sub>47</sub> (‰)', hoverformat = '.3f'))
            
            return(yorkfig)
          })
          
          addWorksheet(wb, "York regression") # Add a blank sheet
          addWorksheet(wb, "York regression CI") # Add a blank sheet 
          
          yorkcalci2 <- yorkcalci
          names(yorkcalci2) <- c("10^6/T^2", "D47_median_est",	"D47_ci_lower_est",	"D47_ci_upper_est")
          
          writeData(wb, sheet = "York regression", yorkcals) # Write regression data
          writeData(wb, sheet = "York regression CI", yorkcalci2)
          
          print(noquote("York regression complete"))
          output$york <- renderPrint({
            do.call(rbind.data.frame,apply(yorkcals, 2, function(x){
              cbind.data.frame(Median= round(median(x), 4), `Lower 95% CI`=round(quantile(x, 0.025), 4),  `Upper 95% CI`=round(quantile(x, 0.975), 4))
            }))
          })
          
        }
        
        if(input$simulateDeming != FALSE) {
          sink(file = "demingmodtext.txt", type = "output")
          demingcals <<- simulateDeming(calData, replicates = replicates)
          sink()
          
          demingci <- RegressionSingleCI(data = demingcals, from = min(calData$Temperature), to = max(calData$Temperature))
          demingcalci <- as.data.frame(demingci)
          
          output$demingcalibration <- renderPlotly({
            demingfig <- plot_ly(data = calibrationData()
            )
            demingfig <- demingfig %>% 
              add_trace(x = ~calibrationData()$Temperature, 
                        y = ~D47,
                        type = 'scatter', 
                        mode = 'markers', 
                        marker = list(color = 'black'),
                        opacity = 0.5,
                        name = 'Raw data',
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
                          line = list(color = '#ffd166'),
                          fillcolor = '#ffd166',
                          opacity = 0.5,
                          name = '95% CI',
                          hovertemplate = paste(
                            "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                            "Δ<sub>47</sub> (‰): %{y}<br>")) %>%
              add_lines(data = demingcalci,
                        x = ~x,
                        y = ~median_est,
                        name = 'Median estimate',
                        line = list(color = "black", dash = 'dash'),
                        hovertemplate = paste(
                          "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                          "Δ<sub>47</sub> (‰): %{y}<br>"))
            demingfig <- demingfig %>% layout(title = '<b> Deming calibration model </b>',
                                              legend=list(title=list(text='Legend')),
                                              xaxis = list(title = 'Temperature (10<sup>6</sup>/T<sup>2</sup>)', hoverformat = '.1f'), 
                                              yaxis = list(title = 'Δ<sub>47</sub> (‰)', hoverformat = '.3f'))
            
            return(demingfig)
          })
          
          addWorksheet(wb, "Deming regression") # Add a blank sheet
          addWorksheet(wb, "Deming regression CI") # Add a blank sheet 
          
          demingcalci2 <- demingcalci
          names(demingcalci2) <- c("10^6/T^2", "D47_median_est",	"D47_ci_lower_est",	"D47_ci_upper_est")
          
          writeData(wb, sheet = "Deming regression", demingcals) # Write regression data
          writeData(wb, sheet = "Deming regression CI", demingcalci2)
          
          print(noquote("Deming regression complete"))
          output$deming <- renderPrint({
            do.call(rbind.data.frame,apply(demingcals, 2, function(x){
              cbind.data.frame(Median= round(median(x), 4), `Lower 95% CI`=round(quantile(x, 0.025), 4),  `Upper 95% CI`=round(quantile(x, 0.975), 4))
            }))
          })
          
        }
        
        #     checkboxInput("linear", "Linear model", FALSE),
        if(input$simulateBLM_measuredMaterial != FALSE) {
          sink(file = "Bayeslinmodtext.txt", type = "output")
          bayeslincals <<- simulateBLM_measuredMaterial(calData, replicates = replicates, isMixed=F, generations=ngenerationsBayes)
          sink()
          
          bayeslincinoerror <- RegressionSingleCI(data = bayeslincals$BLM_Measured_no_errors, from = min(calData$Temperature), to = max(calData$Temperature))
          bayeslincalcinoerror <- as.data.frame(bayeslincinoerror)
          bayeslinciwitherror <- RegressionSingleCI(data = bayeslincals$BLM_Measured_errors, from = min(calData$Temperature), to = max(calData$Temperature))
          bayeslincalciwitherror <- as.data.frame(bayeslinciwitherror)
          
          output$bayeslincalibration <- renderPlotly({
            bayeslinfig <- plot_ly(data = calibrationData()
            )
            bayeslinfig <- bayeslinfig %>% 
              add_trace(x = ~calibrationData()$Temperature, 
                        y = ~D47,
                        type = 'scatter', 
                        mode = 'markers', 
                        marker = list(color = 'black'),
                        opacity = 0.5,
                        name = 'Raw data',
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
                          line = list(color = '#ffd166'),
                          fillcolor = '#ffd166',
                          opacity = 0.5,
                          name = '95% CI - no error',
                          hovertemplate = paste(
                            "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                            "Δ<sub>47</sub> (‰): %{y}<br>")) %>%
              add_lines(data = bayeslincalcinoerror,
                        x = ~x,
                        y = ~median_est,
                        name = 'Median estimate - no error',
                        line = list(color = "black", dash = 'dash'),
                        hovertemplate = paste(
                          "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                          "Δ<sub>47</sub> (‰): %{y}<br>")) %>%
              add_ribbons(data = bayeslincalciwitherror,
                          x = ~x,
                          y = ~median_est,
                          ymin = ~ci_lower_est,
                          ymax = ~ci_upper_est,
                          line = list(color = '#446455'),
                          fillcolor = '#446455',
                          opacity = 0.5,
                          name = '95% CI - with error',
                          hovertemplate = paste(
                            "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                            "Δ<sub>47</sub> (‰): %{y}<br>")) %>%
              add_lines(data = bayeslincalciwitherror,
                        x = ~x,
                        y = ~median_est,
                        name = 'Median estimate - with error',
                        line = list(color = "#446455", dash = 'dash'),
                        hovertemplate = paste(
                          "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                          "Δ<sub>47</sub> (‰): %{y}<br>"))
            bayeslinfig <- bayeslinfig %>% layout(title = '<b> Bayesian linear calibration model </b>',
                                                  legend=list(title=list(text='Legend')),
                                                  xaxis = list(title = 'Temperature (10<sup>6</sup>/T<sup>2</sup>)', hoverformat = '.1f'), 
                                                  yaxis = list(title = 'Δ<sub>47</sub> (‰)', hoverformat = '.3f'))
            
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
          
          writeData(wb, sheet = "Bayesian model no errors", bayeslincals$BLM_Measured_no_errors) # Write regression data
          writeData(wb, sheet = "Bayesian model no errors CI", bayeslincalcinoerror2)
          
          writeData(wb, sheet = "Bayesian model with errors", bayeslincals$BLM_Measured_errors) # Write regression data
          writeData(wb, sheet = "Bayesian model with errors CI", bayeslincalciwitherror2)
          
          print(noquote("Bayesian linear model complete"))
          
          output$blinnoerr <- renderPrint({
            
            do.call(rbind.data.frame,apply(bayeslincals$BLM_Measured_no_errors, 2, function(x){
              cbind.data.frame(Median= round(median(x), 4), `Lower 95% CI`=round(quantile(x, 0.025), 4),  `Upper 95% CI`=round(quantile(x, 0.975), 4))
            }))
            
          })
          
          output$blinwerr <- renderPrint({
            
            do.call(rbind.data.frame,apply(bayeslincals$BLM_Measured_errors, 2, function(x){
              cbind.data.frame(Median= round(median(x), 4), `Lower 95% CI`=round(quantile(x, 0.025), 4),  `Upper 95% CI`=round(quantile(x, 0.975), 4))
            }))
            
          })
          
        }
        
        if(length(unique(calData$Material)) < 2 & input$simulateBLMM_measuredMaterial != FALSE) {
          print(noquote("Bayesian linear mixed models require multiple materials"))
        }
        
        if(input$simulateBLMM_measuredMaterial != FALSE) { #length(unique(calData$Material)) >= 2 &
          calData$MaterialName <<- calData$Material
          calData$Material <<- as.factor(as.numeric(as.factor(calData$Material)))
          
          sink(file = "Bayesmixmodtext.txt", type = "output")
          bayesmixedcals <- simulateBLM_measuredMaterial(data=calData, replicates = replicates, isMixed = T, generations=ngenerationsBayes)
          sink()
          
          bayeslmminciwitherror <- RegressionSingleCI(data = bayesmixedcals$BLMM_Measured_errors, from = min(calData$Temperature), to = max(calData$Temperature))
          bayeslmmincalciwitherror <- as.data.frame(bayeslmminciwitherror)

          addWorksheet(wb, "Bayesian mixed w errors") # Add a blank sheet
          addWorksheet(wb, "Bayesian mixed w errors CI") # Add a blank sheet 
          
          bayeslmmincalciwitherror2 <- bayeslmmincalciwitherror
          names(bayeslmmincalciwitherror2) <- c("10^6/T^2", "D47_median_est",	"D47_ci_lower_est",	"D47_ci_upper_est")
          
          writeData(wb, sheet = "Bayesian mixed w errors", bayesmixedcals$BLMM_Measured_errors) # Write regression data
          writeData(wb, sheet = "Bayesian mixed w errors CI", bayeslmmincalciwitherror2)
          
          output$bayesmixedcalibration <- renderPlotly({
            bayesmixedfig <- plot_ly(data = calibrationData()
            )
            bayesmixedfig <- bayesmixedfig %>% 
              add_trace(x = ~calibrationData()$Temperature, 
                        y = ~D47,
                        type = 'scatter', 
                        mode = 'markers', 
                        marker = list(color = 'black'),
                        opacity = 0.5,
                        name = 'Raw data',
                        text = as.character(calibrationData()$Sample.Name),
                        hovertemplate = paste(
                          "<b>Sample: %{text}</b><br><br>",
                          "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                          "Δ<sub>47</sub> (‰): %{y}<br>",
                          "Mineralogy: ", as.character(calibrationData()$Mineralogy),"<br>",
                          "Type: ", as.character(calibrationData()$Material),
                          "<extra></extra>")) %>%
              add_ribbons(data = bayeslmmincalciwitherror,
                          x = ~x,
                          y = ~median_est,
                          ymin = ~ci_lower_est,
                          ymax = ~ci_upper_est,
                          line = list(color = '#446455'),
                          fillcolor = '#446455',
                          opacity = 0.5,
                          name = '95% CI - with error',
                          hovertemplate = paste(
                            "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                            "Δ<sub>47</sub> (‰): %{y}<br>")) %>%
              add_lines(data = bayeslmmincalciwitherror,
                        x = ~x,
                        y = ~median_est,
                        name = 'Median estimate - with error',
                        line = list(color = "#446455", dash = 'dash'),
                        hovertemplate = paste(
                          "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                          "Δ<sub>47</sub> (‰): %{y}<br>"))
            bayesmixedfig <- bayesmixedfig %>% layout(title = '<b> Bayesian mixed model </b>',
                                                  legend=list(title=list(text='Legend')),
                                                  xaxis = list(title = 'Temperature (10<sup>6</sup>/T<sup>2</sup>)', hoverformat = '.1f'), 
                                                  yaxis = list(title = 'Δ<sub>47</sub> (‰)', hoverformat = '.3f'))
            
            return(bayesmixedfig)
          })
          
          print(noquote("Bayesian mixed model complete"))
          
          output$blinmwerr <- renderPrint(
            ddply(bayesmixedcals$BLMM_Measured_errors, .( material), 
                  function(x) cbind.data.frame('Median'=round(median(x$intercept),4),
                            'Lower 95% CI'=round(quantile(x$intercept, c(0.025)),4),
                            'Upper 95% CI'=round(quantile(x$intercept, c(0.975)),4),
                            'medianSlope'=round(median(x$slope),4),
                            'lwrSlope'=round(quantile(x$slope, c(0.025)),4),
                            'uprSlope'=round(quantile(x$slope, c(0.975)),4)
            ))
        ) 

      
        }
        
      })
      
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
                           type = 'scatter', 
                           mode = 'lines+markers', 
                           linetype = ~as.factor(Material), 
                           color = ~as.factor(Mineralogy),
                           colors = viridis_pal(option = "D", end = 0.9)(minlength),
                           opacity = 0.6,
                           error_y = ~list(array = ~D47error, color = '#000000'),
                           error_x = ~list(array = ~TempError, color = '#000000'),
                           text = as.character(calibrationData()$Sample.Name),
                           hovertemplate = paste(
                             "<b>Sample: %{text}</b><br><br>",
                             "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                             "Δ<sub>47</sub> (‰): %{y}<br>",
                             "Mineralogy: ", as.character(calibrationData()$Mineralogy),"<br>",
                             "Type: ", as.character(calibrationData()$Material),
                             "<extra></extra>"))
      rawcalfig <- rawcalfig %>% layout(title = '<b> Raw calibration data from user input </b>',
                                        legend=list(title=list(text='Material and mineralogy')),
                                        xaxis = list(title = 'Temperature (10<sup>6</sup>/T<sup>2</sup>)'), 
                                        yaxis = list(title = 'Δ<sub>47</sub> (‰)', hoverformat = '.3f'))
      
      return(rawcalfig)
    })
    }else{
      output$rawcaldata <- renderPlotly({
        rawcalfig <- plot_ly(calibrationData(), 
                             x = ~Temperature, 
                             y = ~D47, 
                             type = 'scatter', 
                             mode = 'lines+markers', 
                             linetype = ~as.factor(Material), 
                             color = ~as.factor(Material),
                             colors = viridis_pal(option = "D", end = 0.9)(length(unique(calibrationData()$Material))),
                             opacity = 0.6,
                             error_y = ~list(array = ~D47error, color = '#000000'),
                             error_x = ~list(array = ~TempError, color = '#000000'),
                             text = as.character(calibrationData()$Sample.Name),
                             hovertemplate = paste(
                               "<b>Sample: %{text}</b><br><br>",
                               "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                               "Δ<sub>47</sub> (‰): %{y}<br>",
                               "Type: ", as.character(calibrationData()$Material),
                               "<extra></extra>"))
        rawcalfig <- rawcalfig %>% layout(title = '<b> Raw calibration data from user input </b>',
                                          legend=list(title=list(text='Material')),
                                          xaxis = list(title = 'Temperature (10<sup>6</sup>/T<sup>2</sup>)'), 
                                          yaxis = list(title = 'Δ<sub>47</sub> (‰)', hoverformat = '.3f'))
        
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
          "Unique samples" = length(unique(reconstructionData()$Sample.Name)),
          "Total replicates" = sum(reconstructionData()$N),
          "Mineralogies" = length(unique(reconstructionData()$Mineralogy)),
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
      
      hasMaterial <- ifelse( is.na(reconstructionData()$Material), FALSE, TRUE )
      
      # Remove existing worksheets from wb2 on 'run' click, if any
      if("Linear w uncertainty" %in% names(wb2) == TRUE) 
      {removeWorksheet(wb2, "Linear w uncertainty") & removeWorksheet(wb2, "Linear w no uncertainty")}
      if("Inverse linear w uncertainty" %in% names(wb2) == TRUE) 
      {removeWorksheet(wb2, "Inverse linear w uncertainty") & removeWorksheet(wb2, "Inverse linear w no uncertainty")}
      if("York w uncertainty" %in% names(wb2) == TRUE) 
      {removeWorksheet(wb2, "York w uncertainty") & removeWorksheet(wb2, "York w no uncertainty")}
      if("Deming w uncertainty" %in% names(wb2) == TRUE) 
      {removeWorksheet(wb2, "Deming w uncertainty") & removeWorksheet(wb2, "Deming w no uncertainty")}
      if("Bayes w uncertainty and error" %in% names(wb2) == TRUE) 
      {removeWorksheet(wb2, "Bayes w uncertainty and error") & removeWorksheet(wb2, "Bayes w uncertainty no error") &
          removeWorksheet(wb2, "Bayes w error no uncertainty") & removeWorksheet(wb2, "Bayes no error no uncertainty")}
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
      
      if(input$confirm == FALSE) { print(noquote("Please confirm that your reference frames match")) }
      if(input$confirm == TRUE) {
        
        withProgress(message = 'Running selected reconstructions, please wait', {
          
          # Misc options
          
          if(input$bayesianPredictions == TRUE) {
            
            ##This function runs only Bayesian predictions
            ##(Only Bayesian simple linear with error for now)
            sink("bayespredictions.txt", type = "output")
            PipCriteria<-read.csv('Data/PipCriteria.csv')
            infTempBayesianC <<- clumpipe(calData=calData,
                                       PipCriteria=PipCriteria, 
                                       targetD47=recData$D47, 
                                       error_targetD47=recData$D47error, 
                                       materials = as.numeric(as.factor(ifelse(is.na(recData$Material), 1,recData$Material))),
                                       nrep=100,
                                       hasMaterial = T,
                                       BayesianOnly=T,
                                       generations=ngenerationsBayesianPredictions)
            sink()
            infTempBayesian_werrors<-infTempBayesianC[[1]][,-1]
            infTempBayesian_werrors$Tc<-sqrt(10^6/infTempBayesian_werrors$Tc)-273.15
            colnames(infTempBayesian_werrors)[c(4:5)]<-c('upr', 'lwr')
            infTempBayesian_werrors$lwr<-sqrt(10^6/infTempBayesian_werrors$lwr)-273.15
            infTempBayesian_werrors$upr<-sqrt(10^6/infTempBayesian_werrors$upr)-273.15
            infTempBayesian_werrors<-infTempBayesian_werrors[,c("D47","D47error", "Tc","lwr","upr")]

          
          df0<-infTempBayesian_werrors
          names(df0) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
          rownames(df0) <- NULL
          
          output$BpredictionsErrors <- renderTable({
            
            df0$`Δ47 (‰)` <- formatC(df0$`Δ47 (‰)`, digits = 3, format = "f")
            df0$`Δ47 (‰) error` <- formatC(df0$`Δ47 (‰) error`, digits = 3, format = "f")
            df0$`Temperature (°C)` <- formatC(df0$`Temperature (°C)`, digits = 1, format = "f")
            df0$`Lower 95% CI` <- formatC(df0$`Lower 95% CI`, digits = 1, format = "f")
            df0$`Upper 95% CI` <- formatC(df0$`Upper 95% CI`, digits = 1, format = "f")
            head(df0)
          },
            caption = "Bayesian predictions (BLM_errors)",
            caption.placement = getOption("xtable.caption.placement", "top"),
          rownames = FALSE,
          spacing = "m",
          align = "c"
            
          )
          
          addWorksheet(wb2, "Bayesian linear model, errors") # Add a blank sheet
          writeData(wb2, sheet = "Bayesian linear model, errors", df0)
          
          ##Without errors
          infTempBayesian<-infTempBayesianC[[2]][,-1]
          infTempBayesian$Tc<-sqrt(10^6/infTempBayesian$Tc)-273.15
          colnames(infTempBayesian)[c(4:5)]<-c('upr', 'lwr')
          infTempBayesian$lwr<-sqrt(10^6/infTempBayesian$lwr)-273.15
          infTempBayesian$upr<-sqrt(10^6/infTempBayesian$upr)-273.15
          infTempBayesian<-infTempBayesian[,c("D47","D47error", "Tc","lwr","upr")]

          
          df0.1<-infTempBayesian
          names(df0.1) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
          rownames(df0.1) <- NULL
          
          output$Bpredictions <- renderTable({
            
            df0.1$`Δ47 (‰)` <- formatC(df0.1$`Δ47 (‰)`, digits = 3, format = "f")
            df0.1$`Δ47 (‰) error` <- formatC(df0.1$`Δ47 (‰) error`, digits = 3, format = "f")
            df0.1$`Temperature (°C)` <- formatC(df0.1$`Temperature (°C)`, digits = 1, format = "f")
            df0.1$`Lower 95% CI` <- formatC(df0.1$`Lower 95% CI`, digits = 1, format = "f")
            df0.1$`Upper 95% CI` <- formatC(df0.1$`Upper 95% CI`, digits = 1, format = "f")
            head(df0.1)
          },
            caption = "Bayesian predictions (BLM without errors)",
            caption.placement = getOption("xtable.caption.placement", "top"),
          rownames = FALSE,
          spacing = "m",
          align = "c"
            
          )
          
          addWorksheet(wb2, "Bayesian linear model") # Add a blank sheet
          writeData(wb2, sheet = "Bayesian linear model", df0.1)
          
          ##BLMM
          infTempBayesianBLMM<-infTempBayesianC[[3]][,-1]
          infTempBayesianBLMM$Tc<-sqrt(10^6/infTempBayesianBLMM$Tc)-273.15
          colnames(infTempBayesianBLMM)[c(4:5)]<-c('upr', 'lwr')
          infTempBayesianBLMM$lwr<-sqrt(10^6/infTempBayesianBLMM$lwr)-273.15
          infTempBayesianBLMM$upr<-sqrt(10^6/infTempBayesianBLMM$upr)-273.15
          infTempBayesianBLMM<-infTempBayesianBLMM[,c("D47","D47error", "Tc","lwr","upr")]

          
          df0.2<-infTempBayesianBLMM
          names(df0.2) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
          rownames(df0.2) <- NULL
          
          output$BpredictionsBLMM <- renderTable({
            
            df0.2$`Δ47 (‰)` <- formatC(df0.2$`Δ47 (‰)`, digits = 3, format = "f")
            df0.2$`Δ47 (‰) error` <- formatC(df0.2$`Δ47 (‰) error`, digits = 3, format = "f")
            df0.2$`Temperature (°C)` <- formatC(df0.2$`Temperature (°C)`, digits = 1, format = "f")
            df0.2$`Lower 95% CI` <- formatC(df0.2$`Lower 95% CI`, digits = 1, format = "f")
            df0.2$`Upper 95% CI` <- formatC(df0.2$`Upper 95% CI`, digits = 1, format = "f")
            head(df0.2)
          },
            caption = "Bayesian predictions under a Bayesian linear mixed model",
            caption.placement = getOption("xtable.caption.placement", "top"),
          rownames = FALSE,
          spacing = "m",
          align = "c"
            
          )
          
          addWorksheet(wb2, "Bayesian linear mixed model") # Add a blank sheet
          writeData(wb2, sheet = "Bayesian linear mixed model", df0.2)
          
        }

          # Run prediction function
          if( !is.null(lmcals) ) {
            
            lmrec <<- predictTcNonBayes(data=cbind(recData$D47,recData$D47error), 
                                        slope=median(lmcals$slope), 
                                        slpcnf=CItoSE(quantile(lmcals$slope, 0.975), quantile(lmcals$slope, 0.025)), 
                                        intercept=median(lmcals$intercept), 
                                        intcnf=CItoSE(quantile(lmcals$intercept, 0.975), quantile(lmcals$intercept, 0.025)))
            
            lmrecwun <- lmrec[lmrec$type == "Parameter uncertainty", -1]
            df1 <- lmrecwun
            
            names(df1) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
            rownames(df1) <- NULL
            
            lmrecwoun <- lmrec[lmrec$type == "No parameter uncertainty", -1]

            df2<-lmrecwoun
            names(df2) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
            rownames(df2) <- NULL
            
            output$lmrecswun <- renderTable({
              
              df1$`Δ47 (‰)` <- formatC(df1$`Δ47 (‰)`, digits = 3, format = "f")
              df1$`Δ47 (‰) error` <- formatC(df1$`Δ47 (‰) error`, digits = 3, format = "f")
              df1$`Temperature (°C)` <- formatC(df1$`Temperature (°C)`, digits = 1, format = "f")
              df1$`Lower 95% CI` <- formatC(df1$`Lower 95% CI`, digits = 1, format = "f")
              df1$`Upper 95% CI` <- formatC(df1$`Upper 95% CI`, digits = 1, format = "f")
              head(df1)
            },
              caption = "Linear model with parameter uncertainty",
              caption.placement = getOption("xtable.caption.placement", "top"),
              rownames = FALSE,
              spacing = "m",
              align = "c"
            )
            
            output$lmrecswoun <- renderTable({
              
              df2$`Δ47 (‰)` <- formatC(df2$`Δ47 (‰)`, digits = 3, format = "f")
              df2$`Δ47 (‰) error` <- formatC(df2$`Δ47 (‰) error`, digits = 3, format = "f")
              df2$`Temperature (°C)` <- formatC(df2$`Temperature (°C)`, digits = 1, format = "f")
              df2$`Lower 95% CI` <- formatC(df2$`Lower 95% CI`, digits = 1, format = "f")
              df2$`Upper 95% CI` <- formatC(df2$`Upper 95% CI`, digits = 1, format = "f")
              head(df2)
              },
              caption = "Linear model without parameter uncertainty",
              caption.placement = getOption("xtable.caption.placement", "top"),
              rownames = FALSE,
              spacing = "m",
              align = "c"
              
            )
            
            addWorksheet(wb2, "Linear w uncertainty") # Add a blank sheet
            addWorksheet(wb2, "Linear w no uncertainty") # Add a blank sheet
            
            lmrecwun2 <- lmrec[lmrec$type == "Parameter uncertainty", -1]
            names(lmrecwun2) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
            
            lmrecwoun2 <- lmrec[lmrec$type == "No parameter uncertainty", -1]
            names(lmrecwoun2) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
            
            writeData(wb2, sheet = "Linear w uncertainty", lmrecwun2) # Write reconstruction data
            writeData(wb2, sheet = "Linear w no uncertainty", lmrecwoun2) # Write reconstruction data
            
            print(noquote("Linear reconstruction complete"))
          }
          
          #Inverse weighted linear model 
          if( !is.null(lminversecals) ) {
            
            lminverserec <<- predictTcNonBayes(data=cbind(recData$D47,recData$D47error), 
                                               slope=median(lminversecals$slope), 
                                               slpcnf=CItoSE(quantile(lminversecals$slope, 0.975), quantile(lminversecals$slope, 0.025)), 
                                               intercept=median(lminversecals$intercept), 
                                               intcnf=CItoSE(quantile(lminversecals$intercept, 0.975), quantile(lminversecals$intercept, 0.025)))
            
            lminverserecwun <- lminverserec[lminverserec$type == "Parameter uncertainty", -1]

            df3<-lminverserecwun
            names(df3) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
            rownames(df3) <- NULL
            
            lminverserecwoun <- lminverserec[lminverserec$type == "No parameter uncertainty", -1]

            df4<-lminverserecwoun
            names(df4) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
            rownames(df4) <- NULL
            
            output$lminverserecswun <- renderTable({
              
              df3$`Δ47 (‰)` <- formatC(df3$`Δ47 (‰)`, digits = 3, format = "f")
              df3$`Δ47 (‰) error` <- formatC(df3$`Δ47 (‰) error`, digits = 3, format = "f")
              df3$`Temperature (°C)` <- formatC(df3$`Temperature (°C)`, digits = 1, format = "f")
              df3$`Lower 95% CI` <- formatC(df3$`Lower 95% CI`, digits = 1, format = "f")
              df3$`Upper 95% CI` <- formatC(df3$`Upper 95% CI`, digits = 1, format = "f")
              head(df3)
            },
              caption = "Inverse weighted linear model with parameter uncertainty",
              caption.placement = getOption("xtable.caption.placement", "top"),
            rownames = FALSE,
            spacing = "m",
            align = "c"
              
            )
            
            output$lminverserecswoun <- renderTable({
              
              df4$`Δ47 (‰)` <- formatC(df4$`Δ47 (‰)`, digits = 3, format = "f")
              df4$`Δ47 (‰) error` <- formatC(df4$`Δ47 (‰) error`, digits = 3, format = "f")
              df4$`Temperature (°C)` <- formatC(df4$`Temperature (°C)`, digits = 1, format = "f")
              df4$`Lower 95% CI` <- formatC(df4$`Lower 95% CI`, digits = 1, format = "f")
              df4$`Upper 95% CI` <- formatC(df4$`Upper 95% CI`, digits = 1, format = "f")
              head(df4)
            },
              caption = "Inverse weighted linear model without parameter uncertainty",
              caption.placement = getOption("xtable.caption.placement", "top"),
            rownames = FALSE,
            spacing = "m",
            align = "c"
              
            )
            
            addWorksheet(wb2, "Inverse linear w uncertainty") # Add a blank sheet
            addWorksheet(wb2, "Inverse linear w no uncertainty") # Add a blank sheet
            
            lminverserecwun <- lminverserec[lminverserec$type == "Parameter uncertainty", -1]
            names(lminverserecwun) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
            
            lminverserecwoun2 <- lminverserec[lminverserec$type == "No parameter uncertainty", -1]
            names(lminverserecwoun2) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
            
            writeData(wb2, sheet = "Inverse linear w uncertainty", lminverserecwun) # Write reconstruction data
            writeData(wb2, sheet = "Inverse linear w no uncertainty", lminverserecwoun2) # Write reconstruction data
            
            print(noquote("Inverse weighted linear reconstruction complete"))
            
          }
          
          # York regression
          if( !is.null(yorkcals) ) {
            
            yorkrec <<- predictTcNonBayes(data=cbind(recData$D47,recData$D47error), 
                                          slope=median(yorkcals$slope), 
                                          slpcnf=CItoSE(quantile(yorkcals$slope, 0.975), quantile(yorkcals$slope, 0.025)), 
                                          intercept=median(yorkcals$intercept), 
                                          intcnf=CItoSE(quantile(yorkcals$intercept, 0.975), quantile(yorkcals$intercept, 0.025)))
            
            yorkrecwun <- yorkrec[yorkrec$type == "Parameter uncertainty", -1]
            
            df5<-yorkrecwun
            names(df5) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
            rownames(df5) <- NULL
            
            yorkrecwoun <- yorkrec[yorkrec$type == "No parameter uncertainty", -1]

            df6<-yorkrecwoun
            names(df6) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
            rownames(df6) <- NULL
            
            output$yorkrecswun <- renderTable({
              
              df5$`Δ47 (‰)` <- formatC(df5$`Δ47 (‰)`, digits = 3, format = "f")
              df5$`Δ47 (‰) error` <- formatC(df5$`Δ47 (‰) error`, digits = 3, format = "f")
              df5$`Temperature (°C)` <- formatC(df5$`Temperature (°C)`, digits = 1, format = "f")
              df5$`Lower 95% CI` <- formatC(df5$`Lower 95% CI`, digits = 1, format = "f")
              df5$`Upper 95% CI` <- formatC(df5$`Upper 95% CI`, digits = 1, format = "f")
              head(df5)
            },
              caption = "York regression with parameter uncertainty",
              caption.placement = getOption("xtable.caption.placement", "top"),
            rownames = FALSE,
            spacing = "m",
            align = "c"
              
            )
            
            output$yorkrecswoun <- renderTable({
              
              df6$`Δ47 (‰)` <- formatC(df6$`Δ47 (‰)`, digits = 3, format = "f")
              df6$`Δ47 (‰) error` <- formatC(df6$`Δ47 (‰) error`, digits = 3, format = "f")
              df6$`Temperature (°C)` <- formatC(df6$`Temperature (°C)`, digits = 1, format = "f")
              df6$`Lower 95% CI` <- formatC(df6$`Lower 95% CI`, digits = 1, format = "f")
              df6$`Upper 95% CI` <- formatC(df6$`Upper 95% CI`, digits = 1, format = "f")
              head(df6)
            },
              caption = "York regression without parameter uncertainty",
              caption.placement = getOption("xtable.caption.placement", "top"),
            rownames = FALSE,
            spacing = "m",
            align = "c"
              
            )
            
            addWorksheet(wb2, "York w uncertainty") # Add a blank sheet
            addWorksheet(wb2, "York w no uncertainty") # Add a blank sheet
            
            yorkrecwun2 <- yorkrec[yorkrec$type == "Parameter uncertainty", -1]
            names(yorkrecwun2) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
            
            yorkrecwoun2 <- yorkrec[yorkrec$type == "No parameter uncertainty", -1]
            names(yorkrecwoun2) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
            
            writeData(wb2, sheet = "York w uncertainty", yorkrecwun2) # Write reconstruction data
            writeData(wb2, sheet = "York w no uncertainty", yorkrecwoun2) # Write reconstruction data
            
            print(noquote("York reconstruction complete"))
            
          }
          
          # Deming regression
          if( !is.null(demingcals) ) {
            
            demingrec <<- predictTcNonBayes(data=cbind(recData$D47,recData$D47error), 
                                            slope=median(demingcals$slope), 
                                            slpcnf=CItoSE(quantile(demingcals$slope, 0.975), quantile(demingcals$slope, 0.025)), 
                                            intercept=median(demingcals$intercept), 
                                            intcnf=CItoSE(quantile(demingcals$intercept, 0.975), quantile(demingcals$intercept, 0.025)))
            
            demingrecwun <- demingrec[demingrec$type == "Parameter uncertainty", -1]

            df7<-demingrecwun
            names(df7) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
            rownames(df7) <- NULL
            
            demingrecwoun <- demingrec[demingrec$type == "No parameter uncertainty", -1]

            df8<-demingrecwoun
            names(df8) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
            rownames(df8) <- NULL
            
            output$demingrecswun <- renderTable({
              
              df7$`Δ47 (‰)` <- formatC(df7$`Δ47 (‰)`, digits = 3, format = "f")
              df7$`Δ47 (‰) error` <- formatC(df7$`Δ47 (‰) error`, digits = 3, format = "f")
              df7$`Temperature (°C)` <- formatC(df7$`Temperature (°C)`, digits = 1, format = "f")
              df7$`Lower 95% CI` <- formatC(df7$`Lower 95% CI`, digits = 1, format = "f")
              df7$`Upper 95% CI` <- formatC(df7$`Upper 95% CI`, digits = 1, format = "f")
              head(df7)
            },
              caption = "Deming regression with parameter uncertainty",
              caption.placement = getOption("xtable.caption.placement", "top"),
            rownames = FALSE,
            spacing = "m",
            align = "c"
              
            )
            
            output$demingrecswoun <- renderTable({
              
              df8$`Δ47 (‰)` <- formatC(df8$`Δ47 (‰)`, digits = 3, format = "f")
              df8$`Δ47 (‰) error` <- formatC(df8$`Δ47 (‰) error`, digits = 3, format = "f")
              df8$`Temperature (°C)` <- formatC(df8$`Temperature (°C)`, digits = 1, format = "f")
              df8$`Lower 95% CI` <- formatC(df8$`Lower 95% CI`, digits = 1, format = "f")
              df8$`Upper 95% CI` <- formatC(df8$`Upper 95% CI`, digits = 1, format = "f")
              head(df8)
            },
              caption = "Deming regression without parameter uncertainty",
              caption.placement = getOption("xtable.caption.placement", "top"),
            rownames = FALSE,
            spacing = "m",
            align = "c"
              
            )
            
            addWorksheet(wb2, "Deming w uncertainty") # Add a blank sheet
            addWorksheet(wb2, "Deming w no uncertainty") # Add a blank sheet
            
            demingrecwun2 <- demingrec[demingrec$type == "Parameter uncertainty", -1]
            names(demingrecwun2) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
            
            demingrecwoun2 <- demingrec[demingrec$type == "No parameter uncertainty", -1]
            names(demingrecwoun2) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
            
            writeData(wb2, sheet = "Deming w uncertainty", demingrecwun2) # Write reconstruction data
            writeData(wb2, sheet = "Deming w no uncertainty", demingrecwoun2) # Write reconstruction data
            
            print(noquote("Deming reconstruction complete"))
            
          }
          
          if( !is.null(bayeslincals) ) {
            
            sink(file = "Bayesrectext.txt", type = "output")
            
            bayesrec <<- predictTcNonBayes(data=cbind(recData$D47,recData$D47error), 
                                           slope=median(bayeslincals[[1]]$slope), 
                                           slpcnf=CItoSE(quantile(bayeslincals[[1]]$slope, 0.975), quantile(bayeslincals[[1]]$slope, 0.025)), 
                                           intercept=median(bayeslincals[[1]]$intercept), 
                                           intcnf=CItoSE(quantile(bayeslincals[[1]]$intercept, 0.975), quantile(bayeslincals[[1]]$intercept, 0.025)))
            
            bayesrecNE <<- predictTcNonBayes(data=cbind(recData$D47,recData$D47error), 
                                             slope=median(bayeslincals[[1]]$slope), 
                                             slpcnf=CItoSE(quantile(bayeslincals[[2]]$slope, 0.975), quantile(bayeslincals[[2]]$slope, 0.025)), 
                                             intercept=median(bayeslincals[[2]]$intercept), 
                                             intcnf=CItoSE(quantile(bayeslincals[[2]]$intercept, 0.975), quantile(bayeslincals[[2]]$intercept, 0.025)))
            sink()
            
            # With parameter uncertainty and errors
            bayesrecwunerr <- bayesrec[bayesrec$type == "Parameter uncertainty" , -c(1, 7, 8)]

            df9<-bayesrecwunerr
            names(df9) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
            rownames(df9) <- NULL
            
            # With parameter uncertainty and no errors
            bayesrecwunnoerr <- bayesrecNE[bayesrecNE$type == "Parameter uncertainty" , -c(1, 7, 8)]

            df10<-bayesrecwunnoerr
            names(df10) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
            rownames(df10) <- NULL
            
            # With errors and no parameter uncertainty
            bayesrecwounerr <- bayesrec[bayesrec$type == "No parameter uncertainty" , -c(1, 7, 8)]

            df11<-bayesrecwounerr
            names(df11) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
            rownames(df11) <- NULL
            
            # With no errors and no parameter uncertainty
            bayesrecwounnoerr <- bayesrecNE[bayesrec$type == "No parameter uncertainty" , -c(1, 7, 8)]

            df12<-bayesrecwounnoerr
            names(df12) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
            rownames(df12) <- NULL
            
            
            output$bayesrecswunerr <- renderTable({
              
              df9$`Δ47 (‰)` <- formatC(df9$`Δ47 (‰)`, digits = 3, format = "f")
              df9$`Δ47 (‰) error` <- formatC(df9$`Δ47 (‰) error`, digits = 3, format = "f")
              df9$`Temperature (°C)` <- formatC(df9$`Temperature (°C)`, digits = 1, format = "f")
              df9$`Lower 95% CI` <- formatC(df9$`Lower 95% CI`, digits = 1, format = "f")
              df9$`Upper 95% CI` <- formatC(df9$`Upper 95% CI`, digits = 1, format = "f")
              head(df9)
            },
              caption = "Bayesian regression with parameter uncertainty and errors",
              caption.placement = getOption("xtable.caption.placement", "top"),
            rownames = FALSE,
            spacing = "m",
            align = "c"
              
            )
            
            output$bayesrecswunnoerr <- renderTable({
              
              df10$`Δ47 (‰)` <- formatC(df10$`Δ47 (‰)`, digits = 3, format = "f")
              df10$`Δ47 (‰) error` <- formatC(df10$`Δ47 (‰) error`, digits = 3, format = "f")
              df10$`Temperature (°C)` <- formatC(df10$`Temperature (°C)`, digits = 1, format = "f")
              df10$`Lower 95% CI` <- formatC(df10$`Lower 95% CI`, digits = 1, format = "f")
              df10$`Upper 95% CI` <- formatC(df10$`Upper 95% CI`, digits = 1, format = "f")
              head(df10)
            },
              caption = "Bayesian regression with parameter uncertainty and no errors",
              caption.placement = getOption("xtable.caption.placement", "top"),
            rownames = FALSE,
            spacing = "m",
            align = "c"
              
            )
            
            output$bayesrecswounerr <- renderTable({
              
              df11$`Δ47 (‰)` <- formatC(df11$`Δ47 (‰)`, digits = 3, format = "f")
              df11$`Δ47 (‰) error` <- formatC(df11$`Δ47 (‰) error`, digits = 3, format = "f")
              df11$`Temperature (°C)` <- formatC(df11$`Temperature (°C)`, digits = 1, format = "f")
              df11$`Lower 95% CI` <- formatC(df11$`Lower 95% CI`, digits = 1, format = "f")
              df11$`Upper 95% CI` <- formatC(df11$`Upper 95% CI`, digits = 1, format = "f")
              head(df11)
            },
              caption = "Bayesian regression with errors and no parameter uncertainty",
              caption.placement = getOption("xtable.caption.placement", "top"),
            rownames = FALSE,
            spacing = "m",
            align = "c"
              
            )
            
            output$bayesrecswounnoerr <- renderTable({
              
              df12$`Δ47 (‰)` <- formatC(df12$`Δ47 (‰)`, digits = 3, format = "f")
              df12$`Δ47 (‰) error` <- formatC(df12$`Δ47 (‰) error`, digits = 3, format = "f")
              df12$`Temperature (°C)` <- formatC(df12$`Temperature (°C)`, digits = 1, format = "f")
              df12$`Lower 95% CI` <- formatC(df12$`Lower 95% CI`, digits = 1, format = "f")
              df12$`Upper 95% CI` <- formatC(df12$`Upper 95% CI`, digits = 1, format = "f")
              head(df12)
            },
              caption = "Bayesian regression without errors or parameter uncertainty",
              caption.placement = getOption("xtable.caption.placement", "top"),
            rownames = FALSE,
            spacing = "m",
            align = "c"
              
            )
            
            addWorksheet(wb2, "Bayes w uncertainty and error") # Add a blank sheet
            addWorksheet(wb2, "Bayes w uncertainty no error") # Add a blank sheet 
            addWorksheet(wb2, "Bayes w error no uncertainty") # Add a blank sheet
            addWorksheet(wb2, "Bayes no error no uncertainty") # Add a blank sheet
            
            # With parameter uncertainty and errors
            bayesrecwunerr2 <- bayesrec[bayesrec$type == "Parameter uncertainty" , -c(1, 7, 8)]
            names(bayesrecwunerr2) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
            
            # With parameter uncertainty and no errors
            bayesrecwunnoerr2 <- bayesrecNE[bayesrecNE$type == "Parameter uncertainty" , -c(1, 7, 8)]
            names(bayesrecwunnoerr2) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
            
            # With error and no parameter uncertainty
            bayesrecwounerr2 <- bayesrec[bayesrec$type == "No parameter uncertainty" , -c(1, 7, 8)]
            names(bayesrecwounerr2) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
            
            # No errors no parameter uncertainty
            bayesrecwounnoerr2 <- bayesrecNE[bayesrecNE$type == "No parameter uncertainty" , -c(1, 7, 8)]
            names(bayesrecwounnoerr2) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
            
            writeData(wb2, sheet = "Bayes w uncertainty and error", bayesrecwunerr2) # Write reconstruction data
            writeData(wb2, sheet = "Bayes w uncertainty no error", bayesrecwunnoerr2) # Write reconstruction data
            writeData(wb2, sheet = "Bayes w error no uncertainty", bayesrecwounerr2) # Write reconstruction data
            writeData(wb2, sheet = "Bayes no error no uncertainty", bayesrecwounnoerr2) # Write reconstruction data
            
            print(noquote("Bayesian reconstruction complete"))
            
          }
          
          if( !is.null(bayesmixedcals) ) {

              sink(file = "Bayesmixedrectext.txt", type = "output")
              
              bayesmixedrec <<- predictTcNonBayes(data=cbind(recData$D47,recData$D47error), 
                                             slope=median(bayesmixedcals[[1]]$slope), 
                                             slpcnf=CItoSE(quantile(bayesmixedcals[[1]]$slope, 0.975), quantile(bayesmixedcals[[1]]$slope, 0.025)), 
                                             intercept=median(bayesmixedcals[[1]]$intercept), 
                                             intcnf=CItoSE(quantile(bayesmixedcals[[1]]$intercept, 0.975), quantile(bayesmixedcals[[1]]$intercept, 0.025)))

              sink()
              
              addWorksheet(wb2, "Bayesian mixed model") # Add a blank sheet
              
              df13<-bayesmixedrec
              names(df13) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "Lower 95% CI", "Upper 95% CI")
              rownames(df13) <- NULL
              
              writeData(wb2, "Bayesian mixed model", bayesmixedrec)
              
              output$bayesrecsmixed <- renderTable({
                
                df13$`Δ47 (‰)` <- formatC(df13$`Δ47 (‰)`, digits = 3, format = "f")
                df13$`Δ47 (‰) error` <- formatC(df13$`Δ47 (‰) error`, digits = 3, format = "f")
                df13$`Temperature (°C)` <- formatC(df13$`Temperature (°C)`, digits = 1, format = "f")
                df13$`Lower 95% CI` <- formatC(df13$`Lower 95% CI`, digits = 1, format = "f")
                df13$`Upper 95% CI` <- formatC(df13$`Upper 95% CI`, digits = 1, format = "f")
                head(df13)
              },
                caption = "Bayesian mixed model with parameter uncertainty and errors",
                caption.placement = getOption("xtable.caption.placement", "top"),
              rownames = FALSE,
              spacing = "m",
              align = "c"
                
              )
          }
        })
      } 
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
  
}
  