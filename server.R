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

  #For parameter estimates
  if(exists("wb")) rm(wb) # Delete any existing workbook in preparation for new results
  wb <- createWorkbook("calibration output") # Prepare a workbook for calibration outputs
  
  #For convergence
  if(exists("wb3")) rm(wb3) # Delete any existing workbook in preparation for new results
  wb3 <- createWorkbook("Bayesian output") # Prepare a workbook for calibration outputs
  
  
  observeEvent(calibrationData(),{
  output$myList <-  renderUI({
    numericInput("samples", min = 3, max = nrow(calibrationData()), 
                label = paste0("Number of observations per bootstrap sample, (max. recommended: ", nrow(calibrationData()) ,")" ), 
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
          #"Mineralogies" = length(unique(calibrationData()$Mineralogy)),
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
    
    priors <<- input$priors
    replicates <<- input$replication
    ngenerationsBayes <<- input$generations
    ngenerationsBayesianPredictions <<- input$generations

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
    
    
    ##Also for the Bayesian sheet
    if("Bayesian model no errors" %in% names(wb3) == TRUE) 
    {removeWorksheet(wb3, "Bayesian model no errors")}
    if("Bayesian model with errors" %in% names(wb3) == TRUE) 
    {removeWorksheet(wb3, "Bayesian model with errors") }
    if("Bayesian mixed w errors" %in% names(wb3) == TRUE)
    {removeWorksheet(wb3, "Bayesian mixed w errors")}
    
    
    calData <<- NULL
    calData <<- calibrationData()

    samples <<- if(is.null(input$samples)){nrow(calData) }else{input$samples}
    samples <<- ifelse(samples==2, 3, samples)
    multicore <<- input$simulateDeming
      
    # Recode NA or 0 error values to dummy value
    calData$D47error[calData$D47error == 0] <<- 0.000001
    calData$TempError[calData$TempError == 0] <<- 0.000001
    calData$D47error[is.na(calData$D47error)] <<- 0.000001
    calData$TempError[is.na(calData$TempError)] <<- 0.000001
    calData$Material[is.na(calData$Material)] <<- 1
    
    
    ##Limits of the CI
    
    minLim <- ifelse(input$range[1]==0, min(calData$Temperature),input$range[1])
    maxLim <- ifelse(input$range[2]==0, max(calData$Temperature),input$range[2])
    
    #minLim <- ifelse(input$MinLim==0, min(calData$Temperature),input$MinLim)
    #maxLim <- ifelse(input$MaxLim==0, max(calData$Temperature),input$MaxLim)
    
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
          lmcals <<- simulateLM_measured(calData, replicates = replicates, samples = samples)
          sink()
          
          lmci <<- RegressionSingleCI(data = lmcals, from = minLim, to = maxLim)
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
          
          
          
          cat("Linear regression complete \n *R^2=", round(unlist(attributes(lmcals)$R2[1],4)),
              " (95% CI, ",round(unlist(attributes(lmcals)$R2[2],4)),"-",round(unlist(attributes(lmcals)$R2[3],4)),")"
          )
          

          output$lmcal <- renderPrint({
            do.call(rbind.data.frame,apply(lmcals, 2, function(x){
              cbind.data.frame(Median= round(median(x), 4), `SE`=round(sd(x)/sqrt(length(x)),7), `Lower 95% CI`=round(quantile(x, 0.025), 4),  `Upper 95% CI`=round(quantile(x, 0.975), 4))
            }))
          })
          
        }
        
        if(input$simulateLM_inverseweights != FALSE) {
          sink(file = "inverselinmodtext.txt", type = "output")
          lminversecals <<- simulateLM_inverseweights(calData, replicates = replicates, samples = samples)
          sink()
          
          lminverseci <- RegressionSingleCI(data = lminversecals, from = minLim, to = maxLim)
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
          
          cat("\nInverse regression complete \n *R^2=", round(unlist(attributes(lminversecals)$R2[1],4)),
              " (95% CI, ",round(unlist(attributes(lminversecals)$R2[2],4)),"-",round(unlist(attributes(lminversecals)$R2[3],4)),")"
          )
          
          output$lminversecal <- renderPrint({
            do.call(rbind.data.frame,apply(lminversecals, 2, function(x){
              cbind.data.frame(Median= round(median(x), 4),`SE`=round(sd(x)/sqrt(length(x)),7), `Lower 95% CI`=round(quantile(x, 0.025), 4),  `Upper 95% CI`=round(quantile(x, 0.975), 4))
            }))
          })
          
        }
        
        if(input$simulateYork_measured != FALSE) {
          sink(file = "yorkmodtext.txt", type = "output")
          yorkcals <<- simulateYork_measured(calData, replicates = replicates, samples = samples)
          sink()
          
          yorkci <- RegressionSingleCI(data = yorkcals, from = minLim, to = maxLim)
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
          
          cat("\nYork regression complete")
          output$york <- renderPrint({
            do.call(rbind.data.frame,apply(yorkcals, 2, function(x){
              cbind.data.frame(Median= round(median(x), 4),`SE`=round(sd(x)/sqrt(length(x)),7), `Lower 95% CI`=round(quantile(x, 0.025), 4),  `Upper 95% CI`=round(quantile(x, 0.975), 4))
            }))
          })
          
        }
        
        if(input$simulateDeming != FALSE) {
          sink(file = "demingmodtext.txt", type = "output")
          demingcals <<- simulateDeming(calData, replicates = replicates, samples = samples, multicore=multicore)
          sink()
          
          demingci <- RegressionSingleCI(data = demingcals, from = minLim, to = maxLim)
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
          
          cat("\nDeming regression complete")
          output$deming <- renderPrint({
            do.call(rbind.data.frame,apply(demingcals, 2, function(x){
              cbind.data.frame(Median= round(median(x), 4),`SE`=round(sd(x)/sqrt(length(x)),7), `Lower 95% CI`=round(quantile(x, 0.025), 4),  `Upper 95% CI`=round(quantile(x, 0.975), 4))
            }))
          })
          
        }
        
        #     checkboxInput("linear", "Linear model", FALSE),
        if(input$simulateBLM_measuredMaterial != FALSE) {
          sink(file = "Bayeslinmodtext.txt", type = "output")
          bayeslincals <<- simulateBLM_measuredMaterial(calData, 
                                                        replicates = replicates, 
                                                        isMixed=F, 
                                                        generations=ngenerationsBayes,
                                                        priors = priors, 
                                                        samples = samples,
                                                        multicore = multicore)
          sink()
          
          bayeslincinoerror <- RegressionSingleCI(data = bayeslincals$BLM_Measured_no_errors, from = minLim, to = maxLim)
          bayeslincalcinoerror <- as.data.frame(bayeslincinoerror)
          bayeslinciwitherror <- RegressionSingleCI(data = bayeslincals$BLM_Measured_errors, from = minLim, to = maxLim)
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
          
          
          ##For the Bayesian sheet
          conv_BLM <- do.call(rbind, lapply(seq_along(attr(bayeslincals, "Conv")), function(x){
            cbind(Replicate=x,attr(bayeslincals, "Conv")[[x]][[1]]) }))
          conv_BLM_errors <- do.call(rbind, lapply(seq_along(attr(bayeslincals, "Conv")), function(x){
            cbind(Replicate=x,attr(bayeslincals, "Conv")[[x]][[2]]) }))
          
          conv_BLM <- cbind.data.frame(parameter=row.names(conv_BLM),conv_BLM)
          conv_BLM_errors <- cbind.data.frame(parameter=row.names(conv_BLM_errors),conv_BLM_errors)
          
          
          addWorksheet(wb3, "Bayesian model no errors") # Add a blank sheet
          addWorksheet(wb3, "Bayesian model with errors") # Add a blank sheet
          writeData(wb3, sheet = "Bayesian model no errors", conv_BLM) # Write regression data
          writeData(wb3, sheet = "Bayesian model with errors", conv_BLM_errors) # Write regression data
          
          
          #print(noquote("Bayesian linear model complete"))
          
          cat( paste0("\nBayesian linear model complete \n *with errors \n   *R^2=", round(attr(bayeslincals,"R2s")[1,2],4),
                               " (95% CI, ",round(attr(bayeslincals,"R2s")[1,3],4),"-",round(attr(bayeslincals,"R2s")[1,4],4),")",
                      "\n   *DIC=", round(attr(bayeslincals,"DICs")[1,1],4),
                      " (95% CI, ",round(attr(bayeslincals,"DICs")[1,2],4),"-",round(attr(bayeslincals,"DICs")[1,2],4),")",
                      "\n *without errors\n   *R^2=", round(attr(bayeslincals,"R2s")[2,2],4),
                               " (95% CI, ",round(attr(bayeslincals,"R2s")[2,3],4),"-",round(attr(bayeslincals,"R2s")[2,4],4),")",
                      "\n   *DIC=", round(attr(bayeslincals,"DICs")[2,1],4),
                      " (95% CI, ",round(attr(bayeslincals,"DICs")[2,2],4),"-",round(attr(bayeslincals,"DICs")[2,2],4),")"
          )
          )
          
          output$blinnoerr <- renderPrint({
            
            do.call(rbind.data.frame,apply(bayeslincals$BLM_Measured_no_errors, 2, function(x){
              cbind.data.frame(Median= round(median(x), 4),`SE`=round(sd(x)/sqrt(length(x)),7), `Lower 95% CI`=round(quantile(x, 0.025), 4),  `Upper 95% CI`=round(quantile(x, 0.975), 4))
            }))
            
          })
          
          output$blinwerr <- renderPrint({
            
            do.call(rbind.data.frame,apply(bayeslincals$BLM_Measured_errors, 2, function(x){
              cbind.data.frame(Median= round(median(x), 4),`SE`=round(sd(x)/sqrt(length(x)),7), `Lower 95% CI`=round(quantile(x, 0.025), 4),  `Upper 95% CI`=round(quantile(x, 0.975), 4))
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
          bayesmixedcals <<- simulateBLM_measuredMaterial(data=calData, 
                                                         replicates = replicates, 
                                                         isMixed = T, 
                                                         generations=ngenerationsBayes,
                                                         priors = priors, 
                                                         samples = samples,
                                                         multicore = multicore)
          sink()
          
          bayeslmminciwitherror <- RegressionSingleCI(data = bayesmixedcals$BLMM_Measured_errors, from = minLim, to = maxLim)
          bayeslmmincalciwitherror <- as.data.frame(bayeslmminciwitherror)

          addWorksheet(wb, "Bayesian mixed w errors") # Add a blank sheet
          addWorksheet(wb, "Bayesian mixed w errors CI") # Add a blank sheet 
          
          bayeslmmincalciwitherror2 <- bayeslmmincalciwitherror
          names(bayeslmmincalciwitherror2) <- c("10^6/T^2", "D47_median_est",	"D47_ci_lower_est",	"D47_ci_upper_est")
          
          writeData(wb, sheet = "Bayesian mixed w errors", bayesmixedcals$BLMM_Measured_errors) # Write regression data
          writeData(wb, sheet = "Bayesian mixed w errors CI", bayeslmmincalciwitherror2)
          
          ##For the Bayesian sheet
          conv_BLMM <- do.call(rbind, lapply(seq_along(attr(bayesmixedcals, "Conv")), function(x){
            cbind(Replicate=x,attr(bayesmixedcals, "Conv")[[x]][[1]]) }))

          conv_BLMM <- cbind.data.frame(parameter=row.names(conv_BLMM),conv_BLMM)
          
          addWorksheet(wb3, "Bayesian mixed w errors") # Add a blank sheet
          writeData(wb3, sheet = "Bayesian mixed w errors", conv_BLMM) # Write regression data

          
          
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
          
          
          cat( paste0("\nBayesian mixed model complete \n *R^2 (conditional)=", round(attr(bayesmixedcals,"R2s")[3,3],4),
                      " (95% CI, ",round(attr(bayesmixedcals,"R2s")[3,4],4),"-",round(attr(bayesmixedcals,"R2s")[3,5],4),")",
                      "\n *R^2 (marginal)=", round(attr(bayesmixedcals,"R2s")[4,3],4),
                      " (95% CI, ",round(attr(bayesmixedcals,"R2s")[4,4],4),"-",round(attr(bayesmixedcals,"R2s")[4,5],4),")",
                      "\n *DIC=", round(attr(bayesmixedcals,"DIC")[3,1],4),
                      " (95% CI, ",round(attr(bayesmixedcals,"DIC")[3,2],4),"-",round(attr(bayesmixedcals,"DIC")[3,3],4),")"
          )
          )
          
          output$blinmwerr <- renderPrint(
            ddply(bayesmixedcals$BLMM_Measured_errors, .( material), 
                  function(x) cbind.data.frame('median (I)'=round(median(x$intercept),4),
                                               'SE (I)'=round(sd(x$slope)/sqrt(length(x$intercept)),7),
                            'Lower 95% CI (I)'=round(quantile(x$intercept, c(0.025)),4),
                            'Upper 95% CI (I)'=round(quantile(x$intercept, c(0.975)),4),
                            
                            'median (S)'=round(median(x$slope),4),
                            'SE (S)'=round(sd(x$slope)/sqrt(length(x$slope)),7),
                            'lwr (S)'=round(quantile(x$slope, c(0.025)),4),
                            'upr (S)'=round(quantile(x$slope, c(0.975)),4))
            )
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
  
  
  output$downloadBayesian <- downloadHandler(
    filename = function() { 
      paste("Bayesian_output_", Sys.time(), ".xlsx", sep="")
    },
    
    content = function(file) {
      saveWorkbook(wb3, file, overwrite = TRUE)
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
      
      if(input$confirm == FALSE) { print(noquote("Please confirm that your reference frames match")) }
      if(input$confirm == TRUE) {
        
        withProgress(message = 'Running selected reconstructions, please wait', {
          
          # Misc options
          
          classicPredictions <- input$classicPreds
          
          recData_byS <- recData
          
        
          if(input$simulateBLM_measuredMaterialRec == TRUE){
            
            if(is.null(bayeslincals) & classicPredictions) { print(noquote("Please run the calibration step for Bayesian linear models first")) }else{
              
            ##This function runs only Bayesian predictions
            ##(Only Bayesian simple linear with error for now)
            sink("bayespredictions.txt", type = "output")
            
            ##Linear models
            infTempBayesianCLinear <- if(classicPredictions){ 
              
             cpreds<- predictTcBayes(calibrationData=calData,
                                     data=cbind(recData_byS$D47,
                                                ifelse(recData_byS$D47error==0,0.00001,recData_byS$D47error),
                                                as.numeric(as.factor(ifelse(is.na(recData_byS$Material), 1,recData_byS$Material)))),
                                     generations=ngenerationsBayes,
                                     hasMaterial=F, bootDataset=T, onlyMedian=T, replicates = replicates, multicore = multicore, priors=priors,
                                     errorsD47=F)
             
            }else{
              
              cpreds<- predictTcBayes(calibrationData=calData,
                              data=cbind(recData_byS$D47,
                              ifelse(recData_byS$D47error==0,0.00001,recData_byS$D47error),
                              as.numeric(as.factor(ifelse(is.na(recData_byS$Material), 1,recData_byS$Material)))),
                              generations=ngenerationsBayes,
                              hasMaterial=F, bootDataset=T, onlyMedian=T, replicates = replicates, multicore = multicore, priors=priors,
                             errorsD47=T)

              
            }
            
          sink()
          infTempBayesian_werrors<- infTempBayesianCLinear[infTempBayesianCLinear[,1]=="BLM1_fit",]

          df0<-infTempBayesian_werrors[,-c(1,4)]
          names(df0) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "SE (1SD) Temperature (°C)")
          rownames(df0) <- NULL
          
          output$BpredictionsErrors <- renderTable({
            
            df0$`Δ47 (‰)` <- formatC(df0$`Δ47 (‰)`, digits = 3, format = "f")
            df0$`Δ47 (‰) error` <- formatC(df0$`Δ47 (‰) error`, digits = 4, format = "f")
            df0$`Temperature (°C)` <- formatC(df0$`Temperature (°C)`, digits = 1, format = "f")
            df0$`SE (1SD) Temperature (°C)` <- formatC(df0$`SE (1SD) Temperature (°C)`, digits = 7, format = "f")
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
          infTempBayesian<-infTempBayesianCLinear[infTempBayesianCLinear[,1]=="BLM1_fit_NoErrors",]
          df0.1<-infTempBayesian[,-c(1,4)]
          names(df0.1) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "SE (1SD) Temperature (°C)")
          rownames(df0.1) <- NULL
          
          output$Bpredictions <- renderTable({
            
            df0.1$`Δ47 (‰)` <- formatC(df0.1$`Δ47 (‰)`, digits = 3, format = "f")
            df0.1$`Δ47 (‰) error` <- formatC(df0.1$`Δ47 (‰) error`, digits = 4, format = "f")
            df0.1$`Temperature (°C)` <- formatC(df0.1$`Temperature (°C)`, digits = 1, format = "f")
            df0.1$`SE (1SD) Temperature (°C)` <- formatC(df0.1$`SE (1SD) Temperature (°C)`, digits = 7, format = "f")
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
          
          }
          }
          
          
          # Run prediction function
          # Need to use a sample-based dataset
          
          
          if( input$simulateLM_measuredRec == TRUE ) {
            
            if(is.null(lmcals) & classicPredictions) { print(noquote("Please run the calibration step for the linear model first")) }else{
              
            
            calData$T2 <<- calData$Temperature
            
            lmrec <<-  do.call(rbind,lapply(1:nrow(recData_byS), function(x){
                a <- predictTclassic(calData, targety=recData_byS$D47[x], model='lm', replicates=replicates, bootDataset=!classicPredictions)
                b <- predictTclassic(calData, targety=recData_byS$D47[x]+recData_byS$D47error[x], model='lm', replicates=replicates, bootDataset=!classicPredictions)
                cbind.data.frame("D47"=recData_byS$D47[x],'D47se'=recData_byS$D47error[x], "Tc"=a$temp, "se"=a$temp-b$temp)
              } ))
            
            
            df1 <- lmrec
            
            
            names(df1) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "SE (1SD) Temperature (°C)")
            rownames(df1) <- NULL
            
           
            output$lmrecswun <- renderTable({
              
              df1$`Δ47 (‰)` <- formatC(df1$`Δ47 (‰)`, digits = 3, format = "f")
              df1$`Δ47 (‰) error` <- formatC(df1$`Δ47 (‰) error`, digits = 4, format = "f")
              df1$`Temperature (°C)` <- formatC(df1$`Temperature (°C)`, digits = 1, format = "f")
              df1$`SE (1SD) Temperature (°C)` <- formatC(df1$`SE (1SD) Temperature (°C)`, digits = 7, format = "f")
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
          if( input$simulateLM_inverseweightsRec == TRUE ) {
            
            if(is.null(lminversecals) & classicPredictions) { print(noquote("Please run the calibration step for the weighted linear model first")) }else{
              
            calData$T2 <<- calData$Temperature
            
            lminverserec <<-  do.call(rbind,lapply(1:nrow(recData_byS), function(x){
                a <- predictTclassic(calData, targety=recData_byS$D47[x], model='wlm', replicates=replicates, bootDataset=!classicPredictions)
                b <- predictTclassic(calData, targety=recData_byS$D47[x]+recData_byS$D47error[x], model='wlm', replicates=replicates, bootDataset=!classicPredictions)
                cbind.data.frame("D47"=recData_byS$D47[x],'D47se'=recData_byS$D47error[x], "Tc"=a$temp, "se"=a$temp-b$temp)
                } ))
              
              colnames(lminverserec)[3] <- 'Tc'
              lminverserec$D47se <- recData_byS$D47error
              lminverserec[,c("D47",'D47se', "Tc", "se")]
            
              
            lminverserecwun <- lminverserec
            df3<-lminverserecwun
            
            names(df3) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "SE (1SD) Temperature (°C)")
            rownames(df3) <- NULL
            
            
            
            output$lminverserecswun <- renderTable({
              
              df3$`Δ47 (‰)` <- formatC(df3$`Δ47 (‰)`, digits = 3, format = "f")
              df3$`Δ47 (‰) error` <- formatC(df3$`Δ47 (‰) error`, digits = 4, format = "f")
              df3$`Temperature (°C)` <- formatC(df3$`Temperature (°C)`, digits = 1, format = "f")
              df3$`SE (1SD) Temperature (°C)` <- formatC(df3$`SE (1SD) Temperature (°C)`, digits = 7, format = "f")
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
          if( input$simulateYork_measuredRec == TRUE ) {
            calData$T2 <<- calData$Temperature
            if(is.null(yorkcals) & classicPredictions) { print(noquote("Please run the calibration step for the York linear model first")) }else{
              
            yorkrec <<-   do.call(rbind,lapply(1:nrow(recData_byS), function(x){
                a <- predictTclassic(calData, targety=recData_byS$D47[x], model='York', replicates=replicates, bootDataset=!classicPredictions)
                b <- predictTclassic(calData, targety=recData_byS$D47[x]+recData_byS$D47error[x], model='York', replicates=replicates, bootDataset=!classicPredictions)
                cbind.data.frame("D47"=recData_byS$D47[x],'D47se'=recData_byS$D47error[x], "Tc"=a$temp, "se"=a$temp-b$temp)
              } ))
              
              colnames(yorkrec)[3] <- 'Tc'
              yorkrec$D47se <- recData_byS$D47error
              yorkrec[,c("D47",'D47se', "Tc", "se")]
            
            
            df5 <- yorkrec
            names(df5) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "SE (1SD) Temperature (°C)")
            rownames(df5) <- NULL
            
            
            output$yorkrecswun <- renderTable({
              
              df5$`Δ47 (‰)` <- formatC(df5$`Δ47 (‰)`, digits = 3, format = "f")
              df5$`Δ47 (‰) error` <- formatC(df5$`Δ47 (‰) error`, digits = 4, format = "f")
              df5$`Temperature (°C)` <- formatC(df5$`Temperature (°C)`, digits = 1, format = "f")
              df5$`SE (1SD) Temperature (°C)` <- formatC(df5$`SE (1SD) Temperature (°C)`, digits = 7, format = "f")
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
          if( input$simulateDemingRec == TRUE ) {
            calData$T2 <<- calData$Temperature
            if(is.null(demingcals) & classicPredictions) { print(noquote("Please run the calibration step for the York linear model first")) }else{
              
            demingrec <<-
              
              demingrec<-  do.call(rbind,lapply(1:nrow(recData_byS), function(x){
                a <- predictTclassic(calData, targety=recData_byS$D47[x], model='Deming', replicates=replicates, bootDataset=!classicPredictions)
                b <- predictTclassic(calData, targety=recData_byS$D47[x]+recData_byS$D47error[x], model='Deming', replicates=replicates, bootDataset=!classicPredictions)
                cbind.data.frame("D47"=recData_byS$D47[x],'D47se'=recData_byS$D47error[x], "Tc"=a$temp, "se"=a$temp-b$temp)
                
                } ))
            
            
              colnames(demingrec)[3] <- 'Tc'
              demingrec$D47se <- recData_byS$D47error
              demingrec[,c("D47",'D47se', "Tc", "se")]
            

            df7 <- demingrec
            names(df7) <- c("Δ47 (‰)", "Δ47 (‰) error", "Temperature (°C)", "SE (1SD) Temperature (°C)")
            rownames(df7) <- NULL

            output$demingrecswun <- renderTable({
              
              df7$`Δ47 (‰)` <- formatC(df7$`Δ47 (‰)`, digits = 3, format = "f")
              df7$`Δ47 (‰) error` <- formatC(df7$`Δ47 (‰) error`, digits = 4, format = "f")
              df7$`Temperature (°C)` <- formatC(df7$`Temperature (°C)`, digits = 1, format = "f")
              df7$`SE (1SD) Temperature (°C)` <- formatC(df7$`SE (1SD) Temperature (°C)`, digits = 7, format = "f")
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
  