# Define server logic
server <- function(input, output, session) { 
  options(shiny.maxRequestSize=800*1024^2) 
  
# Calibration tab
  
  output$BayClump_cal_temp <- downloadHandler(
       filename = "BayClump_template.csv",
       content = function(file) {
         write.csv(BayClump_calibration_template, file, row.names = FALSE)
       }
     )
  
  read_batch_with_progress = function(file_path,nrows,no_batches){
    progress = Progress$new(session, min = 1,max = no_batches)
    progress$set(message = "Reading calibration data")
    seq_length = ceiling(seq.int(from = 2, to = nrows-2,length.out = no_batches+1))
    seq_length = seq_length[-length(seq_length)]
    
    #read the first line
    df = read.csv(file_path,skip = 0,nrows = 1)
    col_names = colnames(df)
    
    for(i in seq_along(seq_length)){
      progress$set(value = i)
      if(i == no_batches) chunk_size = -1 else chunk_size = seq_length[i+1] - seq_length[i]
      
      df_temp = read.csv(file_path, skip = seq_length[i], nrows = chunk_size,header = FALSE,stringsAsFactors = FALSE)
      colnames(df_temp) = col_names
      df = rbind(df,df_temp)
    }
    
    progress$close()
    return(df)
  }
  
  
  calibrationData = reactive({
    switch(input$calset,
           'model1' = return(Petersen),
           'model2' = return(Anderson),
           'model1and2' = return(PetersenAnderson),
           'mycal' = reactiveValues({
    req(input$calibrationdata)
    n_rows = length(count.fields(input$calibrationdata$datapath))
    df_out = read_batch_with_progress(input$calibrationdata$datapath,n_rows,10)
    return(df_out)
  }),
            'all' = reactiveValues({
    req(input$calibrationdata)
    n_rows = length(count.fields(input$calibrationdata$datapath))
    df_out = read_batch_with_progress(input$calibrationdata$datapath,n_rows,10)
    alldat <- rbind(df_out, PetersenAnderson)
    return(alldat)
            })
    )
  })
  
  if(exists("wb")) rm(wb) # Delete an existing workbook in preparation for new results
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
       
       # Update the number of bootstrap replicates to run based on user selection
       replicates <- ifelse(input$replication == "50", 50,
                            ifelse(input$replication == "100", 100,
                                   ifelse(input$replication == "500", 500,
                                          ifelse(input$replication == "1000", 1000, NA))))

       # Remove existing worksheets from wb on 'run' click, if any
       if("Linear regression" %in% names(wb) == TRUE) 
         {removeWorksheet(wb, "Linear regression") & removeWorksheet(wb, "Linear regression CI")}
       if("Inverse linear regression" %in% names(wb) == TRUE) 
         {removeWorksheet(wb, "Inverse linear regression") & removeWorksheet(wb, "Inverse linear regression CI")}
       if("York regression" %in% names(wb) == TRUE) 
         {removeWorksheet(wb, "York regression") & removeWorksheet(wb, "York regression CI")}
       if("Deming regression" %in% names(wb) == TRUE) 
         {removeWorksheet(wb, "Deming regression") & removeWorksheet(wb, "Deming regression CI")}
       if("Bayesian linear model" %in% names(wb) == TRUE) 
         {removeWorksheet(wb, "Bayesian linear model") & removeWorksheet(wb, "Bayesian linear model CI")}       

       calData <- NULL
       calData <- calibrationData()
       
       ##If temperature is in degree celsius
       
       calData$T2 <- (10^6)/(calData$Temperature + 273.15)^2 
       calData$Temp_Error <- (10^6)/(calData$Temp_Error + 273.15)^2
       
       if(input$uncertainties == "usedaeron") { # Placeholder for Daeron et al. uncertainties
         calData$Temp_Error <- 1
       }
       
#       if(input$misc == "scale") {
#         calData$Temperature <- scale(calData$Temperature)
#         calData$T2 <- scale(calData$T2)
#         calData$Temp_Error <- scale(calData$Temp_error)
#         calData$D47 <- scale(calData$D47)
#         calData$D47error <- scale(calData$D47error)
#       }

       withProgress(message = 'Running selected models, please wait', {
        # if(is.null(output$contents)) {
        #   print("Please upload calibration data")
        # }
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
         
         
         
         if(input$simulateLM_measured != FALSE) {
           sink(file = "linmodtext.txt", type = "output")
           lmcals <- simulateLM_measured(calData, replicates = replicates)
           sink()
           
           lmci <- RegressionSingleCI(data = lmcals, from = min(calData$T2), to = max(calData$T2))
           lmcalci <- as.data.frame(lmci)
           
           output$lmcalibration <- renderPlotly({
             lmfig <- plot_ly(calibrationData()
             )
             lmfig <- lmfig %>%
                 add_trace(x = ~(10^6)/(calibrationData()$Temperature + 273.15)^2, 
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
             head(lmcalci2)
           })
           
         }
         
         if(input$simulateLM_inverseweights != FALSE) {
           sink(file = "inverselinmodtext.txt", type = "output")
           lminversecals <- simulateLM_inverseweights(calData, replicates = replicates)
           sink()
           
           lminverseci <- RegressionSingleCI(data = lminversecals, from = min(calData$T2), to = max(calData$T2))
           lminversecalci <- as.data.frame(lminverseci)
           
           output$lminversecalibration <- renderPlotly({
             lminversefig <- plot_ly(data = calibrationData()
             )
             lminversefig <- lminversefig %>% 
               add_trace(x = ~(10^6)/(calibrationData()$Temperature + 273.15)^2, 
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
             head(lminversecalci2)
           })
           
         }
         
         if(input$simulateYork_measured != FALSE) {
           sink(file = "yorkmodtext.txt", type = "output")
          yorkcals <- simulateYork_measured(calData, replicates = replicates)
          sink()
          
          yorkci <- RegressionSingleCI(data = yorkcals, from = min(calData$T2), to = max(calData$T2))
          yorkcalci <- as.data.frame(yorkci)
          
          output$yorkcalibration <- renderPlotly({
            yorkfig <- plot_ly(data = calibrationData()
            )
            yorkfig <- yorkfig %>% 
              add_trace(x = ~(10^6)/(calibrationData()$Temperature + 273.15)^2, 
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
            head(yorkcalci2)
          })
        
         }

         if(input$simulateDeming != FALSE) {
           sink(file = "demingmodtext.txt", type = "output")
           demingcals <- simulateDeming(calData, replicates = replicates)
            sink()
            
            demingci <- RegressionSingleCI(data = demingcals, from = min(calData$T2), to = max(calData$T2))
            demingcalci <- as.data.frame(demingci)
            
            output$demingcalibration <- renderPlotly({
              demingfig <- plot_ly(data = calibrationData()
              )
              demingfig <- demingfig %>% 
                add_trace(x = ~(10^6)/(calibrationData()$Temperature + 273.15)^2, 
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
             head(demingcalci2)
           })
           
           }

    #     checkboxInput("linear", "Linear model", FALSE),
         if(input$simulateBLM_measuredMaterial != FALSE) {
           sink(file = "Bayeslinmodtext.txt", type = "output")
           bayeslincals <- simulateBLM_measuredMaterial(calData, generations=1000, replicates = replicates, isMixed=F)
           sink()
           
           bayeslinci <- RegressionSingleCI(data = bayeslincals, from = min(calData$T2), to = max(calData$T2))
           bayeslincalci <- as.data.frame(bayeslinci)
           
           output$bayeslincalibration <- renderPlotly({
             bayeslinfig <- plot_ly(data = calibrationData()
             )
             bayeslinfig <- bayeslinfig %>% 
               add_trace(x = ~(10^6)/(calibrationData()$Temperature + 273.15)^2, 
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
               add_ribbons(data = bayeslincalci,
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
               add_lines(data = bayeslincalci,
                         x = ~x,
                         y = ~median_est,
                         name = 'Median estimate',
                         line = list(color = "black", dash = 'dash'),
                         hovertemplate = paste(
                           "Temperature (10<sup>6</sup>/T<sup>2</sup>): %{x}<br>",
                           "Δ<sub>47</sub> (‰): %{y}<br>"))
             bayeslinfig <- bayeslinfig %>% layout(title = '<b> Bayesian linear calibration model </b>',
                                                     legend=list(title=list(text='Legend')),
                                                     xaxis = list(title = 'Temperature (10<sup>6</sup>/T<sup>2</sup>)', hoverformat = '.1f'), 
                                                     yaxis = list(title = 'Δ<sub>47</sub> (‰)', hoverformat = '.3f'))
             
             return(bayeslinfig)
           })

           addWorksheet(wb, "Bayesian linear model") # Add a blank sheet
           addWorksheet(wb, "Bayesian linear model CI") # Add a blank sheet 
           
           bayeslincalci2 <- bayeslincalci
           names(bayeslincalci2) <- c("10^6/T^2", "D47_median_est",	"D47_ci_lower_est",	"D47_ci_upper_est")
           
           writeData(wb, sheet = "Bayesian linear model", cbind((bayeslincals)[i], bayeslincals[[i]])) # Write regression data
           writeData(wb, sheet = "Bayesian linear model CI", bayeslincalci2)
           
           
           print(noquote("Bayesian linear model complete"))
           output$blin <- renderPrint({
             
             for (i in 1:length(bayeslincals)) {
               print(names(bayeslincalci2)[i])
               print(head(bayeslincalci2[[i]]))
             }
             
           })
           
         }
         
 
    #     calibrationData2 <<- calibrationData()
         
         #sapply(list.files('Functions', full.names = T), source)
           #        sapply(input$selectmods, do.call, args = list())
                   #source(fitClumpedRegressions(x))  
       })
       
       
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
    output$rawcaldata <- renderPlotly({
  rawcalfig <- plot_ly(calibrationData(), 
                       x = ~Temperature, 
                       y = ~D47, 
                       type = 'scatter', 
                       mode = 'lines+markers', 
                       linetype = ~Material, 
                       color = ~Mineralogy,
                       colors = viridis_pal(option = "D", end = 0.9)(minlength),
                       opacity = 0.6,
                       error_y = ~list(array = ~D47error, color = '#000000'),
                       error_x = ~list(array = ~Temp_Error, color = '#000000'),
                       text = as.character(calibrationData()$Sample.Name),
                       hovertemplate = paste(
                         "<b>Sample: %{text}</b><br><br>",
                         "Temperature (°C): %{x}<br>",
                         "Δ<sub>47</sub> (‰): %{y}<br>",
                         "Mineralogy: ", as.character(calibrationData()$Mineralogy),"<br>",
                         "Type: ", as.character(calibrationData()$Material),
                         "<extra></extra>"))
  rawcalfig <- rawcalfig %>% layout(title = '<b> Raw calibration data from user input </b>',
                             legend=list(title=list(text='Material and mineralogy')),
                             xaxis = list(title = 'Temperature (°C)'), 
                             yaxis = list(title = 'Δ<sub>47</sub> (‰)', hoverformat = '.3f'))

    return(rawcalfig)
})
  })

# Reconstruction tab

  output$BayClump_reconstruction_template.csv <- downloadHandler(
    filename = "BayClump_reconstruction_template.csv",
    content = function(file) {
      write.csv(BayClump_reconstruction_template, file, row.names = FALSE)
    }
  )
  
  # Create dummy data
  trace_3 <- rnorm(100, mean = 3)
  trace_4 <- rnorm(100, mean = 0.3)
  trace_5 <- rnorm(100, mean = -3)
  x2 <- c(1:100)
  
  data2 <- data.frame(x2, trace_3, trace_4, trace_5)
  
  output$examplefig2 <- renderPlotly({
    
    fig2 <- plot_ly(data2, x = ~x2)
    fig2 <- fig2 %>% add_trace(y = ~trace_3, name = 'trace 3',mode = 'lines')
    fig2 <- fig2 %>% add_trace(y = ~trace_4, name = 'trace 4', mode = 'lines+markers')
    fig2 <- fig2 %>% add_trace(y = ~trace_5, name = 'trace 5', mode = 'markers')
    fig2 <- fig2 %>% layout(legend=list(title=list(text='<b> This is more example code </b>')))
    
    return(fig2)
  })
  output$table1 <- renderDataTable({
    as.data.frame(data2) %>% 
      mutate_if(is.numeric, round, digits = 3)
},

  rownames=FALSE, options = list(paging=FALSE, scrollX = TRUE)
)
  
  
# Petersen et al tab
  
  output$Petersendat <- renderDataTable({
    Petersen2 <- as.data.frame(Petersen) %>% 
      mutate_if(is.numeric, round, digits = 3)
    return(Petersen2)
  },
  
  rownames=FALSE, options = list(paging=FALSE, scrollX = TRUE))
  
}

  