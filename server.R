# Define server logic
server <- function(input, output, session) { 
  options(shiny.maxRequestSize=800*1024^2) 
  
# Calibration tab
  
  output$BayClump_cal_temp <- downloadHandler(
       filename = "BayClump_calibration_template.csv",
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
       hasMaterial <<- ifelse( is.na(calibrationData()$Material), FALSE, TRUE )
       
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
       if("Bayesian model no errors" %in% names(wb) == TRUE) 
         {removeWorksheet(wb, "Bayesian model no errors") & removeWorksheet(wb, "Bayesian model no errors CI")}
       if("Bayesian model with errors" %in% names(wb) == TRUE) 
         {removeWorksheet(wb, "Bayesian model with errors") & removeWorksheet(wb, "Bayesian model with errors CI")}

       calData <<- NULL
       calData <<- calibrationData()
       
       ##If temperature is in degree celsius
       
       calData$T2 <<- (10^6)/(calData$Temperature + 273.15)^2 
       #calData$TempError <- (10^6)/(calData$TempError + 273.15)^2
       
       if(input$uncertainties == "usedaeron") { # Placeholder for Daeron et al. uncertainties
         calData$TempError <<- 1
       }
       
#       if(input$misc == "scale") {
#         calData$Temperature <- scale(calData$Temperature)
#         calData$T2 <- scale(calData$T2)
#         calData$TempError <- scale(calData$TempError)
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
           lmcals <<- simulateLM_measured(calData, replicates = replicates)
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
            do.call(rbind.data.frame,apply(lmcals, 2, function(x){
               cbind.data.frame(Median= round(median(x), 4), `Lower 95% CI`=round(quantile(x, 0.025), 4),  `Upper 95% CI`=round(quantile(x, 0.975), 4))
             }))
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
             do.call(rbind.data.frame,apply(lminversecals, 2, function(x){
               cbind.data.frame(Median= round(median(x), 4), `Lower 95% CI`=round(quantile(x, 0.025), 4),  `Upper 95% CI`=round(quantile(x, 0.975), 4))
             }))
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
            do.call(rbind.data.frame,apply(yorkcals, 2, function(x){
              cbind.data.frame(Median= round(median(x), 4), `Lower 95% CI`=round(quantile(x, 0.025), 4),  `Upper 95% CI`=round(quantile(x, 0.975), 4))
            }))
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
             do.call(rbind.data.frame,apply(demingcals, 2, function(x){
               cbind.data.frame(Median= round(median(x), 4), `Lower 95% CI`=round(quantile(x, 0.025), 4),  `Upper 95% CI`=round(quantile(x, 0.975), 4))
             }))
           })
           
           }

    #     checkboxInput("linear", "Linear model", FALSE),
         if(input$simulateBLM_measuredMaterial != FALSE) {
           sink(file = "Bayeslinmodtext.txt", type = "output")
           bayeslincals <- simulateBLM_measuredMaterial(calData, generations=1000, replicates = replicates, isMixed=F)
           sink()
           
           bayeslincinoerror <- RegressionSingleCI(data = bayeslincals$BLM_Measured_no_errors, from = min(calData$T2), to = max(calData$T2))
           bayeslincalcinoerror <- as.data.frame(bayeslincinoerror)
           bayeslinciwitherror <- RegressionSingleCI(data = bayeslincals$BLM_Measured_errors, from = min(calData$T2), to = max(calData$T2))
           bayeslincalciwitherror <- as.data.frame(bayeslinciwitherror)
           
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
                       error_x = ~list(array = ~TempError, color = '#000000'),
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
  
  read_batch_with_progress = function(file_path,nrows,no_batches){
    progress = Progress$new(session, min = 1,max = no_batches)
    progress$set(message = "Reading reconstruction data")
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
  
  
  reconstructionData = reactive({
             req(input$reconstructiondata)
             n_rows = length(count.fields(input$reconstructiondata$datapath))
             df_out = read_batch_with_progress(input$reconstructiondata$datapath,n_rows,10)
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
    hasMaterial <- ifelse( is.na(reconstructionData()$Material), FALSE, TRUE )
    
    # Update the number of bootstrap replicates to run based on user selection
  #  replicates <- ifelse(input$replication == "50", 50,
  #                       ifelse(input$replication == "100", 100,
   #                             ifelse(input$replication == "500", 500,
    #                                   ifelse(input$replication == "1000", 1000, NA))))
    
    # Remove existing worksheets from wb on 'run' click, if any
#    if("Linear regression" %in% names(wb) == TRUE) 
 #   {removeWorksheet(wb, "Linear regression") & removeWorksheet(wb, "Linear regression CI")}
  #  if("Inverse linear regression" %in% names(wb) == TRUE) 
#    {removeWorksheet(wb, "Inverse linear regression") & removeWorksheet(wb, "Inverse linear regression CI")}
 #   if("York regression" %in% names(wb) == TRUE) 
  #  {removeWorksheet(wb, "York regression") & removeWorksheet(wb, "York regression CI")}
#    if("Deming regression" %in% names(wb) == TRUE) 
#    {removeWorksheet(wb, "Deming regression") & removeWorksheet(wb, "Deming regression CI")}
#    if("Bayesian model no errors" %in% names(wb) == TRUE) 
#    {removeWorksheet(wb, "Bayesian model no errors") & removeWorksheet(wb, "Bayesian model no errors CI")}
#    if("Bayesian model with errors" %in% names(wb) == TRUE) 
#    {removeWorksheet(wb, "Bayesian model with errors") & removeWorksheet(wb, "Bayesian model with errors CI")}
    
    recData <- NULL
    recData <- reconstructionData()

    withProgress(message = 'Running selected reconstructions, please wait', {
      # if(is.null(output$contents)) {
      #   print("Please upload calibration data")
      # }
      
      # Run prediction function
      if( !is.null(lmcals) ) {
        
        #output$linready <- renderPrint(noquote("Linear calibration model complete"))
        

        lmrec <- predictTcNonBayes(data=cbind(recData$D47,recData$D47error), 
                          slope=median(lmcals$slope), 
                          slpcnf=CItoSE(quantile(lmcals$slope, 0.975), quantile(lmcals$slope, 0.025)), 
                          intercept=median(lmcals$intercept), 
                          intcnf=CItoSE(quantile(lmcals$intercept, 0.975), quantile(lmcals$intercept, 0.025)))
        
        print(noquote("Linear reconstruction complete"))
        
        output$lmrecs <- renderPrint(
          return(as.data.frame(lmrec))
        )
      }

})
    
  })
  
  
  output$recresults <- renderPrint({
    recresult()
  })

}
  