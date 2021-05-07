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
    req(input$calibrationdata)
    n_rows = length(count.fields(input$calibrationdata$datapath))
    
    df_out = read_batch_with_progress(input$calibrationdata$datapath,n_rows,10)
    
    return(df_out)
  })
  
  
  
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

  
    wb <- createWorkbook("calibration output") # Prepare a workbook for calibration outputs
    
     modresult <- eventReactive(input$runmods, {

       calData <- NULL
       calData <- calibrationData()
       
       calData$T2 <- calData$T2 + 273.15 # Convert temperature in degrees C to Kelvin 
       calData$Temp_Error <- calData$Temp_Error + 273.15 # Convert error in degrees C to Kelvin 
       
       if(input$uncertainties == "usedaeron") { # Placeholder for Daeron et al. uncertainties
         calData$Temp_Error <- 1
       }
       
 #      if(input$misc == "scale") {
#         calData$T2 <- scale(calData$T2)
#         calData$Temp_Error <- scale(calData$Temp_error)
#         calData$D47 <- scale(calData$D47)
#         calData$D47_SD <- scale(calData$D47_SD)
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
           lmcals <- simulateLM_measured(calData)
           sink()
           #lmcal_df = data.frame(lmcals) #Class = prediction
           print(noquote("Linear regression complete"))
           output$lmcal <- renderPrint({
             head(lmcals)
           })
           sink(file = "Linear regression calibration.txt")
           print(lmcals)
           sink()
           
         }
         
         if(input$simulateLM_inverseweights != FALSE) {
           sink(file = "inverselinmodtext.txt", type = "output")
           lminversecals <- simulateLM_inverseweights(calData)
           sink()
           #lminversecals_df = data.frame(lminversecals) #Class = prediction
           print(noquote("Inverse linear regression complete"))
           output$lminversecal <- renderPrint({
             head(lminversecals)
           })
           sink(file = "Inverse linear regression calibration.txt")
           print(lminversecals)
           sink()
           
         }
         
         if(input$simulateYork_measured != FALSE) {
           sink(file = "yorkmodtext.txt", type = "output")
          yorkcals <- simulateYork_measured(calData)
          sink()
          #yorkcal_df = data.frame(yorkcals) #Class = prediction
          print(noquote("York regression complete"))
          output$york <- renderPrint({
           head(yorkcals)
          })
          sink(file = "York regression calibration.txt")
          print(yorkcals)
          sink()
        
         }

         if(input$simulateDeming != FALSE) {
           sink(file = "demingmodtext.txt", type = "output")
           demingcals <- simulateDeming(calData)
            sink()
           print(noquote("Deming regression complete"))
           output$deming <- renderPrint({
             head(demingcals)
           })
           
            sink(file = "Deming calibration.txt")
            print(demingcals)
            sink()
           }

    #     checkboxInput("linear", "Linear model", FALSE),
         if(input$simulateBLM_measuredMaterial != FALSE) {
           sink(file = "Bayeslinmodtext.txt", type = "output")
           bayeslincals <- simulateBLM_measuredMaterial(calData, generations=1000, replicates=5, isMixed=F)
           sink()
           
           print(noquote("Bayesian linear model complete"))
           output$blin <- renderPrint({
             
             for (i in 1:length(bayeslincals)) {
               print(names(bayeslincals)[i])
               print(head(bayeslincals[[i]]))
             }

           })
           
           sink(file = "Bayesian linear model calibration.txt")
           print(bayeslincals)
           sink()
         }
         
         if(input$simulateBLMM_measuredMaterial != FALSE) {
           sink(file = "BayeslinMixmodtext.txt", type = "output")
           bayeslincals <- simulateBLM_measuredMaterial(calData, generations=1000, replicates=5, isMixed=T)
           sink()
           
           print(noquote("Bayesian linear simple and mixed models complete"))
           output$blin <- renderPrint({
             
             for (i in 1:length(bayeslincals)) {
               print(names(bayeslincals)[i])
               print(head(bayeslincals[[i]]))
             }
             
           })
           
           sink(file = "Bayesian linear model calibration.txt")
           print(bayeslincals)
           sink()
         }
         
 
    #     calibrationData2 <<- calibrationData()
         
         #sapply(list.files('Functions', full.names = T), source)
           #        sapply(input$selectmods, do.call, args = list())
                   #source(fitClumpedRegressions(x))  
       })
       
      # saveWorkbook(wb, file = "tempcalworkbook.xlsx") # Save to a temporary workbook in the environment
       
                 })
     
     output$modresults <- renderPrint({
        modresult()
     })
     
     output$downloadcalibrations <- downloadHandler(
       filename = function() { 
         paste("Calibration_output_", Sys.Date(), ".zip", sep="")
       },
       content = function(file) {

         files <- ls()
         
         #create an empty list to receive the reads
         cal_list<-list()
         
         # Loop counter
         j <- 1
         
         for (i in files){
           
           #set the name of the file I'm looking for
           f <- paste0(paste(gsub('\\s+','_',i),"calibration.txt", sep=""))
           
           #check its existence and it does exist, get it to the list
           if (file.exists(f)) {
             # Store stack
             cal_list[[j]]<-stack(f)
             
             # Increment loop counter
             j <- j + 1
           }
         }
         
         zip(cal_list, file)
       })
     
    
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
                       y = ~Temperature, 
                       x = ~D47, 
                       type = 'scatter', 
                       mode = 'lines+markers', 
                       linetype = ~Material, 
                       color = ~Mineralogy,
                       colors = viridis_pal(option = "D")(minlength),
                       error_x = ~list(array = ~D47_SD, color = '#000000'),
                       error_y = ~list(array = ~Temp_Error, color = '#000000'),
                       text = as.character(calibrationData()$Sample.Name),
                       hovertemplate = paste(
                         "<b>Sample: %{text}</b><br><br>",
                         "Temperature (°C): %{y}<br>",
                         "Δ<sub>47</sub> (‰): %{x}<br>",
                         "Mineralogy: ", as.character(calibrationData()$Mineralogy),"<br>",
                         "Type: ", as.character(calibrationData()$Material),
                         "<extra></extra>"))
 # rawcalfig <- rawcalfig %>% add_trace(, name = 'Δ47', mode = 'lines+markers')
  rawcalfig <- rawcalfig %>% layout(title = '<b> Raw calibration data from user input </b>',
                             legend=list(title=list(text='Material and mineralogy')),
                             yaxis = list(title = 'Temperature (°C)'), 
                             xaxis = list(title = 'Δ<sub>47</sub> (‰)', hoverformat = '.4f'))

    return(rawcalfig)
})
  })

     # Calibration models plots
     
  observe({
    output$calibrationmodels <- renderPlot({
      plot(yorkcals)
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

  