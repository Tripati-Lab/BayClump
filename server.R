# Define server logic
server <- function(input, output, session) { 

# Calibration tab
  
  output$BayClump_template.csv <- downloadHandler(
       filename = "BayClump_template.csv",
       content = function(file) {
         write.csv(BayClump_template, file, row.names = FALSE)
       }
     )
  
  d1 <- rnorm(300, mean = 4, sd = 2)
  d2 <- rnorm(300, mean = 5, sd = 1.8)
  
  output$dummyplot <- renderPlot(
    plot(d1, d2, main = "Summary stats here")
  )
  
  output$dummytext <- renderUI({
    HTML("What output options do we want?")})
  
  #Fit the regression models under default parameters
  RF_improper<-fitClumpedRegressions(calibrationData=Petersen,
                                     hasMaterial = F,n.iter= 1000, burninFrac=.5)
  RP_improper<-regressionParameters(CompleteModelFit=RF_improper)
  
  output$calibration_table <- DT::renderDataTable({
    datatable(RP_improper) 
  })
  

# Calibration plots tab

# Create dummy data
  trace_0 <- rnorm(100, mean = 5)
  trace_1 <- rnorm(100, mean = 0)
  trace_2 <- rnorm(100, mean = -5)
  x <- c(1:100)

  data <- data.frame(x, trace_0, trace_1, trace_2)

  output$examplefig <- renderPlotly({

  fig <- plot_ly(data, x = ~x)
  fig <- fig %>% add_trace(y = ~trace_0, name = 'trace 0',mode = 'lines')
  fig <- fig %>% add_trace(y = ~trace_1, name = 'trace 1', mode = 'lines+markers')
  fig <- fig %>% add_trace(y = ~trace_2, name = 'trace 2', mode = 'markers')
  fig <- fig %>% layout(legend=list(title=list(text='<b> This is example code </b>')))

    return(fig)
})

# Reconstruction tab
  
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
    as.data.frame(data2)%>% 
      mutate_if(is.numeric, round, digits = 3)
},

options = list(paging=FALSE, scrollX = TRUE)
)
  
  
}
  