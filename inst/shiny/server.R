library(shiny)
package.embryogrowth <- require('HelpersMG')

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  # Expression that generates a histogram. The expression is
  # wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should re-execute automatically
  #     when inputs change
  #  2) Its output type is a plot
  
  # output$resultsInfo <- renderPrint({print(input$MonthRef)})
  
 
  output$AW <- renderText({
    buttonP <- input$goButton
    
    if (buttonP !=0) {
      dataIn <- input$data
      if (dataIn == "") dataIn <- "12.4 13.3 17.2;12.9 19.20 32.9"
      
      dataIn <- gsub("  ", " ", dataIn)
      dataIn <- gsub("; ", ";", dataIn)
      dataIn <- gsub(" ;", ";", dataIn)
      
      l <- strsplit(dataIn, ";")[[1]]
      table.list <- list()
      for (i in 1:length(l)) {
        k <- as.numeric(strsplit(l[i], " ")[[1]])
        table.list <- c(table.list, list(k))
      }
      table <- t(as.data.frame(table.list))
      colnames(table) <- NULL
      rownames(table) <- NULL
      output$DataOut <- renderTable(table)
      
    if (input$type == "Data") {
      outaw <- data.comparison(table, 
                               var.equal=ifelse(input$var.equal=="TRUE", TRUE, FALSE), 
                               criterion = input$Criterion)
    } else {
      outaw <- table.comparison(table)
    }
      outtext <- ""
      for (i in 1:length(outaw)) outtext <- paste0(outtext, names(outaw[i]), "=", outaw[i], "   ")
      outtext
    }
})


})