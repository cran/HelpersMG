library(shiny)

# require('HelpersMG')

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  # Expression that generates a histogram. The expression is
  # wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should re-execute automatically
  #     when inputs change
  #  2) Its output type is a plot
  
  # output$resultsInfo <- renderPrint({print(input$MonthRef)})



  output$typeC <- renderUI({
    if (input$type == 3)
    textInput(inputId="probs", label="Probability for each category (values separated by space)", 
              value="",
              width = "100%", 
              placeholder=c("0.1 0.2 0.7"))
    
  })
  
  output$typeD <- renderUI({
    if (input$type == 1) {
  textInput(inputId="data", label="Data to be analyzed (values separated by space and series by \";\")", 
            value="",
            width = "100%", 
            placeholder=c("12.4 13.3 17.2;12.9 19.20 32.9"))
    } else {
      textInput(inputId="data", label="Data to be analyzed (values separated by space and series by \";\")", 
                value="",
                width = "100%", 
                placeholder=c("10 12 63;20 24 70"))
    }
  })
  
  output$varData <- renderUI({
    if (input$type == 1)
    radioButtons("var.equal", "Do the variances of the series should be supposed equal?",
                 list("Yes"=1, "No"=2), selected=1, inline = TRUE)

  })
  

  
  outtextfn <- eventReactive(input$goButton, {
      dataIn <- input$data
      if (input$type == 3) {
        probs <- input$probs
        if (probs == "") probs <- "0.1 0.2 0.7"
        probs <- gsub("  ", " ", probs)
        probs <- as.numeric(strsplit(probs, " ")[[1]])
        probs <- probs/sum(probs)
      }
      
      
      if (dataIn == "") if (input$type == 1) {
        dataIn <- "12.4 13.3 17.2;12.9 19.20 32.9"
      } else {
        dataIn <- "10 12 63;20 24 70"
      }
      
      
      dataIn <- gsub("  ", " ", dataIn)
      dataIn <- gsub("; ", ";", dataIn)
      dataIn <- gsub(" ;", ";", dataIn)
      
      l <- strsplit(dataIn, ";")[[1]]
      table.list <- list()
      elements <- NULL
      for (i in 1:length(l)) {
        k <- as.numeric(strsplit(l[i], " ")[[1]])
        table.list <- c(table.list, list(k))
        elements <- c(elements, length(k))
      }
      
      if (all(elements == elements[1])) {
        table <- t(as.data.frame(table.list))
      } else {
        lm <- max(elements)
        for (i in seq_along(elements)) {
          if (elements[i] != lm) {
            table.list[[i]] <- c(table.list[[i]], rep(NA, lm-elements[i]))
          }
        }
        table <- t(as.data.frame(table.list))
      }
      colnames(table) <- NULL
      rownames(table) <- paste("Series", 1:nrow(table))
      output$DataOut <- renderTable(table, rownames=TRUE, colnames = FALSE)
      
    if (input$type == 1) {
      outaw <- series.compare(table, 
                               var.equal=ifelse(input$var.equal==1, TRUE, FALSE), 
                               criterion = c("BIC"))
    } else {
      if (input$type == 2) {
        if (all(elements == elements[1])) {
        outaw <- contingencyTable.compare(table, criterion = c("BIC"))
        } else {
          outaw <- NULL
          outtext <- "The contingency table must have all series with the same number of columns."
        }
      } else {
        # je suis dans un test de conformité à probs
        if (all(elements == length(probs))) {
          outaw <- contingencyTable.compare(table, criterion = c("BIC"), 
                                            probs=probs)
        } else {
          outaw <- NULL
          outtext <- "The contingency table must have all series with the same number of columns than the probabilities to be tested."
        }
      }
    }
      if (!is.null(outaw)) {
        outtext <- "w-value="
        for (i in 1:length(outaw)) outtext <- paste0(outtext, names(outaw[i]), "=", outaw[i], "\n")
      }
      outtext

})

output$AW <- renderText({
  outtextfn()
})
  
})