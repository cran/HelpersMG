# library(shiny); library("embryogrowh"); runApp("/Users/marc/Documents/Espace_de_travail_R/shiny/tsd")
# library(shiny); runApp("http://max3.ese.u-psud.fr:3838/phenology/")


library(shiny)
library(HelpersMG)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  
  wellPanel(
    h1("Simple replacement of chi2 and t-tests", align = "center")
    , h2("that do not use p-value", align = "center")
    , p("This web server version v. 1.00 is a simplified version of the complete tools available ", 
        a("here."
          , href="https://cran.r-project.org/package=HelpersMG"
          , target="_blank"))
    , p("HelpersMG package is developped by "
        , a("Marc Girondot"
            , href="http://max2.ese.u-psud.fr/epc/conservation/Girondot/Publications/Marc.html"
            , target="_blank"))
  ),
    
    # Show a plot of the generated distribution
    mainPanel(
      
      selectInput("type", "type"
                  , choices=as.list(c("Contingency table", "Data"))
                  , selected = "Data", multiple = FALSE,
                  selectize = TRUE, width = NULL, size = NULL),
      selectInput("var.equal", "var.equal"
                  , choices=as.list(c(FALSE, TRUE))
                  , selected = TRUE, multiple = FALSE,
                  selectize = TRUE, width = NULL, size = NULL),
      selectInput("Criterion", "Criterion"
                  , choices=as.list(c("AIC", "AICc", "BIC"))
                  , selected = "AIC", multiple = TRUE,
                  selectize = TRUE, width = NULL, size = NULL),
      textInput("data", "Data to be analyzed (values separated by space and lines by ;)", "", width = "100%", 
                placeholder=c("12.4 13.3 17.2;12.9 19.20 32.9")),
      actionButton("goButton", "Compare"), 
      verbatimTextOutput(outputId="AW"), 
      p("The value represents the probability that a single model for all the series is sufficient 
        to model the data. The lower the value (close to 0), the most likely the series are different. 
        The higher the value (close to 1) the most likely the series are similar."),
      tableOutput(outputId="DataOut")
    )
  )
  )