# library(shiny); library("embryogrowth"); runApp("/Users/marcgirondot/Documents/Espace_de_travail_R/_shiny/compare")
# library(shiny); runApp("http://134.158.74.46/compare/")


library(shiny)
package.HelpersMG <- require('HelpersMG')


# Define UI for application that draws a histogram
fluidPage(
  
  # Application title
  
  wellPanel(
    h1(HTML("The <em>w</em>-value: An alternative to <em>t</em>- and <em>&chi;<sup>2</sup></em>-tests"), align = "center")
    , h2(HTML("that does not use <em>p</em>-value"), align = "center")
    , p("This web server version v. 1.04 is a simplified version of the complete tools available ", 
        a("in the HelpersMG R package."
          , href="https://cran.r-project.org/package=HelpersMG"
          , target="_blank"))
    , p(HTML("The <em>w</em>-value methodology has been developped by 
             <a href=\"http://max2.ese.u-psud.fr/epc/conservation/Girondot/Publications/Marc.html\">Marc Girondot</a> and 
             <a href=\"http://www.ese.u-psud.fr/article199.html\">Jean-Michel Guillon</a>."))
    , p("Ecologie, Systématique, Evolution - Université Paris-Sud, CNRS, AgroParisTech, Université Paris Saclay, France.")
    , p("It has been published in: ")
    , p(HTML("<a href=\"https://www.researchgate.net/publication/323847382_The_w-value_An_Alternative_to_t-and_ch2_Tests\">Girondot, M., Guillon, J.-M., 2018. The <em>w</em>-value: An alternative to <em>t</em>- and &chi;<sup>2</sup> tests. Journal of Biostatistics & Biometrics 1(1): 1-4.</a>."))
    , p("")
  ),
    
    # Show a plot of the generated distribution
    mainPanel(
      
      radioButtons("type", "What test should be done ?", 
                   list("Data series"=1, "Homogeneity contingency table"=2, "Conformity contingency table"=3), selected=1, inline = TRUE),
      # selectInput("type", "type"
      #            , choices=as.list(c("Contingency table", "Data"))
      #            , selected = "Data", multiple = FALSE,
      #            selectize = TRUE, width = NULL, size = NULL),
      # selectInput("var.equal", "var.equal"
      #             , choices=as.list(c(FALSE, TRUE))
      #             , selected = TRUE, multiple = FALSE,
      #             selectize = TRUE, width = NULL, size = NULL),
      uiOutput("varData"), 
      # radioButtons("var.equal", "Do the variances of the series should be supposed equal? (for data series test)", 
      #              list("Yes"=1, "No"=2), selected=1, inline = TRUE),
      # selectInput("Criterion", "Criterion"
      #             , choices=as.list(c("AIC", "AICc", "BIC"))
      #             , selected = "AIC", multiple = TRUE,
      #             selectize = TRUE, width = NULL, size = NULL),
      # textInput(inputId="data", label="Data to be analyzed (values separated by space and series by \";\")", 
      #           value="",
      #           width = "100%", 
      #           placeholder=c("12.4 13.3 17.2;12.9 19.20 32.9")),
      uiOutput("typeD"), 
      uiOutput("typeC"), 
      actionButton("goButton", "Compare"), 
      # submitButton("Compare"),
      h4("Results"),
      verbatimTextOutput(outputId="AW"), 
      p("The value \"w-value=BICw(identical)=\" represents the probability that, taking into the available information, a single model for all the series is sufficient 
        to modelize the data. The lower the value (closer to 0), the most likely the series are obtained from different distributions. 
        The higher the value (closer to 1) the most likely the series are obtained from a single distribution."),
      h4("Data table"), 
      tableOutput(outputId="DataOut")
    )
  )
  