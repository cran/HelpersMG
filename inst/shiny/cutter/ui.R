# library(shiny); runApp("/Users/marcgirondot/Documents/Espace_de_travail_R/_shiny/BoneProfileR")
# library(shiny); runApp("http://134.158.74.46/BoneProfileR/")


library(shiny)
if (!requireNamespace("shinyWidgets", quietly = TRUE)) {
  progressbar <- FALSE
} else {
  progressbar <- TRUE
  library("shinyWidgets")
}

package.HelpersMG <- require('HelpersMG')
version <- "5.1 build 1166"


# Define UI for application that draws a histogram
fluidPage(
  titlePanel(h1("Cutte",
                img(src="Rlogo.png", height=40, width=40), align = "center"), 
             windowTitle = "CutteR"), 
  p(HTML("<b><a href=\"https://max2.ese.u-psud.fr/epc/conservation/index.html\">Marc Girondot</a></b> - Laboratoire Ecologie, Systématique, Evolution"), align = "center"),
  p(HTML("Université Paris-Saclay, CNRS, AgroParisTech, France."), align = "center"), 
  
  wellPanel(
    
    p(HTML("<strong>CutteR is a scientific method and a software used 
  to model left and right censored and truncated data.</strong>"), align = "left"),
    p("It used a Bayesian MCMC distribution fit using conditional gamma, lognormal, normal, Weibull, or generalized gamma distribution. The priors are based on Gaussian distribution of parameters obtained from maximum likelihood with large standard deviation. The posterior distribution is obtained with 5000 MCMC iterations."), 
    p(), 
    p(HTML("<strong>Copy data in the box below, select the required options and 
         click 'Run the analysis' button.</strong>"), align = "left")
  ), 
  
  # Show a plot of the generated distribution
  wellPanel(
    textAreaInput(inputId="Data", label="Copy the data here, separated by space, semicolon, tabulation, or return; unknown values must be LDL (Lower the Detection Limit) or UDL (Upper the Detection Limit). Do not use comma to separate values as commas are used as synonyme for decimal point.", 
                  value="",
                  rows = 10, 
                  width="100%", 
                  placeholder=""), 
    p("To try the method, copy these values and set a Lower limit of detection to 10:"), 
    p("89.77 12.75 44.18 60.66 32.68 34.57 28.09 91.11 13.98 56.95 59.64 76.22 84.64 61.45 27.50 17.27 16.39 22.10 15.81 62.86 86.39 50.86 26.33 29.43 53.95 40.94 64.13 19.49 40.26 44.03 66.70 22.10 21.07 36.57 70.29 50.20 54.13 13.11 22.76 40.04 28.49 LDL 13.98 LDL 62.24 55.50 15.30 11.29 14.39 40.88 26.31 33.04 19.02 47.35 56.16 24.00 18.66 20.52 61.26 LDL 15.45 41.65 77.76 71.66 22.13 LDL 46.61 39.65 58.96 LDL 28.01 83.33 64.28 24.31 22.52 LDL 64.47 62.77 11.45 32.91 39.11 71.29 56.73 67.44 32.42 45.63 32.60 60.96 16.38 LDL 12.26 143.26 18.80 13.28 29.26 50.78 15.98 45.35 46.43 89.13")
  ), 
  
  splitLayout(
    textInput(inputId="LDL", label="Lower limit of detection", value="")
    , textInput(inputId="UDL", label="Upper limit of detection", value="")
    , cellWidths=c("50%", "50%")
  ), 
  p("Let a box being empty if the option is not necessary."),
  splitLayout(
    selectInput("censored", label="Choose which cut will be fitted "
                , choices=list("Censored"="censored", 
                               "Truncated"= "truncated")
                , selected = "censored", multiple = FALSE,
                selectize = FALSE, size = NULL),
    selectInput("mixture", label="How many mixture should be used "
                , choices=list("1" = "1", 
                               "2" = "2", 
                               "3" = "3")
                , selected = "1", multiple = FALSE,
                selectize = FALSE, size = NULL),
    selectInput("distribution", label="Choose which distribution will be fitted "
                , choices=list("Auto" ="auto", 
                               "Gamma"="gamma", 
                               "Lognormal"= "lognormal", 
                               "Normal"= "normal", 
                               "Weibull" = "weibull", 
                               "Generalized gamma" = "generalized.gamma")
                , selected = "auto", multiple = FALSE,
                selectize = FALSE, size = NULL),
    cellWidths=c("33%", "33%", "34%")
  ), 
  
  actionButton(inputId="goButton", label="Run the analysis", width="30%", 
               icon("paper-plane"), 
               style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
  HTML("<i>Be patient, the analysis requires between 30 seconds to 5 minutes to be completed.</i>"), 
  p(""), 
  if (progressbar) {progressBar(
    id = "pb2",
    value = 0,
    total = 10000,
    title = "",
    display_pct = TRUE
  )},
  wellPanel(
    
    textOutput(outputId="Result"), 
    tags$style(type="text/css", "#Result {white-space: pre-wrap;}"), 
    plotOutput("Plot")
  ), 
  wellPanel(
    p(paste0("This web server version v. ",version, " is a simplified version of the complete tools available in the HelpersMG R package.")), 
    HTML('<img src="http://perso0.free.fr/cgi-bin/wwwcount.cgi?df=Cutter.dat&dd=C&ft=0">'), 
    p(HTML("<small><i><font color='#006699'>The Virtual Data initiative, run by LABEX P2IO and supported by Université Paris-Saclay, is thanked for providing computing resources on its cloud infrastructure.</font></i></small>"))
  )
  
)
