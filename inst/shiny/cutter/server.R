library(shiny)


# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  if (!requireNamespace("shinyWidgets", quietly = TRUE)) {
    progressbar <- FALSE
  } else {
    progressbar <- TRUE
    library("shinyWidgets")
  }
  
  outplotfn <- eventReactive(eventExpr=input$goButton, ignoreNULL = FALSE, valueExpr={
    
    data <- isolate(input$Data)
    
    if (data == "") {
      
      oldpar <- par(no.readonly = TRUE)    # code line i
      on.exit(par(oldpar))            # code line i + 1 
      
      
      par(mar=c(0, 0, 0, 0))
      plot(x=c(0, 1), y=c(0, 1), axes=FALSE,
           xaxt="n", yaxt="n", main="",
           xlab = "", ylab = "",
           xaxs="i", yaxs="i", type="n")
      text(x = 0.5, y=0.6, labels = "Load data to be analyzed",
           col="red", cex = 1.6)
    } else {
      
      LDL <- isolate(input$LDL)
      UDL <- isolate(input$UDL)
      method <- isolate(input$censored)
      distribution <- isolate(input$distribution)
      n.mixture <- as.numeric(isolate(input$mixture))
      
      if (LDL == "") LDL <- NA
      LDL <- as.numeric(LDL)
      if (UDL == "") UDL <- NA
      UDL <- as.numeric(UDL)
      
      # 89.7707600916124 12.7528137253536 44.1854453769165 60.669359673936 32.6843705727645 34.5730704570828 28.090749695105 91.1180106617791 13.9854115059374 56.9528385259751 59.6424331262563 76.2201823264057 84.6428469240878 61.4539947006212 27.5032991492022 17.2751640631116 16.3951178618038 22.1066596825853 15.8104734635694 62.860136780007 86.3920347016796 50.8633719498414 26.3350623257843 29.4309413858114 53.9577371460467 40.9435118742182 64.1395543447893 19.4952849783914 40.2695044921001 44.0330480738343 66.7024428576167 22.1085943874056 21.0741015828443 36.5769530586322 70.2957796454758 50.2016515754352 54.1383590002363 13.1126129957134 22.7606870305514 40.0486960995248 28.4901935386691 LDL 13.9830972917087 LDL 62.2462103726513 55.5096384295298 15.3024423712639 11.2900250294085 14.3952333513986 40.8815827926869 26.3100410588769 33.0489156162339 19.0232074637007 47.3519350475661 56.1662979035662 24.0095554995613 18.6680724684696 20.5267938077808 61.2603324259142 LDL 15.4503639027022 41.6557100140753 77.7671619325677 71.6647052427532 22.139625491507 LDL 46.6179606969968 39.6541156998326 58.9601278755938 LDL 28.0180573785966 83.3386058631985 64.2807584194209 24.31864943644 22.5232960339106 LDL 64.4776116851174 62.776603594233 11.4504640878683 32.9139695459943 39.1124720316336 71.295614787684 56.7337043921124 67.4445163896557 32.4214974228418 45.6384034482342 32.6078013843732 60.9605529048265 16.3864347564686 LDL 12.2659233909316 143.260422377747 18.8002300449792 13.2826042818602 29.2635178927549 50.7831016585992 15.9864114168415 45.3509802007882 46.4315969264461 89.133165449313
      
      
      # UDL <- NULL
      # LDL <- 10
      # 
      # DL <- 10
      # # Generate 100 random data from a gamma distribution
      # data_x <- rgamma(100, scale=20, shape=2)
      # # remove the data below the detection limit
      # data_x[data_x < DL] <- -Inf
      
      # data <- "10 20 30 40 LDL LDL 5"
      data <- gsub("\\n", " ", data)
      data <- gsub("\\r", " ", data)
      data <- gsub("\\t", " ", data)
      data <- gsub(";", " ", data)
      data <- gsub(",", ".", data)
      data <- gsub(" +", " ", data)
      data <- gsub(" +$", "", data)
      data <- gsub("^ +", "", data)
      
      data_x <- strsplit(data, " ")[[1]]
      data_x <- ifelse(data_x == "UDL", "Inf", data_x)
      data_x <- ifelse(data_x == "LDL", "-Inf", data_x)
      data_x <- as.numeric(data_x)
      
      if (distribution == "auto") {
        
        result_normal <- cutter(observations = data_x,
                                lower_detection_limit = LDL,
                                upper_detection_limit = UDL,
                                cut_method = method,
                                distribution = "normal", 
                                n.iter=NULL, n.mixture=n.mixture, 
                                progress.bar=FALSE)
        
        result_lognormal <- cutter(observations = data_x,
                                   lower_detection_limit = LDL,
                                   upper_detection_limit = UDL,
                                   cut_method = method,
                                   distribution = "lognormal", 
                                   n.iter=NULL, n.mixture=n.mixture, 
                                   progress.bar=FALSE)
        
        result_gamma <- cutter(observations = data_x,
                               lower_detection_limit = LDL,
                               upper_detection_limit = UDL,
                               cut_method = method,
                               distribution = "gamma", 
                               n.iter=NULL, n.mixture=n.mixture, 
                               progress.bar=FALSE)
        
        result_Weibull <- cutter(observations = data_x,
                                 lower_detection_limit = LDL,
                                 upper_detection_limit = UDL,
                                 cut_method = method,
                                 distribution = "weibull", 
                                 n.iter=NULL, n.mixture=n.mixture, 
                                 progress.bar=FALSE)
        
        result_generalized.gamma <- cutter(observations = data_x,
                                           lower_detection_limit = LDL,
                                           upper_detection_limit = UDL,
                                           cut_method = method,
                                           distribution = "generalized.gamma", 
                                           n.iter=NULL, n.mixture=n.mixture, 
                                           progress.bar=FALSE)
        
        C <- compare_AICc(normal=result_normal, 
                          lognormal=result_lognormal, 
                          gamma=result_gamma, 
                          weibull=result_Weibull, 
                          generalized.gamma=result_generalized.gamma, 
                          silent = TRUE)
        
        distribution <- rownames(C)[which.min(as.numeric(C[, "AICc"]))]
        selected <- paste0("The selected distribution model is ", distribution, " with an Akaike weight based on AICc of ", 
                           specify_decimal(C[distribution, "Akaike_weight"]), ".\n")
      } else {
        selected <- ""
      }
      
      if (progressbar) {
        result <- cutter(observations = data_x,
                         lower_detection_limit = LDL,
                         upper_detection_limit = UDL,
                         cut_method = method,
                         distribution = distribution, 
                         n.mixture=n.mixture, session=session)
      } else {
        result <- cutter(observations = data_x,
                         lower_detection_limit = LDL,
                         upper_detection_limit = UDL,
                         cut_method = method,
                         distribution = distribution, 
                         n.mixture=n.mixture, session=NULL)
      }
      
      textout <- print(result, silent=TRUE)
      
      textout <- paste0(selected, textout)
      
      output$Result <- renderText(textout)
      
      
      plot(result)
      
      # plot(1, 1)
    }
  }
  )
  
  
  
  
  output$Plot <- renderPlot({
    
    outplotfn()
    
  })
  
  
  
})
