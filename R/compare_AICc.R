#' compare_AICc compares the AICc of several outputs obtained with the same data.
#' @title Compares the AICc of several outputs
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return A list with DeltaAICc and Akaike weight for the models.
#' @param ... Successive results to be compared as lists.
#' @param factor.value The $value of the list object is multiplied by factor.value to calculate BIC.
#' @param silent If TRUE, nothing is displayed.
#' @param FUN Function used to show values
#' @description This function is used to compare the AICc of several outputs obtained with the same data but with different set of parameters.\cr
#' Each object must have associated \code{logLik()} method with df and nobs attributes.\cr
#' AICc for object x will be calculated as \code{2*factor.value*logLik(x)+(2*attributes(logLik(x))$df*(attributes(logLik(x))$df+1)/(attributes(logLik(x))$nobs-attributes(logLik(x))$df-1)}.\cr
#' @family AIC
#' @examples
#' \dontrun{
#' library("HelpersMG")
#' # Here two different models are fitted
#' x <- 1:30
#' y <- rnorm(30, 10, 2)+log(x)
#' plot(x, y)
#' d <- data.frame(x=x, y=y)
#' m1 <- lm(y ~ x, data=d)
#' m2 <- lm(y ~ log(x), data=d)
#' compare_BIC(linear=m1, log=m2, factor.value=-1)
#' # Here test if two datasets can be modeled with a single model
#' x2 <- 1:30
#' y2 <- rnorm(30, 15, 2)+log(x2)
#' plot(x, y, ylim=c(5, 25))
#' plot_add(x2, y2, col="red")
#' d2 <- data.frame(x=x2, y=y2)
#' m1_2 <- lm(y ~ x, data=d2)
#' x_grouped <- c(x, x2)
#' y_grouped <- c(y, y2)
#' d_grouped <- data.frame(x=x_grouped, y=y_grouped)
#' m1_grouped <- lm(y ~ x, data=d_grouped)
#' compare_AICc(separate=list(m1, m1_2), grouped=m1_grouped, factor.value=-1)
#' # Or simply
#' compare_AICc(m1=list(AICc=100), m2=list(AICc=102))
#' }
#' @export


compare_AICc <- function(..., factor.value=-1, silent=FALSE, 
                         FUN=function(x) specify_decimal(x, decimals=2)) {

  # factor.value=-1
  # silent=FALSE 
  # FUN=function(x) specify_decimal(x, decimals=2)
  # result <- list()
  
  result <- list(...)
  
  if (is.list(result) & length(result)==1) result <- unlist(result, recursive=FALSE)
  
  if (!is.null(result)) {
    if (!is.list(result) || (is.null(names(result))) || (any(names(result)==""))) {
      stop("The results must be included within a list with names; see example.")
    } else {
      out<-NULL
      
      l <- length(result)
      
      if (l<2) {
        stop("A least two results must be provided.")
      } else {
        aicc <- NULL
        name <- names(result)
        for (i in 1:l) {
          
          encours <- result[i]
          # Je vais être un peu plus strict que pour l'AIC
          # Il faut que logLik(encours[[1]]) existe avec les attributs corrects
          # nall, nobs (nb d'observations) et df (nb de parametres)
          
          if (!is.null(encours[[1]]$AICc)) {
            aicc <- c(aicc, encours[[1]]$AICc)
          } else {
            
            aaec <- try(logLik(encours[[1]]), silent=TRUE)
           t <- (inherits(aaec, "try-error"))
          if (t) encours <- encours[[1]]
          
          sumL <- 0
          sumdf <- 0
          sumnobs <- 0

          for (j in 1:length(encours)) {
            encours2 <- encours[[j]]
            aaec <- try(logLik(encours2), silent=TRUE)
            t <- (inherits(aaec, "try-error"))
            if (t) {
              stop(paste("Object", name[i], "has not the required format"))
            }
            sumL <- sumL + logLik(encours2)
            sumdf <- sumdf + attributes(logLik(encours2))$df
            sumnobs <- sumnobs + attributes(logLik(encours2))$nobs
          }
          aicc <- c(aicc, 2*factor.value*sumL+2*sumdf+(2*sumdf*(sumdf+1))/(sumnobs-sumdf-1))
          }
          
        }
        
        bestaicc <- min(aicc)
        ser<-which.min(aicc)
        deltaaicc<-aicc-bestaicc
        aw<-exp(-0.5*deltaaicc)
        saw=sum(aw)
        aw<-aw/saw
        
        out<-data.frame(cbind(AICc=FUN(aicc), DeltaAICc=FUN(deltaaicc), Akaike_weight=FUN(aw)), row.names=name)
        if (!silent) print(paste("The lowest AICc (",sprintf("%.3f", bestaicc) ,") is for series ", name[ser], " with Akaike weight=", sprintf("%.3f", aw[ser]), sep=""))
        
        return(out)
      }
    }
  }
}
