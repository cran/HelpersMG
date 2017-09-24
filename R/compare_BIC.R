#' compare_BIC compares the BIC of several outputs obtained with the same data.
#' @title Compares the BIC of several outputs
#' @author Marc Girondot
#' @return A list with DeltaBIC and Akaike weight for the models.
#' @param ... Successive results to be compared as lists.
#' @param factor.value The $value of the list object is multiplied by factor.value to calculate BIC.
#' @description This function is used to compare the BIC of several outputs obtained with the same data but with different set of parameters.\cr
#' Each object must have associated \code{logLik()} method with df and nobs attributes.\cr
#' BIC for object x will be calculated as \code{2*factor.value*sum(logLik(x))+sum(attributes(logLik(x))$df)*log(attributes(logLik(x))$nobs))}.\cr
#' When several data (i..n) are included, the global BIC is calculated as:\cr
#' \code{2*factor.value*sum(logLik(x)) for i..n+sum(attributes(logLik(x))$df) for i..n*log(attributes(logLik(x))$nobs for i..n))}
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
#' compare_BIC(separate=list(m1, m1_2), grouped=m1_grouped, factor.value=-1)
#' }
#' @export


compare_BIC <- function(..., factor.value=-1) {

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
        bic <- NULL
        name <- names(result)
        for (i in 1:l) {
          
          encours <- result[i]
          # Je vais être un peu plus strict que pour l'AIC
          # Il faut que logLik(encours[[1]]) existe avec les attributs corrects
          # nall, nobs (nb d'observations) et df (nb de parametres)
          
          t <- (class(try(logLik(encours[[1]]), silent=TRUE))=="try-error")
          if (t) encours <- encours[[1]]
          
          sumL <- 0
          sumdf <- 0
          sumnobs <- 0

          for (j in 1:length(encours)) {
            encours2 <- encours[[j]]
            t <- (class(try(logLik(encours2), silent=TRUE))=="try-error")
            if (t) {
              stop(paste("Object", name[i], "has not the required format"))
            }
            sumL <- sumL + logLik(encours2)
            sumdf <- sumdf + attributes(logLik(encours2))$df
            sumnobs <- sumnobs + attributes(logLik(encours2))$nobs
          }
          bic <- c(bic, 2*factor.value*sumL+sumdf*log(sumnobs))	
        }
        
        bestbic<-min(bic)
        ser<-which.min(bic)
        deltabic<-bic-bestbic
        aw<-exp(-0.5*deltabic)
        saw=sum(aw)
        aw<-aw/saw
        
        out<-data.frame(cbind(BIC=bic, DeltaBIC=deltabic, Akaike_weight=aw), row.names=name)
        print(paste("The lowest BIC (",sprintf("%.3f", bestbic) ,") is for series ", name[ser], " with Akaike weight=", sprintf("%.3f", aw[ser]), sep=""))
        
        return(out)
      }
    }
  }
}
