#' FormatCompareAIC formats data to be used with compare_AIC()
#' @title Format data to be used with compare_AIC()
#' @author Marc Girondot
#' @return An object to be used with compare_AIC()
#' @param logLik The log likelihood
#' @param nobs Number of observations
#' @param df Number of parameters
#' @description Format data to be used with compare_AIC(), compare_AICc() and compare_BIC().\cr
#' Note that logLik is supposed to not be -logLik.
#' @family AIC
#' @examples
#' \dontrun{
#' ED <- FormatCompareAIC(logLik=-140, nobs=100, df=3)
#' L <- FormatCompareAIC(logLik=-145, nobs=100, df=4)
#' compare_AIC(L=L, ED=ED)
#' compare_AICc(L=L, ED=ED)
#' compare_BIC(L=L, ED=ED)
#' }
#' @export


FormatCompareAIC <- function(logLik, nobs, df) {
  l <- list(logLik=logLik, 
            AIC=-2*logLik+2*df, 
            AICc=(-2*logLik+2*df)+(2*df*(df+1))/(nobs-df-1), 
            BIC=-2*logLik+log(nobs)*df)
  attributes(l) <- modifyList(attributes(l), list(nall=nobs , nobs=nobs, 
                        df=df, class="compareAIC"))
  return(l)
}
