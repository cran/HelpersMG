#' compare runs a shiny application for basic functions of comparison
#' @title Run a shiny application for basic functions of comparison
#' @author Marc Girondot
#' @return Nothing
#' @description Run a shiny application for basic functions of comparison.
#' @references Girondot, M., Guillon, J.-M., 2018. The w-value: An alternative to t- and X2 tests. Journal of Biostatistics & Biometrics 1, 1-3.
#' @family w-value functions
#' @examples
#' \dontrun{
#' library(HelpersMG)
#' compare()
#' }
#' @export


compare <- function() {
  
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("shiny package is absent; Please install it first")
  }
  
getFromNamespace("runApp", ns="shiny")(appDir = system.file("shiny", package="HelpersMG"), 
                                       launch.browser =TRUE)

}
