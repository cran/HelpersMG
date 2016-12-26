#' compare runs a shiny application for basic functions of comparison
#' @title Run a shiny application for basic functions of comparison
#' @author Marc Girondot
#' @return Nothing
#' @description Run a shiny application for basic functions of comparison
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
