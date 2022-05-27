#' RM_duplicate duplicates a results managment within an object
#' @title Duplicate a results managment within an object.
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return The original object with a duplicated results managment.
#' @param x The object to duplicate a results managment
#' @param RM The name of results managment stored
#' @param RMnamefrom The name of the results managment to be duplicated
#' @param RMnameto The new name of the results managment
#' @description RM_duplicate duplicates a results managment within an object.
#' @family Results Managment
#' @examples
#' \dontrun{
#' library("HelpersMG")
#' # Let an object of class objclass being created
#' obj <- list(A=100, name="My object")
#' class(obj) <- "objclass"
#' # And now I create a RM to this object
#' obj <- RM_add(x=obj, RMname="NewAnalysis1")
#' RM_list(obj)
#' obj <- RM_duplicate(x=obj, RMnamefrom="NewAnalysis1", RMnameto="NewAnalysis2")
#' RM_list(obj)
#' }
#' @export

RM_duplicate <- function (x=stop("An object with results managment must be provided"), 
                     RM = "RM", RMnamefrom=1, RMnameto=2) {
  
  xa <- attributes(x)
  rx <- xa[[RM]][RMnamefrom]
  rxnew <- rx
  rxnew[[RMnamefrom]]$name <- as.character(RMnameto)
  names(rxnew) <- as.character(RMnameto)
  rxnew[[RMnameto]]$timestamp <- timestamp(quiet = TRUE)

  ix <- modifyList(xa[[RM]], rxnew)
  
  rx <- list(RM=ix)
  names(rx) <- RM
  
  attributes(x) <- modifyList(xa, rx)
  return(x)
}
