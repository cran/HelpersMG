#' RM_list returns the list of results managment of an object
#' @title Return the list of results managment of an object.
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return A list with the names of results stored in an object
#' @param x The object to add a results managment
#' @param RM The name of results managment stored
#' @param silent Should the results be shown ?
#' @param max.level If TRUE, will return all list element of the objects
#' @description RM_list returns the list of results managment of an object.
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
#' obj <- RM_add(x=obj, RMname="NewAnalysis2")
#' RM_list(obj)
#' obj <- RM_add(x=obj, RMname="NewAnalysis2", valuename="V1", value=100)
#' RM_get(x=obj, RMname="NewAnalysis2", valuename="V1")
#' obj <- RM_add(x=obj, RMname="NewAnalysis2", valuename="V1", value=200)
#' RM_get(x=obj, RMname="NewAnalysis2", valuename="V1")
#' obj <- RM_add(x=obj, RMname="NewAnalysis2", valuename="V2", value=300)
#' RM_get(x=obj, RMname="NewAnalysis2", valuename="V2")
#' RM_list(obj)
#' rmlist <- RM_list(obj, max.level=TRUE)
#' rmlist
#' }
#' @export

RM_list <- function (x=stop("An object with results managment must be provided"), 
                     RM = "RM", silent=FALSE, max.level = FALSE) {
  x <- attributes(x)[[RM]]
  rx <- list()
  for (i in names(x)) {
    ix <- list(name = i, timestamp = x[[i]]$timestamp)
    if ((length(x[[i]]) > 2) & max.level) {
      jx <- as.list(names(x[[i]])[3:length(x[[i]])])
      names(jx) <- names(x[[i]])[3:length(x[[i]])]
      ix <- c(ix, jx)
    }
    rx <- c(rx, list(ix))
    if (!silent) cat(i, " at ", x[[i]]$timestamp, "\n")
  }
  names(rx) <- names(x)
  return(invisible(rx))
}
