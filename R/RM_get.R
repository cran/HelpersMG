#' RM_get gets a value in results managment to an object
#' @title Get a value in a results managment to an object
#' @author Marc Girondot
#' @return Return a value in a results managment object
#' @param x The object in which to get a result in a results managment
#' @param RM The name of results managment stored
#' @param RMname The name of the results managment to be read
#' @param valuename The name of the value to be read
#' @description Return the value valuename of the results managment RMname.
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
#' }
#' @export

RM_get <- function(x=stop("An object with results managment must be provided"), 
                   RM="RM", 
                   RMname=stop("A results managment name must be provided"), 
                   valuename=NULL) {
  xa <- attributes(x)
  # dans rx j'ai le RM
  rx <- xa[[RM]][[RMname]][[valuename]]
  return(rx)
}
