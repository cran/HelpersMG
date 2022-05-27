#' RM_delete deletes a results managment or a result within a results managment from an object
#' @title Delete a results managment or a result within a results managment from an object
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return The original object with the deleted results managment
#' @param x The object to delete a results managment
#' @param RM The name of results managment stored
#' @param RMname The name of the result that will be deleted or its rank
#' @param valuename The name of the result that will be deleted
#' @description Return the original object with the deleted results managment or result.
#' @family Results Managment
#' @examples
#' \dontrun{
#' library("HelpersMG")
#' # Let an object of class objclass being created
#' obj <- list(A=100, name="My object")
#' class(obj) <- "objclass"
#' # And now I create a RM to this object
#' obj <- RM_add(x=obj, RMname="NewAnalysis1")
#' obj <- RM_add(x=obj, RMname="NewAnalysis2")
#' RM_list(obj)
#' obj <- RM_delete(x=obj, RMname="NewAnalysis1")
#' RM_list(obj)
#' obj <- RM_delete(x=obj, RMname=1)
#' RM_list(obj)
#' obj <- RM_add(x=obj, RMname="NewAnalysis1", valuename="V1", value=100)
#' RM_list(obj)
#' RM_get(x=obj, RMname="NewAnalysis1", valuename="V1")
#' obj <- RM_add(x=obj, RMname="NewAnalysis1", valuename="V2", value=200)
#' RM_get(x=obj, RMname="NewAnalysis1", valuename="V2")
#' obj <- RM_delete(x=obj, RMname="NewAnalysis1", valuename="V1")
#' RM_get(x=obj, RMname="NewAnalysis1", valuename="V1")
#' RM_get(x=obj, RMname="NewAnalysis1", valuename="V2")
#' }
#' @export

RM_delete <- function(x=stop("An object with results managment must be provided"), 
                      RM="RM", 
                      RMname=stop("A name must be provided"), 
                      valuename=NULL) {
  xa <- attributes(x)
  rx <- xa[[RM]]
  if (is.null(rx)) {
    # Il n'y avait pas d'object RM
    return(x)
  }
  
  # Si le name est un numéro, je prends la liste et le nom
  if (is.numeric(RMname)) {
    RMname <- names(rx)[RMname]
  }
  
  # il y en avait un
  irx <- which(names(rx)==RMname)
  # Si le nom n'existe pas
  if (identical(irx, integer(0))) {
    return(x)
  }
  
  if (is.null(valuename)) {
    #  Je retire tout le RNname
    rx <- rx[-irx]
  } else {
    if (!identical(integer(0), which(names(rx[[irx]]) == valuename))) {
      rx[[irx]] <- rx[[irx]][-which(names(rx[[irx]]) == valuename)]
    }
  }
  xa[[RM]] <- rx
  attributes(x) <- xa
  
  return(x)
  
}
