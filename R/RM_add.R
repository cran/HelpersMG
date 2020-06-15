#' RM_add adds a results managment or a value in results managment to an object
#' @title Create a results managment or add a value in a results managment to an object
#' @author Marc Girondot
#' @return The original object with a new value in a results managment object or a new results managment
#' @param x The object to add a results managment or a result in a results managment
#' @param RM The name of results managment stored
#' @param RMname The name of the results managment to be modified or created
#' @param valuename The name of the new value to be added
#' @param value The value to be added
#' @description Return original object with a new value or a new results managment.
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
#' }
#' @export

RM_add <- function(x=stop("An object with results managment must be provided"), 
                   RM="RM", 
                   RMname=stop("A results managment name must be provided"), 
                   valuename=NULL, 
                   value=NULL) {
  xa <- attributes(x)
  # dans rx j'ai le RM
  rx <- xa[[RM]]
  ix <- rx
  
  if (is.null(rx)) {
    # Il n'y avait pas d'object RM
    ix <- list(list(name=RMname, timestamp=timestamp(quiet = TRUE)))
    names(ix) <- RMname
  } else {
    # il y en avait un
    # Je le crée s'il n'existait pas
    if (!any(names(rx) %in% RMname)) {
      ix <- list(list(name=RMname, timestamp=timestamp(quiet = TRUE)))
      names(ix) <- RMname
      ix <- c(rx, ix)
    }
  }
  
  # Est ce que je dois ajouter une valeur
  if ((!is.null(valuename)) & !is.null(value)) {
    if (!any(names(ix[[RMname]]) %in% valuename)) {
      # Le nom n'existait pas
      newvalue <- list(value)
      names(newvalue) <- valuename
      ix[RMname] <- list(modifyList(ix[[RMname]], newvalue))
    } else {
      ix[[RMname]][[valuename]] <- value  
    }
  }
  
  rx <- list(RM=ix)
  names(rx) <- RM
  
  xa[RM] <- rx
  
  attributes(x) <- xa
  return(x)
}
