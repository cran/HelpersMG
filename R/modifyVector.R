#' modifyVector modifies elements of a vector
#' @title Modifies Elements of a Vector
#' @author Marc Girondot
#' @return A modified version of x, with the elements of val replacing the elements of x
#' @param x A named vector.
#' @param val A named vector with components to replace corresponding components in x.
#' @param add If FALSE, only existing elements of x are returned.
#' @description Modifies a vector by changing a subset of elements to match a second vector.
#' @examples
#' library("HelpersMG")
#' e <- c(M=10, L=20, J=30)
#' modifyVector(e, c(U=10, M=30))
#' modifyVector(e, c(U=10, M=30), add=FALSE)
#' @export

modifyVector <- function(x, val, add=TRUE) {
  r <- unlist(modifyList(as.list(x), as.list(val)))
  if (!add) {
    r <- r[names(x)]
  }
  return(r)
}

