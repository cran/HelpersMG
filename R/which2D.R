#' which2D returns the rows and columns at which a value is observed in a data.frame or a matrix
#' @title Returns the rows and columns at which a value is observed
#' @author Marc Girondot
#' @return A matrix with the rows and columns
#' @param x A data.frame or a matrix
#' @param type Can be min, max or =, ==, >, <, >=, <=, <> or !=
#' @param objective If type is different from min or max, gives the value for comparison
#' @param by Does the search must be done by rows or by columns
#' @description This function return the rows and columns of the maximum or minimum or a specific value is observed.
#' @examples
#' library("HelpersMG")
#' e <- data.frame(L=c(20, 20, 10), M=c(60, 10, 20))
#' which2D(e)
#' which2D(e, by="row")
#' which2D(e, type="max")
#' which2D(e, type="=", objective=20)
#' which2D(e, type="<", objective=50)
#' which2D(e, type="=", objective=20, by="row")
#' @export

which2D <- function(x, type="min", objective=NULL, by="column") {
  if (mode(x)=="numeric" | mode(x)=="logical") return(which.min(x))
  if (class(x)=="data.frame") x <- as.matrix(x)
  
  if (by=="row") x <- t(x)
  Dv <- as.vector(x, mode="any")
  if (type=="min")   pos=which(Dv==min(Dv))-1  else   
    if (type=="max") pos=which(Dv==max(Dv))-1 else 
      if (type=="=" | type=="==") pos=which(Dv==objective) else
        if (type==">") pos=which(Dv>objective) else
          if (type==">=") pos=which(Dv>=objective) else 
            if (type=="<") pos=which(Dv<objective) else
              if (type=="<=") pos=which(Dv<=objective) else
                if (type=="!=" | type=="<>") pos=which(Dv!=objective)
  r <- t(sapply(pos, function(y) c(y%%nrow(x)+1, floor(y/nrow(x))+1)))
  if (by=="row") r <- r[, c(2,1)]
  colnames(r) <- c("rows", "columns")
  return(r)
}

