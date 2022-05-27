#' addS3Class add a S3 class to an object
#' @title Add a S3 class to an object. 
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return The same object with the new class as first class
#' @param x The object to add class.
#' @param class The class to add.
#' @description Add a S3 class as first class to an object.
#' @examples
#' print.CF <- function(x) {cat("print.CF ", x)}
#' result <- "Je suis donc je pense"
#' result <- addS3Class(result, class="CF")
#' class(result)
#' print(result)
#' result <- addS3Class(result, class=c("ECF", "OCF"))
#' class(result)
#' print(result)
#' @export

addS3Class <- function(x, class=NULL) {
  class(x) <- c(class, class(x)[!(class(x) %in% class)])
  return(x)
}
