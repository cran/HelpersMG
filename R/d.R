#' d Write an ASCII Representation of an Object
#' @title Write an ASCII Representation of an Object
#' @author Marc Girondot
#' @return A string
#' @param x A named vector object
#' @param file either a character string naming a file or a connection. "" indicates output to the console.
#' @param control character vector indicating deparsing options. See .deparseOpts for their description.
#' @param collapse Characters used to separate values.
#' @description Writes an ASCII text representation of an R object.\cr
#' It can be used as a replacement of dput() for named vectors.\cr
#' The controls "keepNA", "keepInteger" and "showAttributes" are utilized for named vectors.
#' @family Characters
#' @examples
#' d(c(A=10, B=20))
#' dput(c(A=10, B=20))
#' @export

d <- function(x, file = "",
              control = c("keepNA", "keepInteger", "showAttributes"), 
              collapse=", \n  ") {
  
  opts <- intToBits(.deparseOpts(control))
  # 1 (1): keepInteger
  # 4 (3): showAttributes
  # 64 (7): keepNA
  
  if ((any(names(attributes(x)) != "names") & !is.null(attributes(x))) | 
    (opts[3] == 0) | 
       ((class(x) != "numeric") & (class(x) != "character") & (class(x) != "integer"))) {
    dput(x, file=file, control=control)
  } else {
    if (is.null(names(x))) {
      dput(x, file=file, control=control)
    } else {
      if (is.character(file)) 
        if (nzchar(file)) {
          file <- file(file, "wt")
          on.exit(close(file))
        }
      else file <- stdout()
      
      
      if (is.numeric(x)) {
        cat("c(", paste0(sapply(seq_along(x), FUN=function(i) {
          x1 <- x[i]
          ifelse(is.integer(x1) & (opts[1]==1),
                 if (is.na(x1) & (opts[7]==1)) {
                   paste0(names(x1), " = NA_integer_")
                 } else {
                 paste0(names(x1), " = ", x1, "L")
                   }, 
                 if (is.na(x1) & (opts[7]==1)) {
                   paste0(names(x1), " = NA_real_")
                 } else {
                 paste0(names(x1), " = ", format(x1, digits = 17, trim = TRUE))
                 }
          )
          }
        ), 
        collapse=collapse), ")", sep="", file = file)
      } else {
        cat("c(", paste0(names(x), " = '", 
                         x, "'", collapse=collapse), ")", sep="", file = file)
      }
    }
  }
}