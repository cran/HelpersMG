#' tnirp reads an ASCII text representation of a named or not vector object
#' @title Read an ASCII text representation of a named or not vector object
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return A vector
#' @param x A string or a vector of strings with value and possibly names.
#' @param named TRUE if names are included.
#' @description Read an ASCII text representation of a named or not vector object.\cr
#' Note that paste0(rev(c("p", "r", "i", "n", "t")), collapse="") = "tnirp"
#' @family Characters
#' @examples
#' A <- structure(runif(26), .Names=letters)
#' text <- capture.output(A)
#' tnirp(text)
#' 
#' tnirp("         mu   mu_season         OTN       p1.09       p1.10       p1.11 
#'  4.63215947 10.78627511  0.36108497  0.08292101 -0.52558196 -0.76430859 
#'        p1.12       p1.13       p1.14       p1.15       p1.16       p1.17 
#'        -0.75186542 -0.57632291 -0.58017174 -0.57048696 -0.56234135 -0.80645122 
#'        p1.18       p1.19       p1.20       p1.21       p1.22       p1.23 
#'        -0.77752524 -0.80909494 -0.56920540 -0.55317302  0.45757298 -0.64155368 
#'        p1.24       p1.25       p1.26       p1.27       p1.28       p1.29 
#'        -0.59119637 -0.66006794 -0.66582399 -0.66772684 -0.67351412 -0.66941992 
#'        p1.30       p1.31       p1.32       p1.33       p1.34       p1.35 
#'        -0.67038245 -0.68938726 -0.68889078 -0.68779016 -0.68604629 -0.68361820 
#'        p1.36       p1.37       p2.09       p2.10       p2.11       p2.12 
#'        -0.67045238 -0.66115613  2.55403149  2.31060620  2.31348160  2.20958757 
#'        p2.13       p2.14       p2.15       p2.16       p2.17       p2.18 
#'        2.14304918  2.19699719  2.30705457  2.18740019  2.32305811  2.31668302 
#'        p2.19       p2.20       p2.21       p2.22       p2.23       p2.24 
#'        1.99424288  2.06613445  2.38092301  2.40551276  2.31987342  2.30344402 
#'        p2.25       p2.26       p2.27       p2.28       p2.29       p2.30 
#'        2.26869058  2.25008836  2.23385204  2.22768782  2.25341904  1.77043360 
#'        p2.31       p2.32       p2.33       p2.34       p2.35       p2.36 
#'        2.21606813  2.21581431  2.21153872  2.21118013  2.21375660  2.21182196 
#'        p2.37 
#'        1.86137833 ")
#' tnirp("   27.89 289.99
#'       90.56", named=FALSE)
#' @export

tnirp <- function(x, named=TRUE) {
  if (length(x) != 1) x <- paste0(x, collapse = "\n")
  par2 <- gsub(" +", " ", x)
  par3 <- strsplit(par2, "\\n")[[1]]
  if (named) {
  nam <- par3[seq(from=1, to=length(par3), by=2)]
  nam <- gsub("^ ", "", nam)
  nam <- gsub(" $", "", nam)
  nam <- strsplit(paste0(nam, collapse = " "), " ")[[1]]
  }
  val <- par3[seq(from=ifelse(named, 2, 1), to=length(par3), by=ifelse(named, 2, 1))]
  val <- gsub("^ ", "", val)
  val <- gsub(" $", "", val)
  val <- as.numeric(strsplit(paste0(val, collapse = " "), " ")[[1]])
  if (named) {
  names(val) <- nam
  }
  return(val)
}

