#' duplicated_packages lists the duplicated packages with their locations
#' @title List the duplicated packages with their locations
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return A data.frame with 4 elements for each duplicated packages:\cr
#' - versions: the version of the packages\cr
#' - libraries: the locations\cr
#' @description A data.frame with the duplicated packages and their locations and version.\cr
#' The columns Lib1 and Version1 should have the oldest version of the packages.
#' @examples
#' \dontrun{
#' library(HelpersMG)
#' duplicated_packages()
#' # To remove the oldest versions of the installed packages, use
#' li <- duplicated_packages()
#' if (nrow(li) != 0)
#'     for (i in 1:nrow(li))
#'         remove.packages(rownames(li)[i], lib=li[i, "Lib1"])
#' }
#' @export


duplicated_packages <- function() {
  li <- matrix(character(0), ncol=4, 
               dimnames = list(NULL, 
                               c("Lib1", "Version1", "Lib2", "Version2")))
  lp <- list.packages()
  nli <- names(lp)
  if (length(nli) > 1) {
    for (i in 1:(length(nli)-1)) {
      for (j in (i+1):length(nli)) {
        
        d1 <- na.omit(match(rownames(lp[[j]]), rownames(lp[[i]])))
        # print(length(d1))
        if (length(d1) != 0) {
          for (k in seq_along(d1)) {
          m <- rownames(lp[[i]])[d1[k]]
          vi <- lp[[i]][d1[k], "Version"]
          vj <- lp[[j]][m, "Version"]
          if (compareVersion(vi, vj)==-1) {
          li <- rbind(li, matrix(c(nli[i], vi, nli[j], vj), ncol=4, 
                                 dimnames = list(m, 
                                                 c("Lib1", "Version1", "Lib2", "Version2"))))
          } else {
          li <- rbind(li, matrix(c(nli[j], vj, nli[i], vi), ncol=4, 
                                   dimnames = list(m, 
                                                   c("Lib1", "Version1", "Lib2", "Version2"))))
            
          }
          }
        }
      }
    }
  }
  return(li)
}
