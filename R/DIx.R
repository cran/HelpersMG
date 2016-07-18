#' DIx returns an index of quantitative asymmetry and complexity
#' @title Return an index of quantitative asymmetry and complexity named Developmental Instability Index (DIx)
#' @author Marc Girondot
#' @return A numeric value
#' @param l1 Set of measures at one side of an organism
#' @param l2 Set of measures at the other side of an organism
#' @param details If TRUE, will show the details of computing
#' @param version Can be 1 or 2; see description
#' @description Return an index of quantitative asymmetry and complexity.\cr
#' Higher is the value, higher is the complexity (number of objects) and diversity (difference between them).\cr
#' The indice is based on the product of the average angular distance of 
#' Edwards (1971) for all permutations of measures for both sides with the 
#' geometric mean of the inverse of Shannon entropy H for both sides.
#' Let p1 and p2 two vectors of relative measures of objects with sum(p1) = 1 and sum(p2)=1 
#' and n1 being the number of objects in p1 and n2 being the number of objects in p2.\cr
#' Edwards distance for all permutations of p1 and p2 objects are computed and the average value E is calculated.\cr
#' The maximun possible Shannon index for identical n1 is max1 = sum((1/n1) * log(1/n1)).\cr
#' Shannon index is v1 = sum(p1 * log(p1)).\cr
#' If version == 2, the complementary of Shannon index for these n1 objects is used: c1 = 2 * max1 - v1\cr
#' If version == 1, the Shannon index is used directly.\cr
#' The geometry mean between both sides defined the measure of diversity within each side: S=sqrt(c1 * c2)\cr
#' The Developmental Instability Index is then S * E
#' @references 
#' Edwards, A.W.F., 1971. Distances between populations on the basis of gene frequencies. Biometrics 27, 873–881.\cr
#' Shannon C.E. 1948 A mathematical theory of communication. Bell System Technical Journal 27(3), 379-423.\cr
#' @examples
#' l1 <- c(0.1, 0.1, 0.05, 0.2, 0.3, 0.25)
#' l2 <- c(0.2, 0.3, 0.5)
#' DIx(l1, l2)
#' 
#' l1 <- c(0.1, 0.1, 0.05, 0.2, 0.3, 0.25)
#' l2 <- c(0.1, 0.1, 0.05, 0.2, 0.3, 0.25)
#' DIx(l1, l2)
#' 
#' l1 <- c(0.2, 0.3, 0.5)
#' l2 <- c(0.2, 0.3, 0.5)
#' DIx(l1, l2)
#' 
#' l1 <- c(0.2, 0.2, 0.2, 0.2, 0.2)
#' l2 <- c(0.2, 0.3, 0.5)
#' DIx(l1, l2)
#' 
#' l1 <- c(0.2, 0.2, 0.2, 0.2, 0.2)
#' l2 <- c(0.3333, 0.3333, 0.3333)
#' DIx(l1, l2)
#' 
#' l1 <- c(0.2, 0.2, 0.2, 0.2, 0.2)
#' l2 <- c(0.2, 0.2, 0.2, 0.2, 0.2)
#' DIx(l1, l2)
#' 
#' l1 <- c(0.3333, 0.3333, 0.3333)
#' l2 <- c(0.3333, 0.3333, 0.3333)
#' DIx(l1, l2)
#' @export

DIx <- function(l1, l2, details=FALSE, version = 1) {
  l1 <- c(l1, rep(0, max(0, length(l2)-length(l1))))
  l2 <- c(l2, rep(0, max(0, length(l1)-length(l2))))
  
  l1 <- l1/sum(l1)
  l2 <- l2/sum(l2)
  
  pp <- getFromNamespace(".permutations", ns="HelpersMG")(v=l1, r=length(l1), n=length(l1), 
                                                          set=FALSE)
  p <- (mean(apply(pp, 1, function(x) {sqrt(1-(sum(x*l2)))})))
  vn1 <- sum(l1 != 0)
  vn2 <- sum(l2 != 0)
  if (version == 2) {
    vmax1 <- -sum(rep(1/vn1, vn1) * log(rep(1/vn1, vn1)))
    vmax2 <- -sum(rep(1/vn2, vn2) * log(rep(1/vn2, vn2)))
    
    p1 <- 2*vmax1+sum(l1[l1 != 0] * log(l1[l1 != 0]))
    p2 <- 2*vmax2+sum(l2[l2 != 0] * log(l2[l2 != 0]))
  } else {
    p1 <- sum(l1[l1 != 0] * log(l1[l1 != 0]))
    p2 <- sum(l2[l2 != 0] * log(l2[l2 != 0]))
  }
  if (details) {
    return(c(Average.Edwards=p, Geometric.mean.Shannon=sqrt(p1*p2), DIx=p*sqrt(p1*p2)))
  } else {
    return(p*sqrt(p1*p2))
  }
}



