#' moon_phase calculates the moon phase based on a date.
#' @title Moon phase based on a date
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @return Return a value describing the moon phase:\cr
#' 0 and 100 are full moon, 50 is new moon, 25 last quarter and 75 first quater
#' @param date A date in known format
#' @param phase If TRUE, a vector of character with NM, FQ, FL LQ will be returned
#' @description The script gives an index (base 100) that represents moon phase. If the value lays between:\cr
#' 0 and 1.6931595 or 98.3068405 and 100, it is full moon,\cr
#' 23.3068405 and 26.6931595, last quarter,\cr
#' 48.3068405 and 51.6931595, new moon,\cr	
#' 73.3068405 and 76.6931595, first quarter\cr
#' When phase is set to TRUE, a character representing the moon phase is returned.
#' @keywords Moon Lunar Lune
#' @examples 
#' library("HelpersMG")
#' moon_phase(as.Date("2001-12-31"))
#' moon_phase(as.Date("14/04/2010", "%d/%m/%Y"))
#' moon_phase(as.Date("22/06/07", "%d/%m/%y"))
#' moon_phase(seq(from=as.Date("2012-03-01"), 
#'		to=as.Date("2012-04-15"), by="days"))
#' moon_phase(seq(from=as.Date("2012-03-01"), 
#' 		to=as.Date("2012-04-15"), by="days"), phase=TRUE)
#' @export

moon_phase<- function(date=NULL, phase=FALSE) {
if (!is.null(date)) {

XRef<-as.Date("2012-03-07")

nb<-as.numeric(date-XRef)

moon<-nb%%29.530589

moon<-100*(moon/29.530589)

if (phase) {

moonT<-rep(NA, length(moon))
moonT[which((moon>98.3068405) | (moon<1.6931595))]<-"FM"
moonT[which((moon>48.3068405) & (moon<51.6931595))]<-"NM"
moonT[which((moon>23.3068405) & (moon<26.6931595))]<-"LQ"
moonT[which((moon>73.3068405) & (moon<76.6931595))]<-"FQ"

return(moonT)

} else {

return(moon)

}

}
}
