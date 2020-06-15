#' tide.info gets the annual tide calendar for one particular location.
#' @title Annual tide calendar for one particular location
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @return Return a data.frame with tide calendar:\cr
#' Level is the tide level, Tide is the High or Low Tide information and Date.Time is
#' the date/time in POSIXlt format.
#' @param file An html file from the site \code{http://tides.mobilegeographics.com/}
#' @param location Code based on \code{http://tides.mobilegeographics.com/}
#' @param year Year to get the calendar
#' @param latitude The latitude of the tide information
#' @param longitude The longitude of the tide information
#' @param tz Timezone
#' @family Periodic patterns of indices
#' @description The script extracts tide information from:\cr
#' \code{http://tides.mobilegeographics.com/} into a data.frame.\cr
#' The presence of XLM package is required for this function.\cr
#' @keywords Tide
#' @examples
#' \dontrun{
#' library("HelpersMG")
#' lat <- 5.74
#' long <- -54
#' Awala2004 <- tide.info(year=2004, longitude=long, latitude=lat, tz="America/Cayenne")
#' with(Awala2004, plot(Date.Time, Level, bty="n", las=1, type="l", 
#' xlab=paste("Year", as.POSIXlt(Date.Time[1])$year+1900), 
#' ylab="Tide level in m"))
#' latitude <- 45.3667
#' longitude <- -64.3833
#' NovaScotia <- tide.info(year=2018, longitude=longitude, 
#'                        latitude=latitude, tz="America/Halifax" )
#' with(NovaScotia, plot(Date.Time, Level, bty="n", las=1, type="l", xaxt="n", 
#'    xlab=paste("January", as.POSIXlt(Date.Time[1])$year+1900), 
#'    xlim=as.numeric(as.POSIXlt(as.Date(c("2018-01-01", "2018-01-31")))), 
#'    ylab="Tide level in m"))
#' axis(1, at=as.numeric(as.POSIXlt(seq(from=as.Date("2018-01-01"), 
#'                                      to=as.Date("2018-01-31"), by="1 day"))), 
#' labels=1:31)
#' segments(x0=as.numeric(as.POSIXlt(as.Date("2018-01-01"))), 
#' x1=as.numeric(as.POSIXlt(as.Date("2018-01-31"))), y0=0, y1=0, lty=2)
#' }
#' @export

tide.info <- function(file=NULL, year=as.POSIXlt(Sys.time())$year+1900, 
                      location=0, latitude=NA, longitude=NA, tz="") {
  
  
  if (!requireNamespace("XML", quietly = TRUE)) {
    stop("XML package is necessary for this function")
  }
  
  # file=NULL; year=as.POSIXlt(Sys.time())$year+1900;location=0;latitude=NA;longitude=NA;tz=""
  # longitude=4.01; latitude=6.4
  if (!is.na(latitude) & !is.na(longitude)) {
    # The file locationTide is available only within the functions of the packages
    locationTide <- locationTide
    location <- locationTide$location[which.min(abs(locationTide$latitude-latitude)+abs(locationTide$longitude-longitude))][1]
  }

  if (is.null(file)) {
    theurl <- paste0("http://tides.mobilegeographics.com/calendar/year/", location, ".html?y=", year, "&m=1&d=1")
  } else {
    theurl <- file
  }
  
  dest <- file.path(tempdir(), "tide.html")
  download.file(theurl, dest, quiet = TRUE)
  
  a <- readLines(con=dest)
  b <- sapply(a, FUN = function(x) return(gsub("&minus;", "-", x)))
  
  tables <- XML::readHTMLTable(doc=b, header=TRUE, stringsAsFactors = FALSE)
  
  # tables <- XML::readHTMLTable(theurl, stringsAsFactors=FALSE)
  n.rows <- lapply(tables, function(t) dim(t)[1])
  n.rows <- lapply(n.rows, function(t) {ifelse(is.null(t), 1, t)})
  n.rows <- unlist(n.rows)
  
  tl <- Sys.getlocale(category = "LC_TIME")
  Sys.setlocale(category = "LC_TIME", locale = gsub(".._..(.+)", "en_US\\1", tl) )
  Tide.Calendar <- data.frame(Date.time=as.POSIXlt(character(0)), 
                              Level=numeric(0), Tide=character(0))
  months <- which(n.rows>27)
  for (month in 1:12) {
    # print(month)
    table <- tables[[months[month]]]
    Date <- paste0(year,"-", ifelse(month<10, "0", ""), as.character(month), "-", gsub(" ", "0", gsub("... (..)", "\\1", table[,1])))
    # Read morning high tide
    for (col in 2:6) {
      # print (col)
      Time <- gsub("([0-9]+:[0-9]+ [AP]M) .+", "\\1", table[,col])
      
      Date.Time <- strptime(paste(Date, Time), format="%Y-%m-%d %I:%M %p", tz=tz)
      metric <- TRUE
      tcol <- table[,col]
      Level2 <- gsub("[0-9]+:[0-9]+ [AP]M .+ ([-\\.0-9]+) m", "\\1", tcol)
      if (identical(table[,col], Level2)) {
        Level2 <- gsub("[0-9]+:[0-9]+ [AP]M .+ ([-\\.0-9]+) ft", "\\1", tcol)
        metric <- FALSE
      }
      Tide.Calendar <- rbind(Tide.Calendar, data.frame(Date.Time=Date.Time, Level=ifelse(Level2=="", NA, as.numeric(Level2)), 
                                                       Tide=c("High Tide", "Low Tide", "High Tide", 
                                                              "Low Tide", "High Tide")[col-1], 
                                                       stringsAsFactors=FALSE))
    }
  }
  Tide.Calendar <- na.omit(Tide.Calendar)
  Tide.Calendar <- Tide.Calendar[order(Tide.Calendar[,1]),]
  if (!metric) Tide.Calendar[, "Level"] <- Tide.Calendar[, "Level"] * 0.3048
  Sys.setlocale(category = "LC_TIME", locale = tl )
  return(Tide.Calendar)
}
