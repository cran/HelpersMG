#' growlnotify send growl notification for MacOs X systems.
#' @title Send growl notification for MacOs X system.
#' @author Marc Girondot
#' @return None
#' @param textinfo Text to display in the growlnotify window
#' @param help If TRUE, an help is displayed
#' @description This function is used to send a notification to MacOS user.
#' @examples 
#' # If growlnotify is used on a non-mac system, it just quits.
#' growlnotify("It works if you are on a Mac with GrowlNotify installed!")
#' @export


growlnotify <-
function(textinfo="", help=FALSE) {
if(help) {
	cat("This function is used to send a notification to Mac user.\n")
	cat("The syntax is growlnotify(textinfo='Done!')\n")

} else {
	if (.Platform$OS.type=="unix") {
		if (substr(.Platform$pkgType, 1, 3)=="mac") {
			if (any(system("ls /usr/local/bin/", intern=TRUE)=="growlnotify")) {
				system(paste("/usr/local/bin/growlnotify -t 'R script info' -m '", textinfo, "' -a 'R'", sep=""))
			}
		}
	}
}
}
