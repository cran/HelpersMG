#' growlnotify sent growl notification for MacOs X systems. No more supported.
#' @title Send growl notification for MacOs X system.
#' @author Marc Girondot
#' @return None
#' @param textinfo Text to display in the growlnotify window
#' @description This function was used to send a notification to MacOS user.\cr
#' Growlnotify being no longer supported in MacOSX, this function will be removed in future releases.
#' @examples 
#' \dontrun{
#' # If growlnotify is used on a non-mac system, it just quits.
#' growlnotify("It works if you are on a Mac with GrowlNotify installed!")
#' }
#' @export


growlnotify <-
function(textinfo="") {return(textinfo)}
