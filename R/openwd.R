#' openwd will open a finder window with current working directory
#' @title Open a finder window with current working directory in MacOS X and windows
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return A vector with the list of files.
#' @param directory The directory you want to open
#' @description This function opens a finder window with directory files 
#' in MacOS X. It has not been fully tested in Windows. In linux, it just returns the 
#' list of files in directory.\cr
#' By defaut, it uses the current working directory.
#' @examples 
#' \dontrun{
#' openwd()
#' }
#' @export


openwd <- function(directory=getwd()) {
  if (.Platform$OS.type == "unix") {
    if (Sys.info()["sysname"] == "Darwin") {
      system(command=sprintf('open %s', shQuote(directory)))
    }
  } else {
    if (.Platform$OS.type == "windows"){
      shell.exec(directory)
    }
  }
  return(list.files(directory))
}
