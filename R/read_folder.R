#' read.folder reads all files present in a folder
#' @title Reads all files present in a folder and creates a list with the content of these files
#' @author Marc Girondot
#' @return Return a list of the data in the files of the folder (directory for windows users)
#' @param folder Where to search for files; can be or a file path or a folder path
#' @param wildcard Define which files are to be read (examples: "*.*", "*.xls", "essai*.txt")
#' @param read Function used to read file. Ex: read.delim or read.xls from gdata package
#' @param ... Parameters send to the read function
#' @description To create a list, the syntax is:\cr
#' datalist<-read.folder(folder=".", read=read.delim, header=FALSE)\cr
#' It returns NULL with a warning if the folder does not exist or is empty.\cr
#' The names of the elements of the list are the filenames.\cr
#' @examples 
#' \dontrun{
#' library(HelpersMG)
#' # Read all the files from a folder/directory
#' Gratiot<-read.folder(folder=".", wildcard="*.csv", read=read.csv2)
#' }
#' @export


read_folder <- function(folder=try(file.choose(), silent=TRUE), 
                        wildcard="*.*", read=read.delim, ...) {

	if (class(folder)!="try-error") {
    fi <- file.info(folder)
    if (is.na(fi$isdir)) {
      warning("Folder/directory does not exist")
      return(invisible())
    }
    if (!fi$isdir) {
      folder <- dirname(folder)
    }
    lf <- Sys.glob(file.path(folder, wildcard))
  	ladd <- list()
	
    	if (length(lf)!=0) {
	      previous<-getwd()
	      setwd(folder)
	    	for (i in 1:length(lf)) {
      		linter <- read(lf[i], ...)
		      ladd <- c(ladd, list(linter))	
      	}
      	names(ladd) <- lf
      	setwd(previous)
      	return(ladd)
	    } else {
	      warning("No selected files in folder/directory")
		    return(invisible())
	    }
	}

}
