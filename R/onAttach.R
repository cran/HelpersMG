.onAttach <- function(libname, pkgname) {
  actual <- utils::packageDescription(pkgname)[["Version"]]
  
  packageStartupMessage(paste("Welcome in package", pkgname, "version", actual))
  
  # conn <- url("https://hebergement.universite-paris-saclay.fr/marcgirondot/CRAN/HelpersMG/version.txt")
  # 
  # version_get <- try(suppressWarnings(
  #   readLines(con=conn)), silent = TRUE
  # )
  # close(con=conn)
  # if (!(is.null(version_get)) & (!inherits(version_get, "try-error"))) {
  #   if (package_version(actual, strict = TRUE) < package_version(version_get, strict = TRUE)) {
  #     packageStartupMessage('An update is available; use:\ninstall.packages("https://hebergement.universite-paris-saclay.fr/marcgirondot/CRAN/HelpersMG.tar.gz", repos=NULL, type="source")')
  #   } else {
  #     packageStartupMessage("No update is available")
  #   }
  # } else {
  #   packageStartupMessage("No internet connection is available to check for presence of update")
  # }
  
}
