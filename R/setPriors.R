#' setPriors is a general function to set priors for MHalgoGen()
#' @title Set priors for MHalgoGen()
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return Return a data.frame with priors
#' @param par Named vector with init value of parameters
#' @param se Named vector with standard error of parameters
#' @param density Named vector with density or single value
#' @param rules List of rules to define priors
#' @param silent If TRUE, do not show warning.
#' @description Set priors for MHalgoGen()\cr
#' @family mcmcComposite functions
#' @examples
#' \dontrun{
#' library(HelpersMG)
#' rules <- rbind(data.frame(Name="^a", Min=0, Max="x*2"), 
#'               data.frame(Name="^b", Min=0, Max=100))
#' par <- c(a0=10, a1=2, b2=20)
#' (p <- setPriors(par=par, se=NULL, density="dgamma", rules=rules))
#' (p <- setPriors(par=par, se=NULL, density="dnorm", rules=rules))
#' (p <- setPriors(par=par, se=NULL, density="dunif", rules=rules))
#' par <- c(a0=10, a1=2, b2=20, b1=-1)
#' (p <- setPriors(par=par, se=NULL, density="dgamma", rules=rules))
#' }
#' @export


setPriors <- function(par=stop("A vector with init values is necessary."), 
                      se=NULL, 
                      density="dunif", 
                      rules=NULL, silent=FALSE) {
  
  # par=c(a0=10, b2=20); se=NULL; density="dunif"; rules=rbind(data.frame(Name="^a", Min=0, Max="x*2"), data.frame(Name="^b", Min=0, Max=100))
  
  if (is.null(se)) se <- abs(par/2)
  if (length(density) == 1) {
    density <- rep(density, length(par))
    names(density) <- names(par)
  }
  if ((length(density) != length(par)) | (length(se) != length(par))) 
    stop("Wrong number of elements in par, se or density.")
  
  df_priors <- data.frame(Density=character(), 
                          Prior1=numeric(), Prior2=numeric(), 
                          SDProp=numeric(), 
                          Min=numeric(), Max=numeric(), 
                          Init=numeric())
  names_rules <- rules[, "Name"]
  
  for (indice.par in seq_along(par)) {
    name_par <- names(par)[indice.par]
    # J'ai besoin de min, max pour chaque paramètre
    rank_rule <- sapply(names_rules, FUN=function(names_rules_ind) grepl(names_rules_ind, names(par)[indice.par]))
    if (any(rank_rule)) {
      vpar <- unname(par[indice.par])
      vse <- unname(se[indice.par])
      Min <- rules[which(rank_rule), "Min", drop=TRUE][1]
      Max <- rules[which(rank_rule), "Max", drop=TRUE][1]
      
      x <- vpar
      # if (is.na(suppressWarnings(as.numeric(Min)))) {
      Min <- eval(parse(text=Min), envir= environment())
      # } else {
      #   Min <- as.numeric(Min)
      # }
      # if (is.na(suppressWarnings(as.numeric(Max)))) {
      Max <- eval(parse(text=Max), envir= environment())
      # } else {
      #   Max <- as.numeric(Max)
      # }
      
      Min <- min(Min, vpar)
      Max <- max(Max, vpar)
      
      df_priors_ec <- NULL
      if (density[name_par]=="dnorm") {
        df_priors_ec <- data.frame(Density="dnorm", 
                                   Prior1=vpar, 
                                   Prior2=vse, 
                                   SDProp=abs(vpar), 
                                   Min=Min, Max=Max, 
                                   Init=vpar, 
                                   row.names = name_par)
      }
      if ((density[name_par]=="dgamma") & (Min < 0) & (!silent)) {
       warning(paste0("Parameter ", name_par, ": dgamma cannot be used with negative values; it is changed to dunif.")) 
      }
      if ((density[name_par]=="dgamma") & (Min >= 0)) {
        # Si Min < 0, je passe en dunif
        r <- vpar/(vse^2)
        a <- vpar*r
        df_priors_ec <- data.frame(Density="dgamma", 
                                   Prior1=a, Prior2=r, 
                                   SDProp=abs(vpar), 
                                   Min=Min, Max=Max, 
                                   Init=vpar, row.names = name_par)
      }
      if ((density[name_par]=="dunif") | is.null(df_priors_ec)) {
        df_priors_ec <- data.frame(Density="dunif", 
                                   Prior1=Min, Prior2=Max, 
                                   SDProp=abs(vpar), 
                                   Min=Min, Max=Max,  
                                   Init=vpar, row.names = name_par)
      }
      df_priors <- rbind(df_priors, df_priors_ec)
    }
  }
  
  df_priors <- addS3Class(df_priors, c("PriorsmcmcComposite", "data.frame"))
  return(df_priors)
}
