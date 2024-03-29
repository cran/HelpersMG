#' LD50_MHmcmc_p generates set of parameters to be used with LD50_MHmcmc()
#' @title Generates set of parameters to be used with LD50_MHmcmc()
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return A matrix with the parameters
#' @param result An object obtained after a LD50 fit
#' @param accept If TRUE, the script does not wait user information
#' @description Interactive script used to generate set of parameters to be 
#' used with LD50_MHmcmc().
#' @family LD50 functions
#' @examples 
#' \dontrun{
#' library("HelpersMG")
#' data <- data.frame(Doses=c(80, 120, 150, 150, 180, 200),
#' Alive=c(10, 12, 8, 6, 2, 1),
#' Dead=c(0, 1, 5, 6, 9, 15))
#' LD50_logistic <- LD50(data, equation="logistic")
#' pmcmc <- LD50_MHmcmc_p(LD50_logistic, accept=TRUE)
#' }
#' @export


LD50_MHmcmc_p<-function(result=stop("An output from LD50() must be provided"), 
                       accept=FALSE) {

# d'abord je sors les paramtres  utiliser

par <- result$par

# "P"

P <- c("dnorm", par["P"], par["P"]/10, 2, par["P"]/2, par["P"]*2, ifelse(is.na(par["P"]), 29.5, par["P"]))
S <- c("dnorm", par["S"], abs(par["S"])/2, 0.5, -abs(par["S"])*2, abs(par["S"])*2, 
       ifelse(is.na(par["S"]), 0.01, par["S"]))
K <- c("dnorm", par["K"], abs(par["K"])/2, 0.5, -abs(par["K"])*2, abs(par["K"])*2, 
       ifelse(is.na(par["K"]), 0.01, par["K"]))
K1 <- c("dnorm", par["K"], abs(par["K1"])/2, 0.5, -abs(par["K1"])*2, abs(par["K1"])*2, 
        ifelse(is.na(par["K1"]), 0.01, par["K1"]))
K2 <- c("dnorm", par["K2"], abs(par["K2"])/2, 0.5, -abs(par["K2"])*2, abs(par["K2"])*2, 
        ifelse(is.na(par["K2"]), 0.01, par["K2"]))


priors <- list(P=P, S=S, K=K, K1=K1, K2=K2)

prencours <- NULL

for (i in 1:length(par)) {
	prencours <- c(prencours, priors[[names(par)[i]]])
}

parametersMCMC <- matrix(prencours, ncol=7, byrow=T)
colnames(parametersMCMC) <- c("Density", "Prior1", "Prior2", "SDProp", "Min", "Max", "Init")
rownames(parametersMCMC)<-names(par)

parametersMCMC <- as.data.frame(parametersMCMC, stringsAsFactors = FALSE)

for (i in 2:7)
  parametersMCMC[,i] <- as.numeric(parametersMCMC[,i])

parameters <- parametersMCMC

if (accept) {

	return(parameters)

} else {

	repeat {

cat("Proposition:\n")
print(parameters)
cat("Name of the parameter to change or Enter to quit:\n")
f<-scan(nmax=1, quiet=TRUE, what=character())

if (length(f)==0) f <- "q"

if (f=="q") {
	return(parameters)
	
} else {

	variable <- which(f==names(par))
	if (length(variable)==0) {
	cat("The parameter does not exist:\n")
	} else {
	print(variable)
	cat(paste("Change for the parameter ",names(par)[variable],":\n",sep=""))

	cat(paste("Distribution of the prior (Enter for default ",parameters[variable, "Density"], "):", sep=""))
	density<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(density)!=0) { parameters[variable, "Density"] <- density } else { density <- parameters[variable, "Density"] }
	
	if (density == "dunif") {
	
	cat(paste("Distribution of the prior, Minimum (Enter for default ", parameters[variable, "Prior1"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Prior1"] <- f
	cat(paste("Distribution of the prior, Maximum (Enter for default ", parameters[variable, "Prior2"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Prior2"] <- f
	
	} else {
	
	if (density == "dnorm") {
	
	cat(paste("Distribution of the prior, Mean (Enter for default ", parameters[variable, "Prior1"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Prior1"] <- f
	cat(paste("Distribution of the prior, Standard deviation (Enter for default ",parameters[variable, "Prior2"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Prior2"] <- f
	
	} else {

	cat(paste("Distribution of the prior, value 1 (Enter for default ",parameters[variable, "Prior1"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Prior1"] <- f
	cat(paste("Distribution of the prior, value 2 (Enter for default ",parameters[variable, "Prior2"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Prior2"] <- f

	}
	}
	
	
	cat(paste("SD of new proposition (Enter for default ",parameters[variable, "SDProp"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "SDProp"] <- f
	cat(paste("Minimum for the parameter (default ",parameters[variable, "Min"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Min"] <- f
	cat(paste("Maximum for the parameter (Enter for default ",parameters[variable, "Max"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Max"] <- f
	cat(paste("Initial value (Enter for default ",parameters[variable, "Init"], "):", sep=""))
	f<-scan(nmax=1, quiet=TRUE, what=character())
	if (length(f)!=0) parameters[variable, "Init"] <- f
	}

}
}

}

}
