#' dnbinom_new returns density for the negative binomial distribution
#' @title Random numbers for the negative binomial distribution. 
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return Random numbers for the negative binomial distribution
#' @param x vector of (non-negative integer) quantiles.
#' @param size target for number of successful trials, or dispersion parameter (the shape parameter of the gamma mixing distribution). Must be strictly positive, need not be integer.
#' @param prob probability of success in each trial. 0 < prob <= 1.
#' @param mu alternative parametrization via mean.
#' @param sd alternative parametrization via standard deviation.
#' @param var alternative parametrization via variance.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @description Density for the negative binomial distribution with parameters mu, sd, var, size or prob. See \code{dnbinom}.
#' @examples
#' \dontrun{
#' library("HelpersMG")
#' set.seed(1)
#' x <- rnbinom_new(n=100, mu=2, sd=3)
#' LnL <- NULL
#' df <- data.frame(mu=seq(from=0.1, to=8, by=0.1), "-LnL"=NA)
#' for (mu in df[, "mu"])
#' LnL <- c(LnL, -sum(dnbinom_new(x=x, mu=mu, sd=3, log=TRUE)))
#' df[, "-LnL"] <- LnL
#' ggplot(data = df, aes(x = .data[["mu"]], y = .data[["-LnL"]])) + geom_line()
#' # Examples of wrong parametrization
#' dnbinom_new(x=x, mu=c(1, 2), sd=3, log=TRUE)
#' }
#' @export

dnbinom_new <- function(x, size=NULL, prob=NULL, mu=NULL, sd=NULL, var=NULL, log=FALSE) {
  # size=NULL; prob=NULL; mu=NULL; sd=NULL; var=NULL; log=FALSE
  
  if (!is.null(size)) if (length(size) > 1) {
    warning("size can be only of size 1. Only the first element is used.")
    size <- size[1]
  }
  if (!is.null(prob))  if (length(prob) > 1) {
    warning("prob can be only of size 1. Only the first element is used.")
    prob <- prob[1]
  }
  if (!is.null(mu))  if (length(mu) > 1) {
    warning("mu can be only of size 1. Only the first element is used.")
    mu <- mu[1]
  }
  if (!is.null(sd))  if (length(sd) > 1) {
    warning("sd can be only of size 1. Only the first element is used.")
    sd <- sd[1]
  }
  if (!is.null(var))  if (length(var) > 1) {
    warning("var can be only of size 1. Only the first element is used.")
    var <- var[1]
  }
  # Les combinatoires possibles sont
  # size=NULL, prob=NULL # OK
  # size=NULL, mu=NULL # OK
  # size=NULL, sd=NULL # ok
  # size=NULL, var=NULL ¨# ok
  # prob=NULL, mu=NULL # OK
  # prob=NULL, sd=NULL # ok
  # prob=NULL, var=NULL # ok
  # mu=NULL, sd=NULL # OK
  # mu=NULL, var=NULL # OK
  # sd=NULL, var=NULL # Impossible
  
  if (!is.null(sd) & !is.null(var)) if (var != sd^2) {
    warning("sd and var are incompatible.")
    return(rep(NA, length(x)))
  }
  # # prob = size/(size+mu)
  # if (is.null(size) & !is.null(prob) & !is.null(mu)) size <- (prob * mu) / (1- prob)
  if (!is.null(sd)) var <- sd^2
  # prob = size/(size+mu)
  if (is.null(mu) & !is.null(prob) & !is.null(size)) mu <- (size/prob)-size
  # # prob = size/(size+mu)
  # if (is.null(prob) & !is.null(size) & !is.null(mu)) prob <- size/(size+mu)
   # var = mu + mu^2 / size
  if (is.null(size) & !is.null(prob) & !is.null(var)) {
    # size <- mu^2 / (var - mu)
    # mu <- ((size/prob)-size)
    # size <- (prob * mu) / (1- prob)
    # size <- mu^2 / (var - mu)
    # size * (var - mu) <- mu^2
    # size * var - size - mu - mu^2 <- 0
    
    # size <- mu^2 / (var - mu)
    # size <- ((size/prob)-size)^2 / (var - ((size/prob)-size))
    # size <- ((size/prob)-(prob*size)/prob)^2 / (var - (size/prob)+size)
    # size <- ((size/prob)^2-2*((prob*size)/prob)*(size/prob)+((prob*size)/prob)^2) / ((var*prob - size + size*prob)/prob)
    # size <- ((size^2/prob^2)-2*((prob*size^2)/prob^2)+((prob^2*size^2)/prob^2)) / ((var*prob - size + size*prob)/prob)
    # size <- ((size^2/prob^2)-2*((prob*size^2)/prob^2)+((prob^2*size^2)/prob^2)) * (prob/(var*prob - size + size*prob))
    # size <- (((size^2)-2*((prob*size^2))+((prob^2*size^2)))/prob^2) * (prob/(var*prob - size + size*prob))
    # size <- (((size^2)-2*((prob*size^2))+((prob^2*size^2)))/prob) * (1/(var*prob - size + size*prob))
    # size <- (((size^2-2*(prob*size^2)+(prob^2*size^2))/prob)/1) * (1/(var*prob - size + size*prob))
    
    # size <- (((size^2-2*(prob*size^2)+(prob^2*size^2))/prob)/(var*prob - size + size*prob))
    
    # size - (((size^2-2*(prob*size^2)+(prob^2*size^2))/prob)/(var*prob - size + size*prob)) <- 0
    # (size*(var*prob - size + size*prob))/(var*prob - size + size*prob) - (((size^2-2*(prob*size^2)+(prob^2*size^2))/prob)/(var*prob - size + size*prob)) <- 0
    # ((size*var*prob - size^2 + size^2*prob) - ((size^2-2*(prob*size^2)+(prob^2*size^2))/prob))/(var*prob - size + size*prob) <- 0
    # ((size*var*prob - size^2 + size^2*prob) - ((size^2-2*(prob*size^2)+(prob^2*size^2))/prob)) <- 0
    
    #  ((size*var*prob^2 - size^2*prob + size^2*prob^2)/prob - ((size^2-2*(prob*size^2)+(prob^2*size^2))/prob)) <- 0
    #  ((size*var*prob^2 - size^2*prob + size^2*prob^2) - ((size^2-2*(prob*size^2)+(prob^2*size^2))))/prob <- 0
    #  size*var*prob^2 - size^2*prob + size^2*prob^2 - size^2+2*(prob*size^2)-prob^2*size^2 <- 0
    #  var*prob^2 - size*prob + size*prob^2 - size+2*prob*size-prob^2*size <- 0
    #  var*prob^2  <- size*prob - size*prob^2 + size -2*prob*size +prob^2*size
    
    #  var*prob^2  <- size*(prob - prob^2 + 1 -2*prob +prob^2)
    size <- (var*prob^2) / (1-prob)
    mu <- (size/prob)-size
  }
  
  if (is.null(size) & !is.null(mu) & !is.null(var)) size <- mu^2 / (var - mu)
  
  if (is.null(mu) & !is.null(size) & !is.null(var)) {
    # size <- mu^2 / (var - mu)
    # size * var -size * mu <- mu^2
    # mu^2 + size * mu - size * var <- 0
    A <- 1
    B <- size
    C <- - size * var
    delta <- B^2-4*A*C
    if (delta < 0) {
      warning("mu cannot be negative nor mu > var")
      return(rep(NA, length(x)))
    }
    m1 <- (-B+sqrt(delta))/(2*A)
    m2 <- (-B-sqrt(delta))/(2*A)
    if (m1 > 0) {
      mu <- m1
    } else {
      mu <- m2
    }
  }
  
  if (is.null(size) | is.null(mu)) {
    warning("Not enough information were provided.")
    return(rep(NA, length(x)))
  }

  if (size < 0) {
    warning("size cannot be negative nor mu > var")
    return(rep(NA, length(x)))
  } else {
    if ((mu == 0) & any(x != 0)) {
      warning("When mu is 0, x cannot be different from 0!")
    }
    return(dnbinom(x=x, size=size, mu=mu, log=log))
  }
}
