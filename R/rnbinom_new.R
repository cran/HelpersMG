#' rnbinom_new returns random numbers for the negative binomial distribution
#' @title Random numbers for the negative binomial distribution. 
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return Random numbers for the negative binomial distribution
#' @param n number of observations.
#' @param size target for number of successful trials, or dispersion parameter (the shape parameter of the gamma mixing distribution). Must be strictly positive, need not be integer.
#' @param prob probability of success in each trial. 0 < prob <= 1.
#' @param mu alternative parametrization via mean.
#' @param sd alternative parametrization via standard deviation.
#' @param var alternative parametrization via variance.
#' @description See \code{rnbinom}.
#' @examples
#' \dontrun{
#' library("HelpersMG")
#' set.seed(1)
#' x <- rnbinom_new(n=1000, prob=6.25/(5+6.25), size=6.25)
#' mean(x)
#' sd(x)
#' set.seed(1)
#' x <- rnbinom_new(n=1000, mu=5, sd=3)
#' mean(x)
#' sd(x)
#' set.seed(1)
#' x <- rnbinom_new(n=1000, mu=5, var=3^2)
#' mean(x)
#' sd(x)
#' set.seed(1)
#' x <- rnbinom_new(n=1000, mu=5, size=6.25)
#' mean(x)
#' sd(x)
#' set.seed(1)
#' x <- rnbinom_new(n=1000, size=6.25, var=3^2)
#' mean(x)
#' sd(x)
#' set.seed(1)
#' x <- rnbinom_new(n=1000, prob=6.25/(5+6.25), var=3^2)
#' mean(x)
#' sd(x)
#' # Example of wrong parametrization
#' set.seed(1)
#' x <- rnbinom_new(n=1000, sd=3, var=3^2)
#' set.seed(1)
#' x <- rnbinom_new(n=1000, mu=10, var=3^2)
#' }
#' @export

rnbinom_new <- function(n, size=NULL, prob=NULL, mu=NULL, sd=NULL, var=NULL) {
  
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
    return(rep(NA, n))
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
      return(rep(NA, n))
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
    return(rep(NA, n))
  }

  if (size < 0) {
    warning("size cannot be negative nor mu > var")
    return(rep(NA, n))
  } else {
    return(rnbinom(n=n, size=size, mu=mu))
  }
}
