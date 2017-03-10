#' r2norm returns random numbers for Gaussian distributions different at left and right
#' @title Random generation for Gaussian distributions different at left and right 
#' @author Marc Girondot
#' @return r2norm returns random numbers
#' @param n number of observations.
#' @param mean vector of means
#' @param sd_low vector of standard deviations below the mean.
#' @param sd_high vector of standard deviations above the mean.
#' @description Random generation for Gaussian distributions different at left and right 
#' @examples
#' \dontrun{
#' n <- r2norm(1000, mean=25, sd_low=2, sd_high=10)
#' 
#' hist(n)
#' }
#' @export

# Random number when the distribution is known by mediane and 0.025 and 0.975 quantiles
# It will be supposed that each half distribution is Gaussian.
# It should be done better but it is not easy

r2norm <- function(n, mean = 0, sd_low = 1, sd_high = 1) {
  rn <- rnorm(n, mean = 0, sd = 1)
  return(ifelse(rn<0, rn*sd_low+mean, rn*sd_high+mean))
}
