#' plot.PriorsmcmcComposite plot a prior
#' @title Plot a prior defined with setPriors function
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return A ggplot object
#' @param x The priors to show
#' @param parameter The name or rank of prior to show
#' @param ... Not used
#' @description Create a ggplot graph with prior.\cr
#' The function makes minimal effort to decorate the plot.
#' @family mcmcComposite functions
#' @importFrom rlang .data
#' @examples
#' \dontrun{
#' library("HelpersMG")
#' par <- c(a0=10, a1=2, b2=20, b1=-1)
#' rules <- rbind(data.frame(Name="^a", Min=0, Max="x*2"), 
#'               data.frame(Name="^b", Min=0, Max=100))
#' p <- setPriors(par=par, se=NULL, density="dgamma", rules=rules)
#' plot(p, parameter="a0")
#' q <- plot(p, parameter="b1")
#' q + geom_line(color = "red") + theme_bw() + 
#' theme(plot.margin=unit(c(2,1,1,1), 'cm'), 
#'       panel.border = element_blank(), 
#'       axis.line.x.bottom = element_line(colour = "black"), 
#'       axis.line.y.left = element_line(colour = "black")) + 
#' labs(title="Parameter: b1") + theme(plot.title = element_text(hjust = 0.5))
#' }
#' @method plot PriorsmcmcComposite
#' @export


plot.PriorsmcmcComposite <- function(x                                 , 
                                     parameter=1                      , 
                                     ...                                ) {
  
    x_axis <- seq(from=x[parameter, "Min"]-(x[parameter, "Max"]-x[parameter, "Min"])*0.01, 
                  to=x[parameter, "Max"]+(x[parameter, "Max"]-x[parameter, "Min"])*0.01, 
                  length.out=1000)
    y_axis <- do.call(what=x[parameter, "Density"], 
                 args=list(x=x_axis, x[parameter, "Prior1"], x[parameter, "Prior2"], log=FALSE))
    df <- data.frame(x=x_axis, Density=y_axis)
    pgg <- ggplot(data = df, 
                  aes(x = .data[["x"]], 
                      y = .data[["Density"]])) + geom_line()
  
  return(pgg)
}
