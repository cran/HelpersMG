#' plot.IconoCorel checks and corrects the dataframe to be used with IC_threshold_matrix
#' @title Clean the dataframe before to be used with IC_threshold_matrix
#' @author Marc Girondot
#' @return A igraph object
#' @param x The correlation matrix to show
#' @param show.legend.direction the position of the legend of direction; FALSE to not show it
#' @param show.legend.strength the position of the legend with intensity of correlation; FALSE to not show it
#' @param vertex.label.color a vector with the colors of labels
#' @param vertex.label a vector with the labels
#' @param vertex.color a vector of colors
#' @param vertex.label.cex a vector of cex
#' @param title the title of the plot
#' @param plot if TRUE, the plot is shown
#' @param ... other options of plot.igraph()
#' @description This function plots the data as a network. It returns an invisible object that can be used with visIgraph from package visNetwork.
#' \code{https://fr.wikipedia.org/wiki/Iconographie_des_corrélations}
#' @family Iconography of correlations
#' @references Lesty, M., 1999. Une nouvelle approche dans le choix des régresseurs de la régression multiple en présence d’interactions et de colinéarités. Revue de Modulad 22, 41-77.
#' @examples
#' \dontrun{
#' library("HelpersMG")
#' es <- structure(list(Student = c("e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8"), 
#'                  Mass = c(52, 59, 55, 58, 66, 62, 63, 69), 
#'                  Age = c(12, 12.5, 13, 14.5, 15.5, 16, 17, 18), 
#'                  Assiduity = c(12, 9, 15, 5, 11, 15, 12, 9), 
#'                  Note = c(5, 5, 9, 5, 13.5, 18, 18, 18)), 
#'                  row.names = c(NA, -8L), class = "data.frame")
#' 
#' es
#' 
#' df <- IC_clean_data(es, debug = TRUE)
#' cor_matrix <- IC_threshold_matrix(data=df, threshold = NULL, progress=FALSE)
#' cor_threshold <- IC_threshold_matrix(data=df, threshold = 0.3)
#' par(mar=c(1,1,1,1))
#' set.seed(4)
#' library("igraph")
#' library("visNetwork")
#' kk <- plot(cor_threshold, vertex.color="red")
#' # it can be shown also with the visNetwork package
#' visIgraph(kk)
#' cor_threshold_Note <- IC_correlation_simplify(matrix=cor_threshold, variable="Note")
#' plot(cor_threshold_Note)
#' 
#' # You can record the position of elements and use them later
#' ly <- layout_nicely(kk)
#' plot(cor_threshold, vertex.color="red", layout=ly)
#' }
#' @method plot IconoCorel
#' @export


plot.IconoCorel <- function(x                                          , 
                            ...                                        , 
                            show.legend.direction="bottomright"        , 
                            show.legend.strength="topleft"             , 
                            title="Correlation iconography"            , 
                            vertex.label.color="black"                 , 
                            vertex.label=NULL                          , 
                            vertex.color="white"                       ,
                            vertex.label.cex=1                         ,
                            plot=TRUE                                  ) {
  
  if ((!is.element("igraph", installed.packages()[,1]))) {
    warning("Packages igraph is absent; the network cannot be shown.")
    plot <- FALSE
  }
  
  
  if (identical(show.legend.direction, TRUE)) show.legend.direction="bottomright"
  if (identical(show.legend.strength, TRUE)) show.legend.strength="topleft"
  
  threshold <- x$threshold
  matrix <- x$thresholded_correlation
  if (is.null(matrix)) {
    matrix <- x$correlation
    threshold <- "No threshold"
  }
  
  if (is.null(nrow(matrix)) | (nrow(matrix)==0)) {
    stop("No plot can be shown")
  }
  
  
  graph <- getFromNamespace("graph.adjacency", ns="igraph")(abs(matrix), weighted=TRUE, mode="lower")
  E <- getFromNamespace("E", ns="igraph")
  V <- getFromNamespace("V", ns="igraph")
  
  if (is.null(vertex.label)) {
    vertex.label <- colnames(matrix)
  }
  
  w <- getFromNamespace("edge_attr", ns="igraph")(graph, "weight", E(graph))
  if (length(w) <= 1) {
    ww <- 1
  } else {
    ww <- 1 + 4*(w-min(w))/(max(w)-min(w))
  }
  # E(g)[idx]$attr <- value is equivalent to g <- set_edge_attr(g, attr, E(g)[idx], value).
    graph <- getFromNamespace("set_edge_attr", ns="igraph")(graph, "width", E(graph), ww)
  graph <- getFromNamespace("set_edge_attr", ns="igraph")(graph, "color", E(graph), sapply(strsplit(attributes(E(graph))$vnames, "\\|"), 
                                                                                           FUN=function(k) ifelse(matrix[k[1], k[2]]>0, "blue", "red")))
  graph <- getFromNamespace("set_edge_attr", ns="igraph")(graph, "lty", E(graph), sapply(strsplit(attributes(E(graph))$vnames, "\\|"), 
                                                                                         FUN=function(k) ifelse(matrix[k[1], k[2]]>0, 1, 2)))
  
  # The assignment form of $ is a shortcut for set_vertex_attr, e.g. V(g)[idx]$attr <- value is equivalent to g <- set_vertex_attr(g, attr, V(g)[idx], value).
  
  graph <- getFromNamespace("set_vertex_attr", ns="igraph")(graph, "label", V(graph), vertex.label)
  graph <- getFromNamespace("set_vertex_attr", ns="igraph")(graph, "color", V(graph), vertex.color)
  graph <- getFromNamespace("set_vertex_attr", ns="igraph")(graph, "label.color", V(graph), vertex.label.color)
  graph <- getFromNamespace("set_vertex_attr", ns="igraph")(graph, "label.cex", V(graph), vertex.label.cex)
  
  
  
  # par(mar=c(2, 2, 4, 2))
  # set.seed(4)
  if (plot) {
    
    getFromNamespace("plot.igraph", ns="igraph")(graph, 
                                                 edge.width = ww, 
                                                 edge.color = sapply(strsplit(attributes(E(graph))$vnames, "\\|"), 
                                                                     FUN=function(k) ifelse(matrix[k[1], k[2]]>0, "blue", "red")), 
                                                 vertex.label = vertex.label, 
                                                 vertex.color = vertex.color, 
                                                 vertex.size = 1,
                                                 edge.lty = sapply(strsplit(attributes(E(graph))$vnames, "\\|"), 
                                                                   FUN=function(k) ifelse(matrix[k[1], k[2]]>0, 1, 2)), 
                                                 vertex.frame.color = NA, 
                                                 vertex.label.color = vertex.label.color, 
                                                 margin=c(0, 0, 0, 0)
                                                 , ...)
    
    title(main = title)
    if (!identical(show.legend.direction, FALSE)) {
      legend(show.legend.direction, legend=c("Positive relationship", "Negative relationship"), 
             col=c("blue", "red"), lty=c(1, 2))
    }
    
    if (!identical(show.legend.strength, FALSE)) {
      if (length(w) <= 1) {
        legend(show.legend.strength, title="Correlation coefficient",
               legend=specify_decimal(w), 
               lty=1, 
               lwd=1)
      } else {
        ml <- round(max(w)+0.1, 1)
        if (ml>1) ml <- 1
        legend(show.legend.strength, title="Correlation coefficient",
               legend=as.character(seq(from=round(threshold, 1), to=ml, by=0.1)), 
               lty=1, 
               lwd=1+4*(seq(from=round(threshold, 1), to=ml, by=0.1)-min(w))/
                 (max(w)-min(w)))
      }
    }
  }
  
  return(invisible(graph))
}
