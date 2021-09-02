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
#' es <- matrix(c("e1", "52", "12", "12", "5",
#' "e2", "59", "12.5", "9", "5",
#' "e3", "55", "13", "15", "9",
#' "e4", "58", "14.5", "5", "5",
#' "e5", "66", "15.5", "11", "13.5",
#' "e6", "62", "16", "15", "18",
#' "e7", "63", "17", "12", "18",
#' "e8", "69", "18", "9", "18"), ncol=5, byrow = TRUE)
#' colnames(es) <- c("Élève", "Poids", "Âge", "Assiduité", "Note")
#' es <- as.data.frame(es, stringsasFactor=FALSE)
#' es[, 2] <- as.numeric(as.character(es[, 2]))
#' es[, 3] <- as.numeric(as.character(es[, 3]))
#' es[, 4] <- as.numeric(as.character(es[, 4]))
#' es[, 5] <- as.numeric(as.character(es[, 5]))
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
  
  graph <- getFromNamespace("graph.adjacency", ns="igraph")(abs(matrix), weighted=TRUE, mode="lower")
  E <- getFromNamespace("E", ns="igraph")
  V <- getFromNamespace("V", ns="igraph")
  
  if (is.null(vertex.label)) {
    vertex.label <- colnames(matrix)
  }
  
  w <- getFromNamespace("edge_attr", ns="igraph")(graph, "weight", E(graph))
  # E(g)[idx]$attr <- value is equivalent to g <- set_edge_attr(g, attr, E(g)[idx], value).
  graph <- getFromNamespace("set_edge_attr", ns="igraph")(graph, "width", E(graph), 1 + 4*(w-min(w))/(max(w)-min(w)))
  graph <- getFromNamespace("set_edge_attr", ns="igraph")(graph, "color", E(graph), sapply(strsplit(attributes(E(graph))$vnames, "\\|"), 
                                                                                           FUN=function(k) ifelse(matrix[k[1], k[2]]>0, "blue", "red")))
  graph <- getFromNamespace("set_edge_attr", ns="igraph")(graph, "lty", E(graph), sapply(strsplit(attributes(E(graph))$vnames, "\\|"), 
                                                                                         FUN=function(k) ifelse(matrix[k[1], k[2]]>0, 1, 2)))
  
  # The assignment form of $ is a shortcut for set_vertex_attr, e.g. V(g)[idx]$attr <- value is equivalent to g <- set_vertex_attr(g, attr, V(g)[idx], value).
  
  graph <- getFromNamespace("set_vertex_attr", ns="igraph")(graph, "label", V(graph), vertex.label)
  graph <- getFromNamespace("set_vertex_attr", ns="igraph")(graph, "color", V(graph), vertex.color)
  graph <- getFromNamespace("set_vertex_attr", ns="igraph")(graph, "label.color", V(graph), vertex.label.color)
  
  
  
  
  # par(mar=c(2, 2, 4, 2))
  # set.seed(4)
  if (plot) {
    
    
    
  getFromNamespace("plot.igraph", ns="igraph")(graph, 
             edge.width = 1 + 4*(w-min(w))/(max(w)-min(w)), 
              edge.color = sapply(strsplit(attributes(E(graph))$vnames, "\\|"), 
                   FUN=function(k) ifelse(matrix[k[1], k[2]]>0, "blue", "red")), 
              vertex.label = vertex.label, 
              vertex.color = vertex.color, 
              vertex.size = 1,
              edge.lty = sapply(strsplit(attributes(E(graph))$vnames, "\\|"), 
                   FUN=function(k) ifelse(matrix[k[1], k[2]]>0, 1, 2)), 
              vertex.frame.color = NA, 
              vertex.label.color = vertex.label.color, 
              ..., 
              margin=c(0, 0, 0, 0))
  
  title(main = title)
  if (!identical(show.legend.direction, FALSE)) {
  legend(show.legend.direction, legend=c("Positive relationship", "Negative relationship"), 
         col=c("blue", "red"), lty=c(1, 2))
  }
  
  if (!identical(show.legend.strength, FALSE)) {
  legend(show.legend.strength, title="Correlation coefficient",
         legend=as.character(seq(from=floor(threshold*10)/10, to=1, by=0.1)), 
         lty=1, 
         lwd=1+4*(seq(from=floor(threshold*10)/10, to=1, by=0.1)-min(w))/
           (max(w)-min(w)))
  }
  }

  return(invisible(graph))
}
