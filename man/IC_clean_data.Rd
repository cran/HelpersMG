% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IC_clean_data.R
\name{IC_clean_data}
\alias{IC_clean_data}
\title{Clean the dataframe before to be used with IC_threshold_matrix}
\usage{
IC_clean_data(data = stop("A dataframe object is required"),
  use = c("pairwise.complete.obs", "everything", "all.obs",
  "complete.obs", "na.or.complete"), method = c("pearson", "kendall",
  "spearman"), variable.retain = NULL, test.partial.correlation = TRUE,
  progress = TRUE, debug = FALSE)
}
\arguments{
\item{data}{The data.frame to be cleaned}

\item{use}{an optional character string giving a method for computing covariances in the presence of missing values. This must be (an abbreviation of) one of the strings "everything", "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs".}

\item{method}{a character string indicating which correlation coefficient (or covariance) is to be computed. One of "pearson" (default), "kendall", or "spearman": can be abbreviated.}

\item{variable.retain}{a vector with the name of columns to keep}

\item{test.partial.correlation}{should the partial correlations be tested ?}

\item{progress}{Show a progress bar}

\item{debug}{if TRUE, information about progression of cleaning are shown}
}
\value{
A dataframe
}
\description{
This function must be used if missing values are present in the dataset.\cr
It ensures that all correlations and partial correlations can be calculated. 
The columns of the dataframe are removed one per one until all can be calculated without error. 
It is possible to say that one or more columns must be retained because they are of particular importance in the analysis. 
The use and method parameters are used by cor() function. The function uses by default a parallel computing in Unix or MacOSX systems. 
If progress is TRUE and the package pbmcapply is present, a progress bar is displayed. If debug is TRUE, some informations are shown during the process.
\code{https://fr.wikipedia.org/wiki/Iconographie_des_corrélations}
}
\details{
IC_clean_data checks and corrects the dataframe to be used with IC_threshold_matrix
}
\examples{
\dontrun{
library("HelpersMG")
es <- matrix(c("e1", "52", "12", "12", "5",
"e2", "59", "12.5", "9", "5",
"e3", "55", "13", "15", "9",
"e4", "58", "14.5", "5", "5",
"e5", "66", "15.5", "11", "13.5",
"e6", "62", "16", "15", "18",
"e7", "63", "17", "12", "18",
"e8", "69", "18", "9", "18"), ncol=5, byrow = TRUE)
colnames(es) <- c("Élève", "Poids", "Âge", "Assiduité", "Note")
es <- as.data.frame(es, stringsasFactor=FALSE)
es[, 2] <- as.numeric(as.character(es[, 2]))
es[, 3] <- as.numeric(as.character(es[, 3]))
es[, 4] <- as.numeric(as.character(es[, 4]))
es[, 5] <- as.numeric(as.character(es[, 5]))

es

df <- IC_clean_data(es, debug = TRUE)
cor_matrix <- IC_threshold_matrix(data=df, threshold = NULL, progress=FALSE)
cor_threshold <- IC_threshold_matrix(data=df, threshold = 0.3)
par(mar=c(1,1,1,1))
set.seed(4)
plot(cor_threshold)
cor_threshold_Note <- IC_correlation_simplify(matrix=cor_threshold, variable="Note")
plot(cor_threshold_Note)
}
}
\references{
Lesty, M., 1999. Une nouvelle approche dans le choix des régresseurs de la régression multiple en présence d’interactions et de colinéarités. Revue de Modulad 22, 41-77.
}
\seealso{
Other Iconography of correlations: \code{\link{IC_correlation_simplify}},
  \code{\link{IC_threshold_matrix}},
  \code{\link{plot.IconoCorel}}
}
\author{
Marc Girondot
}
\concept{Iconography of correlations}