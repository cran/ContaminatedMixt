\name{wine}
\docType{data}
\alias{wine}
\title{Wine Data Set}
\description{
These data are the results of a chemical analysis of wines grown in the same region in Italy but derived from three different cultivars. 
The analysis determined the quantities of 13 constituents found in each of the three types of wine: Barolo, Grignolino, Barbera. 
The data set is used to evaluate the ability of the \code{CNmixt()} function in clustering the data assuming unknown their cultivars. 
}
\usage{
data(wine)
}
\format{
This data frame contains 178 rows, each corresponding to a different cultivar of wine produced in Piedmont (Italy), and 14 columns. 
The first column is the type of wine (\code{Type}), a factor variable with the following levels: \code{Barolo}, \code{Grignolino}, \code{Barbera}. 
The variables measured on the three types of wines are the following: 
\code{Alcohol}, \code{Malic acid}, \code{Ash}, \code{Alcalinity}, \code{Magnesium}, \code{Phenols}, \code{Flavanoids}, \code{Nonflavanoids}, \code{Proanthocyanins}, \code{Color intensity}, \code{Hue}, \code{OD280.OD315Dilution}, \code{Proline}. 
All variables but the label class are continuous.
}
\details{
The original data set comprises 27 variables. 
Here a subset of 14 variables only has been included.
}
\source{This dataset is from the UCI machine learning repository and it is available at \url{http://archive.ics.uci.edu/ml/datasets/Wine}.}
\references{

Forina M., Lanteri S. Armanino C., Casolino C., Casale M., Oliveri, P. (2008). 
V-PARVUS. \emph{An Extendible Pachage of programs for esplorative data analysis, classification and regression analysis}. 
Dip. Chimica e Tecnologie Farmaceutiche ed Alimentari, Universita' di Genova.}

\seealso{
\code{\link{ContaminatedMixt-package}}, \code{\link{CNmixt}}
}

\keyword{datasets}

\examples{

data("wine")

group <- wine[, 1]
pairs(wine[, -1], cex = 0.6, pch = c(2, 3, 1)[group], col = c(3, 4, 2)[group], gap = 0,
cex.labels = 0.6)

res3 <- CNmixt(wine[, -1], G = 3, model = "EEE", initialization = "random.post",
  seed = 5, parallel = FALSE)
agree(res3, givgroup = group)
pairs(res3, cex = 0.6, gap = 0, cex.labels = 0.6)

}
