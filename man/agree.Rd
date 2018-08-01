\name{agree}
\alias{agree}
\title{Agreement Between Partitions}
\description{
Evaluates the agreement of a given partition with respect to the partition arising from the mixture of multivariate contaminated normal distributions.
If the mixture has been fitted for classification purposes, the agreement will be based on the unlabeled observations only.
}
\usage{
agree(object, givgroup, criterion = "BIC")  
}
                 
%- maybe also 'usage' for other objects documented here.
\arguments{                      
  \item{object}{
An object of class \code{\link{ContaminatedMixt}}
}
  \item{givgroup}{
vector, of the same dimension of the number of observations used to fit the model in \code{object}, representing a given partition
}
  \item{criterion}{an optional character string with the information criterion to consider; supported values are: \code{"AIC"}, \code{"AICc"}, \code{"AICu"}, \code{"AIC3"}, \code{"AWE"}, \code{"BIC"}, \code{"CAIC"}, \code{"ICL"}. 
  Default value is \code{"BIC"}.  
}
}
\value{
A contingency table.
}

\references{
Punzo A., Mazza A. and McNicholas P. D. (2018). \pkg{ContaminatedMixt}: An R Package for Fitting Parsimonious Mixtures of Multivariate Contaminated Normal Distributions. \emph{Journal of Statistical Software}, \bold{85}(10), 1--25.

Punzo A. and McNicholas P. D. (2016). Parsimonious mixtures of multivariate contaminated normal distributions. \emph{Biometrical Journal}, \bold{58}(6), 1506--1537.
}
\author{
Antonio Punzo, Angelo Mazza, Paul D. McNicholas 
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{                  
\code{\link{ContaminatedMixt-package}}, 
\code{\link{CNmixt}}
}