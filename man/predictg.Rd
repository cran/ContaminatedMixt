\name{CNpredict}
\alias{CNpredict}
\alias{predict.ContaminatedMixt}
\title{Cluster Prediction}
\description{
Cluster prediction for multivariate observations based on uncontaminated/contaminated normal mixture models 
}
\usage{
CNpredict(newdata, prior, mu, invSigma, eta=NULL, alpha=NULL)
\method{predict}{ContaminatedMixt}(object, newdata, ...)
}
                 
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{newdata}{
a \code{dim=c(n,p)} matrix representing the coordinates of \code{n} new data point(s)
}
  \item{object}{an object of class \code{ContaminatedMixt} resulting from a call to \code{\link{CNmixt}}. 
	When several models have been estimated, \code{\link{getBestModel}} is used to  select one of them}
   \item{\dots}{
Options to be passed to \code{\link{getBestModel}}
}
 \item{prior}{
a vector with \code{length=G}, where \code{G} is the number of components of the mixture model. 
Its \eqn{k}th component is the mixing proportion for the \eqn{k}th component
}
  \item{mu}{
a \code{dim=c(p,G)} matrix with mean values for each component of the mixture model
}
  \item{invSigma}{
an array with \code{dim=c(p,p,G)} whose element \code{invSigma[,,k]} is the inverse covariance matrix for the \eqn{k}th component of the mixture model.
}
  \item{alpha}{
a vector of \code{length=G} with the proportions of good observations; it must be a number between 0 and 1. 
Use \code{NULL} for uncontaminated models
}
  \item{eta}{
a vector of \code{length=G} with the degree of contamination; it should be a number greater than 1. 
Use \code{NULL} for uncontaminated models
}
  
}
\value{
a vector with group membership
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
\code{\link{ContaminatedMixt-package}}}

\examples{

point <- c(0,0,0)
mu <- c(1,-2,3)
Sigma <- diag(3)
alpha <- 0.8
eta <- 5
f <- dCN(point, mu, Sigma, alpha, eta)
x <- rCN(10, mu, Sigma, alpha, eta)

}
