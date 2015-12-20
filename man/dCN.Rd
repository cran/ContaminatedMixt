\name{dCN}
\alias{dCN}
\alias{rCN}
\title{Multivariate Contaminated Normal Distribution}
\description{
Probability density function and random number generation for the multivariate contaminated normal distribution.
}
\usage{
dCN(x, mu = rep(0,p), Sigma, alpha = 0.99, eta = 1.01)
rCN(n, mu = rep(0,p), Sigma, alpha = 0.99, eta = 1.01)  
}
                 
%- maybe also 'usage' for other objects documented here.
\arguments{                      
  \item{x}{
either a vector of length \code{p} or a matrix with \code{p} columns, being \code{p = ncol(Sigma)}, representing the coordinates of the point(s) where the density must be evaluated
}
  \item{mu}{
either a vector of length \code{p}, representing the mean value, or (except for \code{rCN}) a matrix whose rows represent different mean vectors; if it is a matrix, its dimensions
must match those of \code{x}
}
  \item{Sigma}{
a symmetric positive-definite matrix representing the scale matrix of the distribution; a vector of length 1 is also allowed (in this case, \code{p = 1} is set)
}
  \item{alpha}{
proportion of good observations; it must be a number between 0 and 1
}
  \item{eta}{
degree of contamination; it should be a number greater then 1
}
  \item{n}{
the number of random vectors to be generated
}
}
\value{
\code{dCN} returns a vector of density values; \code{rCN} returns a matrix of \code{n} rows of random vectors
}

\references{
Punzo, A. and McNicholas, P. D. (2013). Robust Clustering via Parsimonious Mixtures of Contaminated Gaussian Distributions. \emph{arXiv.org} e-print \bold{1305.4669}, available at: \url{http://arxiv.org/abs/1305.4669}.
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
