\name{m.step}
\alias{m.step}
\title{M-step of the EM algorithm for Parsimonious Normal Mixtures}
\description{
Carries out the M-step for EM algorithm
}
\usage{
m.step(X, modelname, z, mtol=1e-10, mmax=10)
}
                 
%- maybe also 'usage' for other objects documented here.
\arguments{                      
  \item{X}{
a matrix such that \eqn{n} rows correspond to observations and \eqn{p} columns correspond to variables. 
}
  \item{modelname}{A three letter sequence indicating the covariance structure.
Possible values are: \code{"EII"}, \code{"VII"}, \code{"EEI"}, \code{"VEI"}, \code{"EVI"}, \code{"VVI"}, \code{"EEE"}, \code{"VEE"}, \code{"EVE"}, \code{"EEV"}, \code{"VVE"}, \code{"VEV"}, \code{"EVV"}, \code{"VVV"}.
}
  \item{z}{
A matrix of weights such that \eqn{n} rows correspond to observations and \eqn{G} columns correspond to groups.
}
\item{mtol}{
The convergence criteria for the M-step if an iterative procedure is necessary.
}
\item{mmax}{
The maximum number of iterations for an iterative procedure.}
}
\value{
A list of the model parameters with the mu, Sigma, invSigma and px for each group.
}

\references{
Punzo, A. and McNicholas, P. D. (2016). Parsimonious mixtures of multivariate contaminated normal distributions. \emph{Biometrical Journal}, DOI: 10.1002/bimj.201500144.
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