\name{agree}
\alias{agree}
\title{Agreement Between Partitions}
\description{
Evaluates the agreement of a given partition with respect to the partition arising from the mixture of multivariate contaminated normal distributions.
If the mixture has been fitted for classification purposes, the agreement will be based on the unlabelled observations only.
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
Punzo, A. and McNicholas, P. D. (2015). Parsimonious mixtures of contaminated Gaussian distributions with application to allometric studies. \emph{arXiv.org} e-print \bold{1305.4669}, available at: \url{http://arxiv.org/abs/1305.4669}.
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