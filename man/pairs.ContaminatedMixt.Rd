\name{pairs.ContaminatedMixt}
\alias{pairs.ContaminatedMixt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scatterplot Matrix for ContaminatedMixt Objects
}
\description{
A matrix of scatterplots, for objects of class \code{\link{ContaminatedMixt}}, is produced.
}
\usage{
\method{pairs}{ContaminatedMixt}(x, criterion = "BIC", \ldots) 
}

\arguments{
  \item{x}{
an object of class \code{\link{ContaminatedMixt}}
}
  \item{criterion}{an optional character string with the information criterion to consider; supported values are: \code{"AIC"}, \code{"AICc"}, \code{"AICu"}, \code{"AIC3"}, \code{"AWE"}, \code{"BIC"}, \code{"CAIC"}, \code{"ICL"}. 
  Default value is \code{"BIC"}.  
}
  \item{\dots}{
Options to be passed to \code{pairs}.
}
}


\author{
Antonio Punzo, Angelo Mazza, Paul D. McNicholas 
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ContaminatedMixt}}
}
