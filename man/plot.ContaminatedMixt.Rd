\name{plot.ContaminatedMixt}
\alias{plot.ContaminatedMixt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scatterplot for ContaminatedMixt Objects
}
\description{
Scatterplot, with optionally superimposed contours, for objects of class \code{\link{ContaminatedMixt}}.
}
\usage{
\method{plot}{ContaminatedMixt}(x, criterion = "BIC", contours = FALSE, xmarg = 1, ymarg = 2,  
                     res = 200, levels = seq(.0001,1,by=0.01), \dots) 
}

\arguments{
  \item{x}{
an object of class \code{\link{ContaminatedMixt}}
}
  \item{criterion}{a string with the information criterion to consider; supported values are: \code{"AIC"}, \code{"AICc"}, \code{"AICu"}, \code{"AIC3"}, \code{"AWE"}, \code{"BIC"}, \code{"CAIC"}, \code{"ICL"}. 
		Default value is \code{"BIC"}. 
} 
  \item{contours}{
if \code{TRUE}, the contours of the mixture density are superimposed on the plot.
Default is \code{FALSE}.
}
  \item{xmarg}{
scalar argument giving the position of the variable to be used on the \eqn{x}-axis.
}
  \item{ymarg}{
scalar argument giving the position of the variable to be used on the \eqn{y}-axis.
}
  \item{res}{
scalar argument giving the resolution for the calculation grid required for the contour plot. 
Default is 200, which results in a \eqn{200 \times 200} grid.
}
  \item{levels}{
Numeric vector giving the levels at which contours should be drawn. 
Default is to draw a contour in 0.01 steps, starting from the contour of height .0001. 
This may result in more/less contours than desired depending on the resulting density.
}
  \item{\dots}{
Options to be passed to \code{plot}.
}
}

\author{
Antonio Punzo, Angelo Mazza, Paul D. McNicholas 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{CNmixt}}
}
