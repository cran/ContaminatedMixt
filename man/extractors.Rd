\name{Extractor functions}
\alias{getBestModel}
\alias{getPosterior}
\alias{getSize}
\alias{getCluster}
\alias{getPar}
\alias{getIC}
\alias{getDetection}
\alias{whichBest}
\alias{summary.ContaminatedMixt}
\alias{print.ContaminatedMixt}

%\alias{print.ContaminatedMixt.IC}



%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extractors for \code{ContaminatedMixt} Class Objects.
}
\description{
These functions extract values from \code{ContaminatedMixt} class objects. 
}
\usage{

getBestModel(object, criterion = "BIC", G = NULL, model = NULL)
getPosterior(object, ...)
getSize(object, ...)
getCluster(object, ...)
getPar(object, ...)
getIC(object,criteria)
getDetection(object,...)
whichBest(object, criteria = NULL, G = NULL, model = NULL)

\method{summary}{ContaminatedMixt}(object, criterion = "BIC", 
  digits = getOption("digits")-2, \dots)
\method{print}{ContaminatedMixt}(x, \dots)

}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object, x}{a class \code{ContaminatedMixt} object.
  }
  \item{criterion}{a string with the information criterion to consider; supported values are: \code{"AIC"}, \code{"AICc"}, \code{"AICu"}, \code{"AIC3"}, \code{"AWE"}, \code{"BIC"}, \code{"CAIC"}, \code{"ICL"}. 
		Default value is \code{"BIC"}. 
}  
  \item{criteria}{a vector of strings  with the names of information criteria to consider. 
  If \code{NULL}, all the supported information criteria are considered.
}
  \item{G}{an optional vector containing the numbers of groups to consider. 
  If not specified, all the estimated models are considered.
  }
\item{model}{an optional vector of character strings indicating the parsimonious models to consider. 
If not specified, all the estimated models are considered. 
  }
  \item{digits}{ 
        integer used for number formatting.
}
  \item{\dots}{   
  additional arguments to be passed to \code{getBestModel} (or to \code{whichBest} for the \code{print} method). 
}
  
}
\details{
When several models have been estimated, these functions consider the best model according to the information criterion in \code{criterion}, among the estimated models having a number of components among those in \code{G} and a parsimonious model among those in \code{model}.
\code{getIC} provides values for the information criteria in \code{criteria}.

The \code{getBestModel} method returns a \code{ContaminatedMixt} object containing the best model only, selected as described above.
}

