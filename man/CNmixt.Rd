\name{CNmixt}
\alias{CNmixt}
\title{Fitting for the Parsimonious Mixtures of Contaminated Normal Distributions}
\description{
Fits, by using the expectation conditional-maximization (ECM) algorithm, parsimonious mixtures of multivariate contaminated normal distributions (with eigen-decomposed scale matrices) to the given data within a clustering paradigm (default) or classification paradigm. 
Can be run in parallel.
Likelihood-based model selection criteria are used to select the parsimonious model and the number of groups.
}
\usage{
CNmixt(X, G, model = NULL, initialization = "mixt",
   alphafix = NULL, alphamin = 0.5, etafix = NULL, etamax = 1000,
   seed = NULL, start.z = NULL, start.v = NULL, start = 0,
   ind.label = NULL, label = NULL, iter.max = 1000, threshold = 1.0e-03, 
   parallel = FALSE, eps = 1e-100)  
}
                 
%- maybe also 'usage' for other objects documented here.
\arguments{                      
  \item{X}{
  a matrix or data frame such that \eqn{n} rows correspond to observations and \eqn{p} columns correspond to variables. 
  Note that this function currently only works with multivariate data (\eqn{p > 1}).
  }
  \item{G}{
  a vector containing the numbers of groups to be tried. 
  }
  \item{model}{
  a vector indicating the model(s) to be fitted.
  Possible values are: \code{"EII"}, \code{"VII"}, \code{"EEI"}, \code{"VEI"}, \code{"EVI"}, \code{"VVI"}, \code{"EEE"}, \code{"VEE"}, \code{"EVE"}, \code{"EEV"}, \code{"VVE"}, \code{"VEV"}, \code{"EVV"}, \code{"VVV"}.
	If \code{NULL}, then all 14 models are fitted.
  }
  \item{initialization}{
  initialization strategy for the ECM algorithm. 
  It can be:
  \itemize{
    \item \code{"mixt"} (default): the initial (\eqn{n \times G}) soft classification matrix (of posterior probabilities of groups membership) arises from a preliminary run of mixtures of multivariate normal distributions as fitted by the \code{gpcm()} function of the \pkg{mixture} package (see \code{\link[mixture:gpcm]{mixture:gpcm}} for details).
    \item \code{"kmeans"}: the initial (\eqn{n \times G}) hard classification matrix arises from a preliminary run of the \eqn{k}-means algorithm;
    \item \code{"random.soft"}: the initial (\eqn{n \times G}) soft classification matrix (of posterior probabilities of groups membership) is randomly generated; 
    \item \code{"random.hard"}: the initial (\eqn{n \times G}) hard classification classification matrix is randomly generated;
    \item \code{"manual"}: the user must specify the initial (\eqn{n \times G}) soft/hard classification matrix, via the argument \code{start.z} and, optionally, the initial (\eqn{n \times G}) matrix of posterior probabilities to be a good observation in each group, via the argument \code{start.v}.    
  	}
  }
  \item{alphafix}{
  a vector of length \eqn{G} with the proportion of good observations in each group. 
  If \code{length(alphafix) != G}, then the first element is replicated \eqn{G} times. 
  Default value is \code{NULL}.
  }
  \item{alphamin}{
  a vector of length \eqn{G} with the minimum proportion of good observations in each group.  
  If \code{length(alphamin) != G}, then the first element is replicated \eqn{G} times.
  Default value is \code{0.5}. 
  }
  \item{etafix}{a vector of length \eqn{G} with the
  values of the contamination parameter to be considered in the estimation phase for each group. 
  If \code{length(etafix) != G}, then the first element is replicated \eqn{G} times. 
	Default value is \code{NULL}.
  }

  \item{etamax}{a vector of length \eqn{G} with the
  maximum value for the contamination parameter to be considered in the estimation phase for each group when \code{etafix} is \code{NULL}. 
	If \code{length(etamax) != G}, then the first element is replicated \eqn{G} times.
	Default value is \code{1000}.
  }
  \item{seed}{
  the seed for the random number generator, when random initializations are used; if \code{NULL}, current seed is not changed. 
  Default value is \code{NULL}.
  }
  \item{start.z}{
  initial \eqn{n \times G} matrix of either soft or hard classification. 
  Default value is \code{NULL}.
  }
  \item{start.v}{
  initial \eqn{n \times G} matrix of posterior probabilities to be a good observation in each group. 
  Default value is a \eqn{n \times G} matrix of ones.
  }
  \item{start}{
  when \code{initialization = "mixt"}, initialization used for the \code{gpcm()} function of the \pkg{mixture} package (see \code{\link[mixture:gpcm]{mixture:gpcm}} for details). 
  }
  \item{ind.label}{
  vector of positions (rows) of the labeled observations. 
  }
  \item{label}{
  vector, of the same dimension as \code{ind.label}, with the group of membership of the observations indicated by \code{ind.label}. 
  }
  \item{iter.max}{
  maximum number of iterations in the ECM algorithm. 
	Default value is \code{1000}. 
  }
  \item{threshold}{
  threshold for Aitken's acceleration procedure. 
	Default value is \code{1.0e-03}. 
  }
  \item{parallel}{
  When \code{TRUE}, the package \code{\link[parallel:parallel-package]{parallel}} is used for parallel computation. 
  When several models are estimated, computational time is reduced. 
  The number of cores to use may be set with the global option \code{cl.cores}; default value is detected using \code{\link[parallel:detectCores]{detectCores()}}.
}  
  \item{eps}{
  an optional scalar. 
  It sets the smallest value for the eigenvalues of the component scale matrices. 
  Default value is \code{1e-100}.
}

}
\details{
The multivariate data contained in \code{X} are either clustered or classified using parsimonious mixtures of multivariate contaminated normal distributions with some or all of the 14 parsimonious models described in Punzo and McNicholas (2015).
Model specification (via the \code{model} argument) follows the nomenclature popularized in other packages such as \pkg{mixture} and \pkg{mclust}.  
Such a nomenclature refers to the decomposition and constraints on the scale matrix (see Banfield and Raftery, 1993, Celeux and Govaert, 1995 and Punzo and McNicholas, 2015 for details): 
\deqn{\Sigma_g = \lambda_g \Gamma_g \Delta_g \Gamma_g'.}
The nomenclature describes (in order) the volume (\eqn{\lambda_g}), shape (\eqn{\Delta_g}), and orientation (\eqn{\Gamma_g}), in terms of \code{"V"}ariable, \code{"E"}qual, or the \code{"I"}dentity matrix. 
As an example, the string \code{"VEI"} would refer to the model where \eqn{\Sigma_g = \lambda_g \Delta}.
Note that for \eqn{G=1}, several models are equivalent (for example, \code{"EEE"} and \code{"VVV"}). 
Thus, for \eqn{G=1} only one model from each set of equivalent models will be run.

The algorithms detailed in Celeux and Govaert (1995) are considered in the first CM-step of the ECM algorithm to update \eqn{\Sigma_g} for all the models apart from \code{"EVE"} and \code{"VVE"}.
For \code{"EVE"} and \code{"VVE"}, majorization-minimization (MM) algorithms (Hunter and Lange, 2000) and accelerated line search algorithms on the Stiefel manifold (Absil, Mahony and Sepulchre, 2009 and Browne and McNicholas, 2014), which are especially preferable in higher dimensions (Browne and McNicholas, 2014), are used to update \eqn{\Sigma_g}; the same approach is also adopted in the \pkg{mixture} package for those models.

Starting values are very important to the successful operation of these algorithms and so care must be taken in the interpretation of results.
All the initializations considered here provide initial quantities for the first CM-step of the ECM algorithm.
}
\value{
An object of class \code{ContaminatedMixt} is a list with components:
	\item{call}{an object of class \code{call}}
	\item{best}{a data frame with the best number of mixture components (first column) and the best model (second column) with respect to the two model selection criteria adopted (BIC and ICL)}
	\item{bestBIC,bestICL}{for the best BIC and ICL models, these are two lists (of the same type) with components:
		\itemize{
		\item \code{modelname}: the name of the best model.
		\item \code{npar}: number of free parameters.
		\item \code{X}: matrix of data.
		\item \code{G}: number of mixture components.
		\item \code{p}: number of variables.
		\item \code{prior}: weights for the mixture components.
		\item \code{priorgood}: weights for the good observations in each of the \code{k} groups.
		\item \code{mu}: component means.
		\item \code{Sigma}: component covariance matrices for the good observations.
		\item \code{eta}: component contamination parameters.	
		\item \code{iter.stop}: final iteration of the ECM algorithm.
		\item \code{z}: matrix with posterior probabilities for the outer groups.
		\item \code{v}: matrix with posterior probabilities for the inner groups.
		\item \code{ind.label}: vector of positions (rows) of the labeled observations.
		\item \code{label}: vector, of the same dimension as \code{ind.label}, with the group of membership of the observations indicated by \code{ind.label}.	
		\item \code{group}: vector of integers indicating the maximum a posteriori classifications for the best model.
		\item \code{loglik}: log-likelihood value of the best model.	
		\item \code{BIC}: BIC value	
		\item \code{ICL:} ICL value	
		\item \code{call}: an object of class \code{call} for the best model.
		}		
	}
}

\references{
Absil, P. A., Mahony, R. and Sepulchre, R. (2009). Optimization Algorithms on Matrix Manifolds. Princeton University Press, Princeton, NJ.

Banfield, J. D. and Raftery A. E. (1993). Model-Based Gaussian and Non-Gaussian Clustering. \emph{Biometrics}, \bold{49}(3), 803--821.

Browne, R. P. and McNicholas, P. D. (2013). Estimating Common Principal Components in High Dimensions. \emph{Advances in Data Analysis and Classification}, \bold{8}(2), 217--226.

Browne, R. P. and McNicholas, P. D. (2014). Orthogonal Stiefel manifold optimization for eigen-decomposed covariance parameter estimation in mixture models. \emph{Statistics and Computing},
\bold{24}(2), 203--210.

Browne, R. P. and McNicholas, P. D. (2015). \pkg{mixture}: Mixture Models for Clustering and Classification. R package version 1.4.

Celeux, G. and Govaert, G. (1995). Gaussian Parsimonious Clustering Models. \emph{Pattern Recognition}. \bold{28}(5), 781--793.

Hunter, D. R. and Lange, K. (2000). Rejoinder to Discussion of ``Optimization Transfer Using Surrogate Objective Functions''. \emph{Journal of Computational and Graphical Statistics}, \bold{9}(1), 52--59.

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
\code{\link{ContaminatedMixt-package}}
}
\examples{

## Note that the example is extremely simplified 
## in order to reduce computation time

# Artificial data from an EEI Gaussian mixture with G = 2 components

library("mnormt")
p <- 2
set.seed(12345)
X1 <- rmnorm(n = 200, mean = rep(2, p), varcov = diag(c(5, 0.5)))
X2 <- rmnorm(n = 200, mean = rep(-2, p), varcov = diag(c(5, 0.5)))
noise <- matrix(runif(n = 40, min = -20, max = 20), nrow = 20, ncol = 2)
X <- rbind(X1, X2, noise)

group <- rep(c(1, 2, 3), times = c(200, 200, 20))
plot(X, col = group, pch = c(3, 4, 16)[group], asp = 1, xlab = expression(X[1]),
ylab = expression(X[2]))

# ---------------------- #
# Model-based clustering #
# ---------------------- #

res1 <- CNmixt(X, model = c("EEI", "VVV"), G = 2, parallel = FALSE)

summary(res1)

agree(res1, givgroup = group)

plot(res1, contours = TRUE, asp = 1, xlab = expression(X[1]), ylab = expression(X[2]))

# -------------------------- #
# Model-based classification #
# -------------------------- #

indlab <- sample(1:400, 20)
lab <- group[indlab]
res2 <- CNmixt(X, G = 2, model = "EEI", ind.label = indlab, label = lab)

agree(res2, givgroup = group)

}

