CNpredict <- function(newdata,prior,mu,invSigma,eta=NULL,alpha=NULL){
  n <- nrow(newdata)    # sample size
  p <- ncol(newdata)    # number of variables
  G <- length(prior)
  if (!all(dim(mu)==c(p,G), dim(invSigma)==c(p,p,G),length(prior)==c(G))) stop("Error in the dimensions of arguments.")
  group <-rep(0,n)
  if(any(sapply(list(newdata,mu,invSigma,prior),function(x) missing(x)))) stop("Required argument missing")
  if(is.null(eta)) temp_em <- .C("RestepU",
              as.integer(group), as.integer(n), as.integer(p), as.integer(G), as.double(newdata), #1
              as.double(t(mu)),as.double(invSigma), as.double(prior),
              PACKAGE="ContaminatedMixt")
  else{
    if (!all(length(eta)==c(p,G), length(alpha)==c(p,p,G))) stop("Error in the dimensions of arguments.")
    temp_em<-.C("RestepC",
              as.integer(group),as.integer(n), as.integer(p), as.integer(G), as.double(newdata), #1
              as.double(t(mu)),as.double(invSigma),as.double(eta),as.double(alpha),
              as.double(prior),
              PACKAGE="ContaminatedMixt")
  }
  return(temp_em[[1]])
}

predict.ContaminatedMixt <- function(object,newdata,...){
  model = getBestModel(object,...)
  with(model$models[[1]], CNpredict(newdata, prior, mu,invSigma,eta,alpha))
}

  