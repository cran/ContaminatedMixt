########################################
## Density of a contaminated Gaussian ##
########################################

dCN <- function(x, mu = rep(0,p), Sigma, alpha = 0.99, eta = 1.01){
  
  if(missing(Sigma))
    stop("Sigma is missing")
  if(alpha<0 | alpha>1)
    stop("alpha must be in (0,1)")
  if(eta<1)
    stop("eta must be greater than 1")
  
  p <- if(is.matrix(Sigma)) 
    ncol(Sigma)
  else 1
  if(p == 1) 
    return(alpha*dnorm(x, mu, sqrt(Sigma), log = FALSE)+(1-alpha)*dnorm(x, mu, sqrt(eta*Sigma), log = FALSE))
  x <- if(is.vector(x)) 
    matrix(x, 1, p)
  else data.matrix(x)
  if(is.vector(mu)) 
    mu <- outer(rep(1, nrow(x)), mu)
  if(is.matrix(mu) && (nrow(mu) != nrow(x) || ncol(mu) != ncol(x))) 
    stop("mismatch of dimensions of 'x' and 'mu'")
  if(is.vector(mu)) 
    mu <- outer(rep(1, nrow(x)), mu)
  X <- t(x - mu)
  
  # good
  
  concgood    <- mnormt::pd.solve(Sigma, silent=FALSE, log.det = TRUE)
  Qgood       <- apply((concgood %*% X) * X, 2, sum)
  log.detgood <- attr(concgood, "log.det")
  logPDFgood  <- as.vector(Qgood + p * logb(2 * pi) + log.detgood)/(-2)
  
  # bad
  concbad     <- mnormt::pd.solve(eta*Sigma, silent=FALSE, log.det = TRUE)
  Qbad        <- apply((concbad %*% X) * X, 2, sum)
  log.detbad  <- attr(concbad, "log.det")
  logPDFbad   <- as.vector(Qbad + p * logb(2 * pi) + log.detbad)/(-2)
  
  PDF <- alpha*exp(logPDFgood)+(1-alpha)*exp(logPDFbad)
  PDFmod <- (PDF<=10^(-323))*10^(-323)+(PDF>10^(-323))*PDF
  
  return(PDFmod)
  
}

.dCN <- function(x, mean = rep(0, d), varcov,alpha=0.99,eta=1.01, log = FALSE){
  d <- if (is.matrix(varcov)) 
    ncol(varcov)
  else 1
  if(alpha<0 | alpha>1)
    stop("alpha must be in (0,1)")
  if(eta<1)
    stop("eta must be greater than 1")
  if(d == 1) 
    return(alpha*dnorm(x, mean, sqrt(varcov), log = FALSE)+(1-alpha)*dnorm(x, mean, sqrt(eta*varcov), log = FALSE))
  x <- if(is.vector(x)) 
    matrix(x, 1, d)
  else data.matrix(x)
  n <- nrow(x)
  pdf= matrix(0, nrow=n, ncol=1)
  temp_em<-.C("dCN",
              as.integer(n), as.integer(d), as.integer(1), as.double(x),as.double(mean), #1
              as.double(varcov), as.double(eta), as.double(alpha), as.double(pdf),
              PACKAGE="ContaminatedMixt")
  pdf= matrix(temp_em[[9]], nrow=n, ncol=1)
  return(pdf)
}
