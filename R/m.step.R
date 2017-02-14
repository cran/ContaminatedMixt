m.step <- function(X,modelname, z, mtol=1e-10, mmax=10){
  n = nrow(X)
  p = ncol(X)
  G = ncol(z)
  sigmar    = matrix(0, nrow=G, ncol=p^2 )
  invsigmar = matrix(0, nrow=G, ncol=p^2 )  
  mu = matrix(0, nrow = G, ncol = p)
  PXgood <- matrix(0, nrow=n, ncol=G)

  temp_em<-.C("mstepU", 
  as.integer(n), as.integer(p), as.integer(G), as.double(z), #1
  as.double(sigmar), as.double(invsigmar), as.double(mu), as.double(mtol), #2
  as.integer(mmax), as.double(X), as.character(modelname), as.double(PXgood),#3
  PACKAGE="ContaminatedMixt")
  Sigma    = array(temp_em[[5]], dim= c(p,p,G),dimnames = list(paste("X.",1:p,sep=""),paste("X.",1:p,sep=""),paste("group ",1:G,sep=""))) 
  invSigma = array(temp_em[[6]], dim= c(p,p,G),dimnames = list(paste("X.",1:p,sep=""),paste("X.",1:p,sep=""),paste("group ",1:G,sep="")) )
  mu       = matrix(temp_em[[7]], nrow=p, ncol=G, byrow=TRUE,dimnames=list(paste("X.",1:p,sep=""),paste("group ",1:G,sep="")))
  px = matrix(temp_em[[12]], nrow = n, ncol = G)
  return(list(mu=mu,Sigma=Sigma,invSigma=invSigma,px=px))
}