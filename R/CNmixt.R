.CNmixtG <- function(X, G, initialization,  modelname,  contamination,alphafix, alphamin, 
                     start.z, start.v, start, label, iter.max, threshold, 
                     eps,AICcond, doCV, k,verbose){

   if (!doCV){
    mlist <- .CNmixtG2(X, G, initialization,  modelname,  contamination, alphafix, alphamin, start.z, start.v, start, label, iter.max, threshold, eps,verbose) 
    IC <- .ComputeIC(mlist)
    if(AICcond & length(unique(label))>1){
      mlist2 <- .CNmixtG2(X, G, initialization,  modelname,  contamination,alphafix, alphamin, start.z, start.v, start, label=rep(0,nrow(X)), iter.max, threshold, eps,verbose) 
      loglik3<- CNllikelihood(X,mlist$prior,mlist$mu,mlist$invSigma,mlist$eta,mlist$alpha)
      IC$BEC=mlist$loglik-mlist2$loglik
      IC$AICcond=2*mlist$loglik-4*mlist$obslll+2*loglik3
      if(verbose) cat("\n")
      }
    mlist$IC=IC
  }
  else{
    e = el = vector("numeric",k)
    xl = label!=0   
    xl = as.factor(xl)
    folds = caret::createFolds(xl,k=k)
    for(i in 1:k){
      if(verbose) cat(paste("\n",i,"fold"))
      folds_i = folds
      folds_i[[i]] = NULL # remove the ith fold
      folds_i = unlist(folds_i, use.names = FALSE) 
      lfolds_i = label[folds_i]
      temp <- .CNmixtG2(X[folds_i,], G, initialization,  modelname,  contamination,alphafix, alphamin, start.z, start.v, start, label=rep(0,length(folds_i)), iter.max, threshold, eps,verbose) 
      ltemp = temp$group[lfolds_i!=0]
      el[i] = length(ltemp)
      e[i] = mclust::classError(ltemp,lfolds_i[lfolds_i!=0])$errorRate
    }
    mlist = list(
      model  = modelname,
      contamination = contamination,
      G = G,
      e=e,
      el=el,
      mean_error = mean(e),
      k=k)
  }
    
  return(mlist)
}

.CNmixtG2 <- function(X, G, initialization,  modelname,  contamination, alphafix, alphamin, 
                     start.z, start.v, start, label, iter.max, threshold, 
                     eps, verbose){
  n <- nrow(X)    # sample size
  p <- ncol(X)    # number of variables
  if (p == 1 ) modelname = paste0(modelname,"II")
  posterior = .postInit(initialization,X,n,p,G,start.z,label,modelname, threshold,start)

  sigmar    = matrix(0, nrow=G, ncol=p^2 )
  invsigmar = matrix(0, nrow=G, ncol=p^2 )  
  mu = matrix(0, nrow = G, ncol = p)
  prior      <- numeric(G)
  group <- numeric(n)
  iteration = 0
  llvalue =0
  obslll = 0
  mtol=1e-10
  mmax=10

  if(contamination){
    v <- array(0.99,c(n,G),dimnames=list(1:n,paste("group ",1:G,sep="")))
    npar <- (G-1) + p*G + .ncovpar(modelname=modelname, p=p, G=G) + G
    if(is.null(alphafix)) npar <- npar + G
    if(!is.null(start.v)) v <- start.v  
    eta = rep(1.01,G)
    alpha = .checkPar(alphafix,G,0.999)
    alphamin = .checkPar(alphamin,G,0.5)

  temp_em<-.C("loopC",
  as.integer(n), as.integer(p), as.integer(G), as.double(posterior), #1
  as.double(sigmar), as.double(invsigmar), as.double(mu), as.double(mtol), #2
  as.integer(mmax), as.double(X), as.integer(label),as.character(modelname), #3
  as.integer(iter.max),as.double(threshold),as.double(prior), as.integer(iteration),#4
  as.double(llvalue),  as.double(obslll),as.integer(group), as.double(v), #5
  as.double(eta), as.double(alpha), as.double(!is.null(alphafix)), as.double(alphamin),#6
  as.integer(verbose),
  PACKAGE="ContaminatedMixt")
  
    v = array(temp_em[[20]], dim= c(n,G),dimnames=list(1:n,paste("group ",1:G,sep="")) )
    eta  = temp_em[[21]]
    alpha = temp_em[[22]]
  }else{
    temp_em<-.C("loopU", 
  as.integer(n), as.integer(p), as.integer(G), as.double(posterior), #1
  as.double(sigmar), as.double(invsigmar), as.double(mu), as.double(mtol), #2
  as.integer(mmax), as.double(X), as.integer(label),as.character(modelname), #3
  as.integer(iter.max),as.double(threshold),as.double(prior), as.integer(iteration),#4
  as.double(llvalue),as.double(obslll),as.integer(group),
  as.integer(verbose),
  PACKAGE="ContaminatedMixt")
    npar <- (G-1) + p*G + .ncovpar(modelname=modelname, p=p, G=G) 
    alpha = eta =v = NULL
  }
  if (p==1) modelname = substr(modelname,1,1)
  posterior= array(temp_em[[4]], dim= c(n,G),dimnames=list(1:n,paste("group ",1:G,sep="")) )
  Sigma    = array(temp_em[[5]], dim= c(p,p,G),dimnames = list(paste("X.",1:p,sep=""),paste("X.",1:p,sep=""),paste("group ",1:G,sep=""))) 
  invSigma = array(temp_em[[6]], dim= c(p,p,G),dimnames = list(paste("X.",1:p,sep=""),paste("X.",1:p,sep=""),paste("group ",1:G,sep=""))) 
  mu       = matrix(temp_em[[7]], nrow=p, ncol=G, byrow=TRUE,dimnames=list(paste("X.",1:p,sep=""),paste("group ",1:G,sep="")))
  prior  = temp_em[[15]]
  iteration = temp_em[[16]]
  loglik = temp_em[[17]]
  obslll = temp_em[[18]]
  group = temp_em[[19]]
  
  # Classification Matrix #
  #group <- apply(posterior,1,which.max)
  innergroup  <- numeric(n)
  if(contamination){
    for(i in 1:n)
      innergroup[i] <- ifelse(v[i,group[i]]<0.5,"bad","*")  
  }
  detection <- data.frame(group=group,innergroup=innergroup)
  z.const  <- (posterior<.Machine$double.xmin)*.Machine$double.xmin+(posterior>.Machine$double.xmin)*posterior   # vincolo per evitare i NaN nel calcolo di tau*log(tau)

  result <- list(
    model  = modelname,
    contamination = contamination,
    npar      = npar,
    X         = X,            
    G         = G,            
    p         = p,            
    n         = n,            
    prior     = prior,
    alpha     = alpha,
    mu        = mu,
    Sigma     = Sigma,
    invSigma = invSigma,
    eta       = eta,
    iter.stop = iteration,
    posterior = posterior,
    v         = v,
    label     = label,                   
    group     = group,
    detection = detection,
    loglik    = loglik,
    obslll = obslll,
    entropy = -sum(z.const[-label==0,]*log(z.const[-label==0,])),
    call      = match.call()
  )
  
  class(result) <- "CNmixt"
  return(result)
}
.ComputeIC <- function(mlist){
  loglik <- mlist$loglik  
  npar <- mlist$npar
  label <- mlist$label
  posterior <- mlist$posterior
  G <- mlist$G
  n <- nrow(mlist$X)
  # Information Criteria #
  IC <- list()
  IC$AIC   <- 2*loglik - npar*2
  IC$BIC   <- 2*loglik - npar*log(n)
  IC$AIC3  <- 2*loglik - npar*3  
  IC$CAIC  <- 2*loglik - npar*(1+log(n))
  IC$AWE   <- 2*loglik - 2*npar*(3/2+log(n))  
  z.const  <- (posterior<.Machine$double.xmin)*.Machine$double.xmin+(posterior>.Machine$double.xmin)*posterior   # vincolo per evitare i NaN nel calcolo di tau*log(tau)
  hard.z   <- (matrix(rep(apply(posterior,1,max),G),n,G,byrow=F)==posterior)*1
  ECM      <- sum(hard.z[-label==0,]*log(z.const[-label==0,]))
  IC$ICL   <- IC$BIC+2*ECM

  if (n-npar-1>0) {
    IC$AICc  <- IC$AIC - (2*npar*(npar+1))/(n-npar-1)
    IC$AICu  <- IC$AICc - n*log(n/(n-npar-1))
  }
  else {
    IC$AICc  <-IC$AICu <- -Inf
  }
  
  return(IC)
}


CNllikelihood <- function(X,prior,mu,invSigma,eta=NULL,alpha=NULL){
  n <- nrow(X)    # sample size
  p <- ncol(X)    # number of variables
  G <- length(prior)
  llvalue=0
  if (!all(dim(mu)==c(p,G), dim(invSigma)==c(p,p,G),length(prior)==c(G))) stop("Error in the dimensions of arguments.")
 
  if(any(sapply(list(X,mu,invSigma,prior),function(x) missing(x)))) stop("Required argument missing")
  if(is.null(eta)) temp_em <- .C("RllikelihoodU",
                                 as.double(llvalue), as.integer(n), as.integer(p), as.integer(G), as.double(X), #1
                                 as.double(t(mu)),as.double(invSigma), as.double(prior),
                                 PACKAGE="ContaminatedMixt")
  else{
    if (!all(length(eta)==c(p,G), length(alpha)==c(p,p,G))) stop("Error in the dimensions of arguments.")
    temp_em<-.C("RllikelihoodC",
                as.double(llvalue),as.integer(n), as.integer(p), as.integer(G), as.double(X), #1
                as.double(t(mu)),as.double(invSigma),as.double(eta),as.double(alpha),
                as.double(prior),
                PACKAGE="ContaminatedMixt")
  }
  return(temp_em[[1]])
}




