plot.ContaminatedMixt <- function(x, criterion="BIC",contours=FALSE, xmarg=1, ymarg=2, res=200, levels=seq(.0001,1,by=0.01), ...){
  criterion <- match.arg(criterion,.ICnames(x$models[[1]]))
  bivres <- getBestModel(x,criterion=criterion)$models[[1]]
  plot(bivres$X[,c(xmarg,ymarg)], col="white", ...)
  if(contours){
    lims <- par()$usr
    xseq <- seq(lims[1], lims[2], length.out=res)
    yseq <- seq(lims[3], lims[4], length.out=res)
    resgood <- resbad <- rescont <- array(0,c(res,res,bivres$G))
    val     <- array(0,c(res,res))
    for(g in 1:bivres$G){
      m=bivres$mu[c(xmarg,ymarg),g]
      S=bivres$Sigma[c(xmarg,ymarg),c(xmarg,ymarg),g]
      resgood[,,g] <- outer(xseq,yseq,.bnorm,m=m,S=S)
      if (bivres$contamination){
        resbad[,,g]  <- outer(xseq,yseq,.bnorm,m=m,bivres$eta[g]*S)
        rescont[,,g] <- bivres$alpha[g]*resgood[,,g]+(1-bivres$alpha[g])*resbad[,,g]
        val          <- val + bivres$prior[g]*rescont[,,g]
        }
      else val <- val + bivres$prior[g]*resgood[,,g]
    }
    contour(x=xseq, y=yseq, z=val, add=TRUE, levels=levels, col=rgb(0.5,0.5,0.5,alpha=0.7))
  }
  labels = if (bivres$contamination) bivres$detection$innergroup else "*"
  text(bivres$X[,c(xmarg,ymarg)], labels=labels, col=bivres$group, ...)
  
}


# #####################################
# ## Bivariate Gaussian Distribution ##
# #####################################
# 
# .bivnorm <- function(x,y,mu,Sigma){
#   muX    <- mu[1]
#   muY    <- mu[2]
#   sigmaX <- sqrt(Sigma[1,1])
#   sigmaY <- sqrt(Sigma[2,2])
#   rho    <- Sigma[1,2]/(sigmaX*sigmaY)
#   mah    <- 1/(1-rho^2)*( (x-muX)^2/sigmaX^2 - 2*rho*(x-muX)*(y-muY)/(sigmaX*sigmaY) + (y-muY)^2/sigmaY^2 )
# 
#   1/(2*pi*sigmaX*sigmaY*sqrt(1-rho^2))*exp(-1/2*mah)
# 
# }

.bnorm <- function(X,Y,m,S){
  mvtnorm::dmvnorm(matrix(c(X,Y),ncol = 2),m,S)
}