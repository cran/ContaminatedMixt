# ------------------------- #
# posteriors initialization #
# ------------------------- #
.postInit <- function(initialization,X,n,p,G,start.z,label,modelname, threshold,start){
if(initialization=="random.post"){
  z  <- array(runif(n*G),c(n,G)) # soft posterior probabilities (no-normalized) (n x G) 
  z  <- z/rowSums(z)             # soft posterior probabilities (n x G)

} 

if(initialization=="random.clas"){
  z <- t(rmultinom(n, size = 1, prob=rep(1/G,G)))  # hard posterior probabilities (n x G)
} 

if(initialization=="manual"){ 
      z  <- start.z 
}

if(initialization=="kmeans"){
  clusters  <- kmeans(x=X, centers=G, nstart = 5)      
  z         <- mclust::unmap(clusters$cluster,1:G)
} 

if(initialization=="mixt"){
  mixturefit <- mixture::gpcm(data=X, G=G, mnames=modelname, start=start, 
                              label=label, veo=TRUE, atol=threshold, pprogress=FALSE)
  z <- mixturefit$z
  
  if(is.null(z) & p>1){
    z <- mclust::Mclust(data=X, G = G, modelNames= modelname)$z
  }
}

if(is.null(z)){
  msg <- paste0("\nModel ",modelname," with G = ",G," was not estimated due to bad initialization.")
  stop(msg)
}

if (min(colSums(z))<=p) z <- (z+0.0000001)/rowSums(z+0.0000001)

# z for labeled observations
if(length(unique(label))>1) z[label!=0] <- mclust::unmap(label[label!=0], G=G)
z
}
.checkPar <- function (parameter,G,value){
  if(is.null(parameter)){
    parameter <- rep(value,G) 
  } else{
    if(length(parameter) == 1) parameter <- rep(parameter,G)
    if(length(parameter)!=G) {
      parameter <- rep(parameter[1],G)
    }
  }
  parameter
}