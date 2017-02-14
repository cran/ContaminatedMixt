CNmixt_main <- function(X,G,contamination,model,initialization,AICcond,alphafix,alphamin,
                        seed,start.z,start.v,start,label,iter.max,threshold,
                        parallel,eps,doCV,k){
  initialization <- match.arg(initialization,c("mixt","kmeans","random.post","random.clas","manual"))
  if(is.data.frame(X)) X <- as.matrix(X) 
  n <- nrow(X)
  p <- ncol(X)    # number of variables
  if (is.null(X))     stop('Hey, we need some data, please! X is null')
  if (!is.matrix(X))  stop('X needs to be in matrix form')
  if (!is.numeric(X)) stop('X is required to be numeric')
  if (n == 1)   stop('nrow(X) is equal to 1')
  if (any(is.na(X)))  stop('No NAs allowed.')
  
  if (is.null(G)) stop('G is NULL')
  G <- as.integer(ceiling(G))
  if (any(G < 1)) stop('G is not a positive integer')
  if(nrow(X)<ncol(X)){
    warning("Dimensionality of data exceeds the sample size: it may result in model-fitting failure") 
  }
  if (is.null(contamination)) contamination = c(TRUE,FALSE)
  else if(!is.logical(contamination)) stop("contamination is not a logical")
  if (is.null(label)) label = rep(0,n)
  else{
    if (any(as.integer(label)!=label) || length(label)!=n) stop ("length must be NULL or an integer vector of length equal to nrow(X)")
    if (length(unique(label))-1>min(G)) {
      warning(paste('models with G<unique(label) will not be estimated'))
      G <- G[G>length(unique(label))]
      }
  }

  modelXnormNames <- if(p==1) c("E","V") else c("EII","VII","EEI","VEI","EVI","VVI","EEE","VEE","EVE","EEV","VVE","VEV","EVV","VVV") 
  model <- if (is.null(model)) modelXnormNames else match.arg(model, modelXnormNames, several.ok = TRUE)
  mm <- expand.grid(list(k=G, model=model,contamination = contamination))

  if (!is.null(model) & 1 %in% G){
  eqmod <-  if(p==1) data.frame(name = modelXnormNames,number=c(1,1)) else eqmod <- data.frame(name = modelXnormNames, number=c(1,1,2,2,2,2,3,3,3,3,3,3,3,3))
    mm <- merge(mm,eqmod,by.x="model",by.y="name")
    mm1 <- mm[mm$k==1,]
    mm2 <- mm1[!duplicated(subset(mm1,select= -model, drop=FALSE)),]
    if (nrow(mm1) > nrow(mm2)){
      cat("With G = 1, some models are equivalent, so only one model from each set of equivalent models will be run.\n")
    }
    mm <- rbind(mm2,mm[mm$k!=1,])
  }
  mm <- mm[order(mm$k),,drop=FALSE]

    job <- function(i){
    cat("\nEstimating model")
    if (!is.null(mm$model[i])) cat(paste0(" ",mm$model[i]))
    cat(ifelse(mm$contamination[i]," contaminated",""))
    cat(paste0(" with G = ",mm$k[i],":"))
     .CNmixtG(
      X=X,  		                      
      G=mm$k[i],                            
      initialization=initialization,      
      modelname=as.character(mm$model[i]), 
      contamination=mm$contamination[i],
      alphafix=alphafix,
      alphamin=alphamin, #rep(alphamin,mm$k[i]),
      seed=seed,
      start.z=start.z,                 		
      start.v=start.v,                   	
      start=start,                      
      label=label,                   
      iter.max=iter.max,                
      threshold=threshold,
      eps=eps,
      AICcond=AICcond,
      doCV=doCV,
      k=k
    )  
  }

  if(parallel){
    cores <- getOption("cl.cores", parallel::detectCores())
    cat(paste("\n Using",cores,"cores\n"))
    cl <- parallel::makeCluster(cores)
    #clusterExport(cl,envir=environment())
    par <- parallel::parLapply(cl=cl,1:nrow(mm),function(i) job(i))
    parallel::stopCluster(cl)
  }
  else {
    par <- lapply(1:nrow(mm),function(i) job(i))
  }
  i<- 1
  cat("\n")
  while (!i > length(par)){
    if (! is.null(par[[i]]$error)){
      cat(paste(par[[i]]$error,"\n"))
      par[[i]] <- NULL
    }
    i<- i + 1
  }
  if (!is.null(par)){
    class = if (doCV)  "ContaminatedMixt.CV" else "ContaminatedMixt"
    res <-
      structure(
        list(
          models = par
        ),              
        class = class
      )
    print(res) 
    invisible(res) 
  } else {
  cat("No model was estimated.\n")
  return(NULL)
  }
}

