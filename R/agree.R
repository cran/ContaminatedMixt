agree <- function(object, givgroup, criterion = "BIC"){
  criterion <- match.arg(criterion,.ICnames(object$models[[1]]))
  best <- getBestModel(object,criterion=criterion)
  res <- best$models[[1]]
  ind.lab <- (1:res$n)[res$label != 0]
  ind.unlab <- (1:res$n)[res$label == 0]
  nunlab <- length(ind.unlab)
  
  #if(length(givgroup) != n) 
  #  stop("'givgroup' must have the same length of the sample on which the PMCGD model was fitted")
  
  groups <- numeric(nunlab)
  cont <- 0
  for(i in ind.unlab){
    cont <- cont+1
    if(res$detection[i,2]=="bad")
      groups[cont] <- "bad points"
    else
      groups[cont] <- paste("",res$group[i],sep=" ")     
  }
  dnn=c("givgroup","groups")
  if(length(ind.lab)==0)
    return(table(givgroup, groups,dnn=dnn))
  else
    return(table(givgroup[-ind.lab], groups,dnn=dnn))
  
}
