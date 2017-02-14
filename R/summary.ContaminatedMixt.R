summary.ContaminatedMixt <-function(object,criterion="BIC", digits = getOption("digits")-2, ...)
{
  criterion <- match.arg(criterion,.ICnames(object$models[[1]]))
  best <- getBestModel(object,criterion=criterion,...)
  obj <- best$models[[1]]
  title1 <- paste0("Best fitted model according to ",criterion)
  
  nch <- nchar(title1)
  cat(rep("-",nch ),"\n",sep="")
  cat(title1, "\n")
  cat(rep("-", nch),"\n",sep="")
  #
  tab <- data.frame("log-likelihood" = obj$loglik, "n" = length(obj$group), 
                    "par" = obj$npar,row.names = "")
  tab[[criterion]] <- obj$IC[[criterion]]
  
  print(tab, digits = digits)
  #
  cat("\nClustering table:")
  print(table(obj$group), digits = digits)
  #
  cat("\nPrior: ")
  cat(paste(names(obj$prior), format(obj$prior,digits=digits), sep = " = ", 
            collapse = ", "), sep = "")
  cat("\n")
  #
  
  if(!is.null(obj$model)){
    cat("Model: ", ifelse(obj$contamination,"Contaminated ","Uncontaminated "), as.character(obj$model), " (", .ModelNames(obj$model)$type, 
        ") with ", obj$G, ifelse(obj$G > 1, " components\n", " component\n"),
        sep="")
    cat("\n")
  }
  #
  
  if(!is.null(obj$mu)){
    cat("Variables")
    cat("\n Means:\n")
    print(obj$mu, digits = digits)
    cat("\n Variance-covariance matrices:\n") 
    for(i in seq_len(obj$G)){
      cat(paste0("  Component ",i,"\n"))
      print(obj$Sigma[,,i], digits = digits) 
    }
    if(obj$contamination){ 
      cat("\n Alpha\n")
      print(best$models[[1]]$alpha, digits = digits)
      cat("\n Eta\n")
      print(best$models[[1]]$eta, digits = digits)
    }
    
    cat("\n")
  }
}