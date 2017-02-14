print.ContaminatedMixt.CV <- function(x,...){
  if (length(x$models) >0) {
    best <- whichBestCV(x,...)
    if (length(x$models)>1){
        m <- paste(names(best)[1],collapse=", ")
        m <- paste("\nBest model according to", x$models[[best]]$k, "fold cross validation is")
      } else m <- "\nCross validation for the"
      m <- paste0(m, ifelse(x$models[[best]]$contamination,
                            " contaminated"," uncontaminated"),", with")
      m <- paste(m, "G =", x$models[[best]]$G,"group(s),")
      if (!is.null(x$models[[best]]$model)) 
        m <- paste(m, "and parsimonious structure",  x$models[[best[[1]]]]$model)
      
      cat(m,"\n")
  }
  else cat("No models have been estimated.")
  invisible(x)
}
