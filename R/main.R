CNmixt <- function(
  X,                            # matrix of data
  G,                            # vector with number of groups to be evaluated
  contamination= NULL, #c("C","U"),
  model=NULL,                   # models to be considered in model selection
  initialization="mixt",        # initialization procedure: "random.post", "random.clas", "manual", or "mixt"
  alphafix=NULL,                   # vector of dimension G with proportion of good observations in each group
  alphamin=0.5,                 # vector of minimum proportions of good data 
  seed=NULL,
  start.z=NULL,                 # (n x k)-matrix of soft or hard classification: it is used only if initialization="manual"    
  start.v=NULL,                 # (n x 2 x k)-array of soft or hard classification in each group: it is used only if initialization="manual"  	
  start=0,                      # initialization for the package mixture
  label=NULL,                   # groups of the labelled observations
  AICcond = FALSE,
  iter.max=1000,                # maximum number of iterations in the EM-algorithm
  threshold=1.0e-03,            # stopping rule in the Aitken rule
  parallel = FALSE,
  eps=1.0e-100
){
  args=mget(names(formals()),sys.frame(sys.nframe()))
  args$doCV=FALSE
  res = do.call("CNmixt_main", args)
  res$call= match.call()
  res
}
CNmixtCV <- function(
  X,                            # matrix of data
  G,                            # vector with number of groups to be evaluated
  contamination= NULL, #c("C","U"),
  model=NULL,                   # models to be considered in model selection
  initialization="mixt",        # initialization procedure: "random.post", "random.clas", "manual", or "mixt"
  k = 10,
  alphafix=NULL,                   # vector of dimension G with proportion of good observations in each group
  alphamin=0.5,                 # vector of minimum proportions of good data 
  seed=NULL,
  start.z=NULL,                 # (n x k)-matrix of soft or hard classification: it is used only if initialization="manual"    
  start.v=NULL,                 # (n x 2 x k)-array of soft or hard classification in each group: it is used only if initialization="manual"  	
  start=0,                      # initialization for the package mixture
  label=NULL,                   # groups of the labelled observations
  iter.max=1000,                # maximum number of iterations in the EM-algorithm
  threshold=1.0e-03,            # stopping rule in the Aitken rule
  parallel = FALSE,
  eps=1.0e-100
){
  args=mget(names(formals()),sys.frame(sys.nframe()))
  args$doCV=TRUE
  do.call("CNmixt_main", args)
  res$call= match.call()
  res
}