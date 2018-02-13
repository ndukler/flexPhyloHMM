#' Computes site probabilities for a binary character model
#'
#' @param tree a phylo object from the ape package
#' @param tree.states a matrix of 0 and 1 values where each column is a species and each row is a site
#' @param inactive.prob The probability that state=0 @ the root of the tree
#' @param rate the mutational rate
#' @param log if TRUE state probabilies are returned in log space
#' @return a vector containing the per-site probability
alignmentProb <- function(tree,tree.states,inactive.prob=0.95,rate=10^-3,log=TRUE){
  if(is.matrix(tree.states) | is.data.frame(tree.states)){
    if(ncol(tree.states)==length(tree$tip.label)){
      tree.states=data.table::as.data.table(tree.states)
    } else if(nrow(tree.states)==length(tree$tip.label)){
      tree.states=data.table::as.data.table(t(tree.states))
    } else {
      stop("tree.states must have the same number of rows or columns as there are leaves on the tree")
    }
  } else if(!is.data.table(tree.states)){
    stop("tree.states must be a matrix, data.frame, or data.table")
  }
  data.table::setnames(tree.states,colnames(tree.states),tree$tip.label)
  ## Create phyDat alignment objects
  pd=phangorn::phyDat(tree.states,type="USER",levels=c(0,1))
  ## Estimate liklihood of each tree
  tree.fit=phangorn::pml(tree,data=pd,rate=rate,bf=c(inactive.prob,1-inactive.prob))
  ## Convert to probability
  if(log){
    return(tree.fit$siteLik)
  } else {
    return(exp(tree.fit$siteLik))
  }    
}