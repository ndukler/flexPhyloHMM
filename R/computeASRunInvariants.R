#' Pre-computes all values needed to compute the probability of the absorbing state
#'
#' @param tree a phylo object
#' @param rate a mutation rate
#' @param base.freq a vector, (p(0),p(1)) of the probabilities of 0/1 at the root node 
#' @param log if TRUE returns log transition matricies
#' @return list of parameter values that are constant for calculating the probability of the absorbing state accross sites
#' @export
computeASRunInvariants <- function(tree,rate,base.freq.zero){
    ## Ensure the tree is sorted correctly
    tree=reorder(tree,"postorder")
    base.freq=c(base.freq.zero,1-base.freq.zero)
    ## Pre-calculate everything needed for an iterative traversal
    tt=createTraversalTable(tree)
    ## Change co-ordinates to 0-based
    tt$traversal=tt$traversal-1
    tt$tips=tt$tips-1
    ## Pre-calculate list of rate matricies for the tree
    Q.list=expMatList(rate,base.freq.zero,tree$edge.length)
    ## Compute the log value of the partition function over all un-enumerated states 
    ## Na.log=log(1-exp(logSumExp(seqProb(tree,sub.states,parms=list(inactive.prob=base.freq[1],rate=rate))$prob)))
    Q.concat=do.call("rbind",Q.list)
    return(list(tt=tt,Q.concat=Q.concat,rate=rate,lbf=log(base.freq)))
}

