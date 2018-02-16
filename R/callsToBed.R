## A method to convert viterbi or marginal paths to a 0-based [,) formated bed file containing state and species labels

#' @export
callsToBed <- function(hmm,marginal=NULL,viterbi=NULL,marginalThresh=0.5){
    if(is.null(marginal) & is.null(viterbi)){
        stop("Either marginal or viterbi path must be passed to function.")
    }
    if(!is.null(marginal) & !is.null(viterbi)){
        stop("Cannot pass both viterbi and marginal to function")
    }
    ## Get which states are the most likely under the specified criteria
    if(!is.null(marginal)){
        states=lapply(marginal,function(x) rle(unlist(apply(x>=marginalThresh,1,which))))
    } else {
        states=lapply(viterbi,rle)
    }
    ## Get chromosomal loci of all bins that are not the full background state
    loci=list()
    for(i in 1:length(states)){
        end=cumsum(states[[i]]$lengths)[which(states[[i]]$values>1)] ## 0-based end bin indicies [,)
        start=end-states[[i]]$lengths[which(states[[i]]$values>1)]
        if(length(start>0)){
            loci[[as.character(i)]]=data.table::data.table(chrom=hmm$emission$invariants$bed[i,]$chrom,
                    start=hmm$emission$invariants$bed[i,]$start + hmm$emission$invariants$binSize * start,
                    end=hmm$emission$invariants$bed[i,]$start + end*hmm$emission$invariants$binSize,
                    state=states[[i]]$values[states[[i]]$values>1])            
        }
    }
    loci=data.table::rbindlist(loci)
    ## Add species that have active elements in each state
    loci[,species:=paste0(hmm$emission$invariants$tree$tip.label[as.logical(hmm$emission$invariants$treeStates[state[1],])],collapse=";"),by="state"]
    loci=loci[order(chrom,start,end)]
    setcolorder(loci, c("chrom", "start", "end","species","state"))
    return(loci)
}
