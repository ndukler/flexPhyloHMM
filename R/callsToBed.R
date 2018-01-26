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
        ## Get sloppy state locations (Where there is a single background state seperating two adjacent non-background states due to poor peak alignment )
        sloppy.states=which(states[[i]]$values==1 & states[[i]]$lengths==1)
        sloppy.states=sloppy.states[sloppy.states>1 & sloppy.states<length(states[[i]]$lengths)] ## ignore bins at edge of chain
        if(length(start>0)){
            loci[[as.character(i)]]=data.table::data.table(chrom=hmm$emission$invariants$bed[i,]$chrom,
                    start=hmm$emission$invariants$bed[i,]$start + hmm$emission$invariants$binSize * start,
                    end=hmm$emission$invariants$bed[i,]$start + end*hmm$emission$invariants$binSize,
                    state=states[[i]]$values[states[[i]]$values>1],
                    temp=which(states[[i]]$values>1),
                    sloppy=FALSE)
            ## If there are any sloppy states, mark them
            if(length(sloppy.states)>0){
                ## Determine number of species on each side of a sloppy state call
                left.calls=apply(hmm$emission$invariants$treeStates[loci[[as.character(i)]][temp %in% (sloppy.states-1)]$state,drop=FALSE],1,sum)
                right.calls=apply(hmm$emission$invariants$treeStates[loci[[as.character(i)]][temp %in% (sloppy.states+1)]$state,drop=FALSE],1,sum)
                ## Get the rle index of the state that covers fewer species
                sloppy.rle.index=sort(c((sloppy.states-1)[left.calls<right.calls],(sloppy.states+1)[left.calls>right.calls]))
                ## Add * to state column of sloppy calls
                loci[[as.character(i)]][temp %in% sloppy.rle.index, sloppy:=TRUE]
            }            
            loci[[as.character(i)]][,temp:=NULL]            
        }        
    }
    loci=data.table::rbindlist(loci)
    ## Add species that have active elements in each state
    loci[,species:=paste0(hmm$emission$invariants$tree$tip.label[as.logical(hmm$emission$invariants$treeStates[state[1],])],collapse=";"),by="state"]
    loci=loci[order(chrom,start,end)]
    loci[sloppy==TRUE,chrom:=paste("#",chrom)]
    loci[,sloppy:=NULL]
    setcolorder(loci, c("chrom", "start", "end","species","state"))
    return(loci)
}
