#' Heuristically find noisy state calls
#'
#' @param thmm the tree hmm the the masking will be performed on
#' @param noise.calls A list of bins that are noisy and the corresponding state
#' @param noise.threshold The minimum value of sf.sum required to not be masked
#' @return a list of transition matricies, one entry for each specified branch length
#' @export 
maskNoiseCalls <- function(thmm,noise.calls,noise.threshold=16){
    ## Iterate over chains
    for(i in 1:length(noise.calls)){
        if(nrow(noise.calls[[i]])>0){
            state.indices=as.numeric(unlist(apply(noise.calls[[i]][sf.sum<noise.threshold], 1, function(x) x[1]:x[2]))) ## unroll state indicies
            one.hot=thmm$transition$invariants$treeStates[rep.int(noise.calls[[i]]$state,times=with(noise.calls[[i]],end-start+1)),,drop=F] ## Unroll states
            ## For each species ...
            for(s in 1:ncol(one.hot)){
                ## Set rows where there is a noisy active element to 0 for both data and scale factors
                thmm$emission$data[[i]][[s]][state.indices[one.hot[,s]>0],]=0
                thmm$emission$invariants$scaleFactors[[i]][[s]][state.indices[one.hot[,s]>0]]=0
            }  
        }
    }
    thmm$emission$computeLeafProbabilities() ## update leaf probabilities to maintain consistency
}
