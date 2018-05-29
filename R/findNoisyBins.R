#' Heuristically find noisy state calls
#'
#' @param state.calls a list of vectors where each element of the vector is a state call and each list index is a chain 
#' @param scale.factors list of scale factors
#' @param tree.states 0/1 matrix of tree states
#' @param max.sep Max distance allowed between two states for them to be considered potentially noisy
#' @return a list of transition matricies, one entry for each specified branch length
#' @export 
findNoiseStates <-  function(state.calls,scale.factors,tree.states,max.sep=1){
    ## Get the run length encoding for the viterbi path
    state.rle=lapply(state.calls, rle)
    ## Convert the rle to index ranges
    element.index=lapply(state.rle, function(x) {
        end=cumsum(x$lengths)
        return(data.table::data.table(start=end-x$lengths+1,end=end,state=x$values)[x$values>1,])
    })    
   
    ## Get indicies for sets of elements that should be compared against each other
    noise.collect=lapply(element.index, function(x,max.sep) {
        noise.runs=rle((x$start[-1]-x$end[-nrow(x)]-1)<=max.sep)
        ## Get the ranges for the collections of elements that make up each noisy region
        end.collect=cumsum(noise.runs$lengths)+1
        start.collect=end.collect-noise.runs$lengths
        return(data.table::data.table(start=start.collect[noise.runs$values],end=end.collect[noise.runs$values]))
    },max.sep=max.sep)

    ## Subset element.index to contain only elements in noisy collections
    noisy.elements=element.index
    for(i in 1:length(element.index)){
        if(nrow(noise.collect[[i]])>0){
            noisy=as.numeric(unlist(apply(noise.collect[[i]],1, function(x) x[1]:x[2])))
            collect=rep.int(1:nrow(noise.collect[[i]]),times=with(noise.collect[[i]],end-start+1))
            noisy.elements[[i]]=noisy.elements[[i]][noisy]
            noisy.elements[[i]][,collection:=collect]
        } else {
            noisy.elements[[i]]=data.table::data.table(start=numeric(0),end=numeric(0),state=numeric(0),sf.sum=numeric(0))
        }
    }

    ## Iterate over each region and get the sum of scale factors in active regions 
    for(i in 1:length(noisy.elements)){
        if(nrow(noisy.elements[[i]])>0){
            ## ## Get sum of scale factors in active elements for each element of noisy collections ## ##
            one.hot=tree.states[rep.int(noisy.elements[[i]]$state,times=with(noisy.elements[[i]],end-start+1)),,drop=F] ## Extract 0/1 encoding for states
            one.hot[one.hot>1]=1
            state.indices=as.numeric(unlist(apply(noisy.elements[[i]], 1, function(x) x[1]:x[2]))) ## unroll state indicies
            ## get relevant scale-factors, zero-out contributions from non-active states, then sum per call
            noisy.elements[[i]][,sf.sum:=unlist(lapply(split(do.call("cbind",lapply(scale.factors[[i]], function(x) x[state.indices] ))*one.hot, rep.int(1:length(noisy.elements[[i]]$state),times=with(noisy.elements[[i]],end-start+1))),sum))]
            
            ## ## For each collection, keep the elements that do not have the  highest score
            noisy.elements[[i]]=noisy.elements[[i]][,noise:=(sf.sum!=max(sf.sum)),by="collection"][noise==TRUE,]
            noisy.elements[[i]][,noise:=NULL]
            noisy.elements[[i]][,collection:=NULL]
        }
    }   
    return(noisy.elements)
}
