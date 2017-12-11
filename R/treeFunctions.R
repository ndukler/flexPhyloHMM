## Modified from http://blog.phytools.org/2012/01/function-to-get-descendant-node-numbers.html
getDescendants<-function(tree,node,curr=NULL,onlyTips=FALSE){
    if(is.null(curr)) curr<-vector()
    daughters<-tree$edge[which(tree$edge[,1]==node),2]
    curr<-c(curr,daughters)
    w<-which(daughters>=length(tree$tip))
    if(length(w)>0) for(i in 1:length(w)) 
        curr<-getDescendants(tree,daughters[w[i]],curr)
    if(onlyTips)
        curr=curr[curr %in% 1:length(tree$tip.label)]
    return(curr)
}

## Produces a list that contains:
## 1) A matrix that shows all possible states given a number of "changepoints" along the tree
## 2) A matrix that records the products of the edge lengths of the branches the changepoints  occur and the
##     corresponding number of  changepoints
enumerateTreeStates <-  function(tree,changepoints=1){
    if(changepoints < 0 | changepoints %% 1 > 0){
        stop("changepoints must be an integer and >= 0")
    }
    if(ape::is.rooted(tree)){
        warning("Tree cannot be rooted for this operation, unrooting tree")
        tree=ape::unroot(tree)
    }
    el=list()
    leafConfig=list()
    ## Create two configurations for 0-mutation events
    leafConfig[[1]]=matrix(c(0,1),nrow=2,ncol=length(tree$tip.label))
    el[[1]]=data.table::data.table(bl.prod=c(0,0),cp.count=c(0,0))
    if(changepoints >0){
        ## Add fully active state
        for(i in 1:changepoints){
            ## Generate all combinations of changepoints
            cp.edge=t(combn(1:nrow(tree$edge),i))
            ## Record all the edge lengths associated with each event
            el[[i+1]]=data.table::data.table(bl.prod=rep(apply(cp.edge,1, function(x) prod(tree$edge.length[x])),2),cp.count=c(i,i))
            ## Iterate over all combinations
            leafConfig[[as.character(i)]]=matrix(0,nrow=nrow(cp.edge)*2,ncol=length(tree$tip.label))
            for(x in 1:nrow(cp.edge)){
                cmd=paste0("table(factor(c(tree$edge[cp.edge[x,],2],",paste0("getDescendants(tree,",tree$edge[cp.edge[x,],2],",onlyTips=TRUE)",collapse=","),"),levels=1:length(tree$tip.label)))")
                switchPerLeaf=as.numeric(eval(parse(text=cmd)))
                leafConfig[[i+1]][2*x-1,]=switchPerLeaf %% 2
                leafConfig[[i+1]][2*x,]=(switchPerLeaf+1) %% 2           
            }
        }        
    }
    leafConfig=do.call("rbind",leafConfig)
    el=data.table::rbindlist(el)
    rorder=order(rowSums(leafConfig))
    el=el[rorder]
    leafConfig=leafConfig[rorder,]
     return(list(config=leafConfig,edgeLengthProd=el))
}

## Compute probability of each tree mutational state
infiniteMutProb <- function(tree,tree.states,parms=list(inactive.prob=0.95,rate=10^-3),log.prob=TRUE){
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
    setnames(tree.states,colnames(tree.states),tree$tip.label)
    ## Create phyDat alignment objects
    pd=phangorn::phyDat(tree.states,type="USER",levels=c(0,1))
    ## Estimate liklihood of each tree
    tree.fit=phangorn::pml(tree,data=pd,rate=parms$rate,bf=c(parms$inactive.prob,1-parms$inactive.prob))
    ## Convert to probability
    state.log.probs=with(tree.fit,siteLik-(max(siteLik)+log(sum(exp(siteLik-max(siteLik))))))
    if(log.prob){
        return(list(tree=tree,states=tree.states,prob=state.log.probs))
    } else {
        return(list(tree=tree,states=tree.states,prob=exp(state.log.probs)))
    }    
}
