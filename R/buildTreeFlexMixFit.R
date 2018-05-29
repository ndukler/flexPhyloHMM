## Build tree hmm using pre-fit data
#' @export
buildTreeFlexMixFit <- function(dataList,tree,preFitEmissions,scaleFactor,sampleNormFactors,bed,binSize=100,maxChangepoints=1,delRadius=10,treeParams=NULL,stationary.constraint=FALSE,threads=1,absorbing.state=FALSE,log.dir="."){
    tree=ape::unroot(tree)
    ## Compute number of states
    tree.states=unique(enumerateTreeStates(tree,maxChangepoints,absorbing.state)$config)
    nstates=nrow(tree.states)
    ## Check if all possible states have been enumerated and set absorbing state to FALSE if they have been
    if(absorbing.state && nstates > 2^length(tree$tip.label)){
        warning("All possible states have been enumerated, setting absorbing.state to FALSE")
        absorbing.state=FALSE
        tree.states=unique(enumerateTreeStates(tree,maxChangepoints,absorbing.state)$config)
        nstates=nrow(tree.states)
    }
    
    ## Build emission parameter list for each species
    param=unlist(lapply(as.list(tree$tip.label),function(x){
        ncomp=nrow(preFitEmissions[[x]][paramName=="mix.weight"])+1
        z=c(8,10,rep(5,ncomp-1),rep(0.5,ncomp-1))
        names(z)=c(paste0("mu.back.",x),paste0("mix.mu.base.",x),rep(paste0("mix.mu.",x),ncomp-1),rep(paste0("mix.weight.",x),ncomp-1))
        return(z)
    }))
    fixed=unlist(lapply(as.list(tree$tip.label),function(x){
        ncomp=nrow(preFitEmissions[[x]][paramName=="mix.weight"])+1
        z=c(TRUE,TRUE,rep(TRUE,ncomp-1),rep(TRUE,ncomp-1))
        names(z)=c(paste0("mu.back.",x),paste0("mix.mu.base.",x),rep(paste0("mix.mu.",x),ncomp-1),rep(paste0("mix.weight.",x),ncomp-1))
        return(z)
    }))
    param=split(param,names(param))
    fixed=split(fixed,names(fixed))
    ## Get number components
    numComponents=nrow(preFitEmissions[[1]][paramName=="mix.weight"])+1
    ## Key sample normalization factors
    setkey(sampleNormFactors,"experiment.id")
    ## Build matrix of parameters for dispersion function, each row corresponds to species (same order as tree leaves), each column
    dispersionParams=matrix(0,nrow=length(tree$tip.label),ncol=2)
    for(i in 1:length(tree$tip.label)){
        dispersionParams[i,1]=preFitEmissions[[tree$tip.label[i]]][paramName=="asymptDisp"]$value
        dispersionParams[i,2]=preFitEmissions[[tree$tip.label[i]]][paramName=="extraPois"]$value
    }
    ## Build emission object - make sure that data is in the same order as the tree
    tree.emis=treeEmission$new(
        data=lapply(dataList, function(y){
            out=lapply(as.list(tree$tip.label),function(x) y[[x]])
            names(out)=tree$tip.label
            return(out)
        }) ,nstates=nstates, params=param,fixed=fixed,
        invariants=c(tree=tree,treeStates=tree.states,
            scaleFactors= lapply(scaleFactor, function(y){
                out=lapply(as.list(tree$tip.label),function(x) y[[x]])
                names(out)=tree$tip.label
                return(out)
            }),
            dispersionParams=dispersionParams,
            numComponents=numComponents,
            sampleNormFactors=lapply(as.list(tree$tip.label), function(x) {
                sampleNormFactors[colnames(dataList[[1]][[x]])]$norm.factors
            }),
            bed=bed, binSize=binSize,absorbing.state=absorbing.state)
        )
    tree.emis$invariants$deletionRanges=computeDeletionRanges(tree.emis$invariants$scaleFactor,delRadius)
    ## Set emission object parameters
    for(p in names(preFitEmissions)){
        subT=preFitEmissions[[p]][pType=="emission"]
        setkeyv(subT,c("paramName","replicate"))
        for(r in unique(subT$replicate)){
            tree.emis$setParamValue(paste0("mu.back.",p),value=subT[.("mu.back",r)]$value,chain=-1,replicate=r)
            tree.emis$setParamValue(paste0("mix.mu.",p),value=subT[.("mix.mu",r)]$value,chain=-1,replicate=r)
            tree.emis$setParamValue(paste0("mix.mu.base.",p),value=subT[.("mix.mu.base",r)]$value,chain=-1,replicate=r)
            tree.emis$setParamValue(paste0("mix.weight.",p),value=subT[.("mix.weight",r)]$value,chain=-1,replicate=r)
        }
    }
    
    ## If a list of tree parameters is not supplied
    if(is.null(treeParams)){
        ## Provide initial estimates of the auto-correlation parameters from pre-fit HMMs
        estim.autocor.inactive=mean(unlist(lapply(preFitEmissions,function(x) x[paramName=="s11"]$value)))
        estim.autocor.active=mean(unlist(lapply(preFitEmissions,function(x) x[paramName=="s22"]$value)))
        ## Use the auto-correlation factors to get the fraction of bases which are active
        rl.inact=1/(1-estim.autocor.inactive)
        rl.act=1/(1-estim.autocor.active)
        estim.inact.freq=rl.inact/(rl.inact+rl.act)
        ## Build objects and return HMM
        if(!stationary.constraint){
            tree.trans=treeTransition$new(params=list(autocor.active=estim.autocor.active,autocor.inactive=estim.autocor.inactive,inactive.freq=estim.inact.freq,rate=10^-4),
                nstates=nstates,
                lowerBound=list(autocor.active=0.1,autocor.inactive=0.1,inactive.freq=0.001,rate=10^-8),
                upperBound=list(autocor.active=0.99999,autocor.inactive=0.99999,inactive.freq=0.999,rate=10^3),
                invariants=list(tree=tree,treeStates=tree.states,absorbing.state=absorbing.state,max.changepoints=maxChangepoints))
        } else {
            tree.trans=treeTransitionSC$new(params=list(autocor.active=estim.autocor.active,inactive.freq=estim.inact.freq,rate=10^-4),nstates=nstates,
                lowerBound=list(autocor.active=0.01,inactive.freq=10^-3,rate=10^-8),
                upperBound=list(autocor.active=0.99999,inactive.freq=0.99999,rate=10^2),
                invariants=list(tree=tree,treeStates=tree.states,absorbing.state=absorbing.state,max.changepoints=maxChangepoints))    
        }        
    } else if(is.list(treeParams)){ ## If a list of tree parameters is supplied
        ## Figure out if the parameter set matches the transition model with or without a stationary constraint
        if(length(treeParams) == 4 & all(names(treeParams) %in% c("autocor.active","autocor.inactive","inactive.freq","rate"))){
            tree.trans=treeTransition$new(params=treeParams,
                nstates=nstates,
                lowerBound=list(autocor.active=0.6,autocor.inactive=0.6,inactive.freq=0.8,rate=10^-8),
                upperBound=list(autocor.active=0.98,autocor.inactive=0.99999,inactive.freq=0.999,rate=10^-1),
                invariants=list(tree=tree,treeStates=tree.states,absorbing.state=absorbing.state,max.changepoints=maxChangepoints))
        } else if(length(treeParams) == 3 & all(names(treeParams) %in% c("autocor.active","inactive.freq","rate"))){
            tree.trans=treeTransitionSC$new(params=treeParams,nstates=nstates,
                lowerBound=list(autocor.active=0.01,inactive.freq=10^-3,rate=10^-8),
                upperBound=list(autocor.active=0.99999,inactive.freq=0.99999,rate=10^2),
                invariants=list(tree=tree,treeStates=tree.states,absorbing.state=absorbing.state,max.changepoints=maxChangepoints))    
        } else {
            write("Incorrect tree parameter specified",stdout())
        }        
    } else {
        write("treeParams must be a list",stdout())
    }
    tree.hmm=flexHMM::HMM$new(tree.emis,tree.trans,threads=threads,log.dir=log.dir)
    return(tree.hmm)
}
