######
## Emission related code
######

###
## All code for a tree emission model
###

## Create tree emission model class
#' @export
treeEmission <- R6::R6Class("treeEmission",inherit=flexHMM::Emission)

treeEmission$set("public","updateEmissionProbabilities",function(){
    if(sum(self$getParameterTable()$replicate!=-1)>0){ stop("model doesn't allow for replicated parameters yet")}
    ## Initialize emission Log prob matrix
    cl=as.numeric(unlist(lapply(self$data,function(x) nrow(x[[1]]))))
    self$emissionLogProb=foreach(l=1:length(cl)) %do% {
        matrix(0,nrow=cl[l],ncol=self$nstates)
    }
    
    ## For each species
    for(chain in 1:length(self$data)){
        ## Compute the mu's for the active state
        si=iterators::iter(self$invariants$tree$tip.label)
        backList=foreach(s=si) %do%{
            ## Transform the mu's to make them identifiable
            return(matrix(self$params[self$getParamIndicies(paste0("mu.back.",s))],ncol=1))
        }
        si=iterators::iter(self$invariants$tree$tip.label)
        meansList=foreach(s=si) %do%{
            ## Transform the mu's to make them identifiable
           return(matrix(cumsum(c(self$params[self$getParamIndicies(paste0("mix.mu.base.",s))],self$params[self$getParamIndicies(paste0("mix.mu.",s))])),ncol=1))
        }
        si=iterators::iter(self$invariants$tree$tip.label)
        weightsList=foreach(s=si) %do% {
            ## Construct the weights via stick-breaking
            return(computeWeights(self$params[self$getParamIndicies(paste0("mix.weight.",s))]))
        }
        self$emissionLogProb[[chain]] = updateTreeEmisProbFlexMixCpp(
                                self$data[[chain]],
                                self$invariants$scaleFactors[[chain]],
                                self$invariants$deletionRanges[[chain]],
                                self$invariants$treeStates,
                                self$nstates,
                                backList,
                                meansList,
                                weightsList,
                                self$invariants$sampleNormFactors,
                                self$invariants$dispersionParams)
    }
},overwrite=TRUE)

treeEmission$set("public","checkInvariantValidity",function(invariants){
    ## Check that the number of species in the data actually equals the nuimber expected in the tree states
    nspec=length(self$data[[1]])
    if(ncol(invariants$treeStates)!=nspec){
        stop(paste("The number of species in the tree is", ncol(invariants$treeStates), "the number in the data is", nspec))
    }
    ## Check that the length of of dispersionParams
    if(!is.numeric(invariants$dispersionParams)){
        stop("dispersionParams must be a numeric vector.")
    }
    if(length(invariants$dispersionParams)!=2){
        stop("Length of dispersionParams must be two.")
    }
})

## Data masking function
treeEmission$set("public","maskData",function(maskSignalThresh=0.05){
    for(chain in 1:length(self$data)){
        for(spec in names(self$data[[chain]])){
            foo=(self$invariants$scaleFactors[[chain]][[spec]]<=maskSignalThresh)
            self$data[[chain]][[spec]][foo]=0
            self$invariants$scaleFactors[[chain]][[spec]][foo]=0
        }
    }
},overwrite=TRUE)

######
## Transition related code
######

## Tree transition with stationarity constraint
#' @export
treeTransitionSC <- R6::R6Class("treeTransitionSC",inherit=flexHMM::Transition)
## Compute tree transition probabilities (creates sparse transition matrix)
treeTransitionSC$set("public","updateTransitionProbabilities",function(){
    ## Create transition matric
    self$transitionLogProb=matrix(0,self$nstates,self$nstates)
     ## Compute the marginal probablilities for all states except the fully inactive state and use it to weight the first row of transition probabilities
    state.probs=infiniteMutProb(self$invariants$tree,self$invariants$treeStates,parms=list(inactive.prob=self$params[self$getParamIndicies("inactive.freq")],rate=self$params[self$getParamIndicies("rate")]),TRUE)
    ## Solve for rho_active parameter as function of marginal state probabilities and p_inactive (probabilities are in log space)
    ## rho.act=1-(exp(state.probs$prob[1]-logSumExp(state.probs$prob[-1])) * (1-self$params[self$getParamIndicies("autocor.inactive")]))
    rho.inact=1-((1-self$params[self$getParamIndicies("autocor.active")])*exp(logSumExp(state.probs$prob[-1])-state.probs$prob[1]))
    ## Set all sub-diagonal values in the first column
    self$transitionLogProb[-1,1]=1-self$params[self$getParamIndicies("autocor.active")]
    ## Set all elements of the diagonal to the appropriate autocorrelation
    diag(self$transitionLogProb)= self$params[self$getParamIndicies("autocor.active")]
    self$transitionLogProb[1,1]=rho.inact
    ## Set all values in the first row of the transition matrix
    sp=exp(state.probs$prob[-1]-logSumExp(state.probs$prob[-1]))
    self$transitionLogProb[1,-1]=(1-rho.inact)*sp
    self$transitionLogProb=log(self$transitionLogProb)
})

## Tree transition with no stationarity constraint
#' @export
treeTransition <- R6::R6Class("treeTransition",inherit=flexHMM::Transition)
## Compute tree transition probabilities (creates sparse transition matrix)
treeTransition$set("public","updateTransitionProbabilities",function(){
    ## Create transition matric
    self$transitionLogProb=matrix(0,self$nstates,self$nstates)
     ## Compute the marginal probablilities for all states except the fully inactive state and use it to weight the first row of transition probabilities
    state.probs=infiniteMutProb(self$invariants$tree,self$invariants$treeStates,parms=list(inactive.prob=self$params[self$getParamIndicies("inactive.freq")],rate=self$params[self$getParamIndicies("rate")]),TRUE)
    ## Set all sub-diagonal values in the first column
    self$transitionLogProb[-1,1]=1-self$params[self$getParamIndicies("autocor.active")]
    ## Set all elements of the diagonal to the appropriate autocorrelation
    diag(self$transitionLogProb)= self$params[self$getParamIndicies("autocor.active")]
    self$transitionLogProb[1,1]=self$params[self$getParamIndicies("autocor.inactive")]
    ## Set all values in the first row of the transition matrix
    sp=exp(state.probs$prob[-1]-logSumExp(state.probs$prob[-1]))
    self$transitionLogProb[1,-1]=(1-self$params[self$getParamIndicies("autocor.inactive")])*sp
    self$transitionLogProb=log(self$transitionLogProb)
})


####
## Misc. functions
####
logSumExp <- function(x){
    a=max(x)
    return(a+log(sum(exp(x-a))))
}

## Constructs set of ranges that should be considered deletions
## Ranges are 0-based and have format [,)
computeDeletionRanges <- function(scaleFactors,emptyRadius=10){
    return(lapply(scaleFactors,function(x){
        lapply(x,function(y){
            foo=rle(RcppRoll::roll_max(c(rep(0,emptyRadius),as.numeric(y),rep(0,emptyRadius)),n=2*emptyRadius+1,align="center"))
            e=cumsum(foo$lengths)[foo$values==0]
            b=e-foo$lengths[foo$values==0]
            return(cbind(b,e))
        })
    }))
}

## Build tree hmm using pre-fit data
#' @export
buildTreeFlexMixFit <- function(dataList,tree,preFitEmissions,scaleFactor,dispersionParams,sampleNormFactors,bed,binSize=100,maxChangepoints=1,delRadius=10,treeParams=NULL,stationary.constraint=FALSE,threads=1){
    tree=ape::unroot(tree)
    ## Compute number of states
    tree.states=unique(enumerateTreeStates(tree,maxChangepoints)$config)
    nstates=nrow(tree.states)
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
            bed=bed, binSize=binSize)
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
                lowerBound=list(autocor.active=0.6,autocor.inactive=0.6,inactive.freq=0.8,rate=10^-8),
                upperBound=list(autocor.active=0.98,autocor.inactive=0.99999,inactive.freq=0.999,rate=10^-1),
                invariants=list(tree=tree,treeStates=tree.states))
        } else {
            tree.trans=treeTransitionSC$new(params=list(autocor.active=estim.autocor.active,inactive.freq=estim.inact.freq,rate=10^-4),nstates=nstates,
                lowerBound=list(autocor.active=0.6,inactive.freq=0.8,rate=10^-5),
                upperBound=list(autocor.active=0.98,inactive.freq=0.999,rate=10^-1),
                invariants=list(tree=tree,treeStates=tree.states))    
        }        
    } else if(is.list(treeParams)){
        if(length(treeParams) == 4 & all(names(treeParams) %in% c("autocor.active","autocor.inactive","inactive.freq","rate"))){
            tree.trans=treeTransition$new(params=treeParams,
                nstates=nstates,
                lowerBound=list(autocor.active=0.6,autocor.inactive=0.6,inactive.freq=0.8,rate=10^-8),
                upperBound=list(autocor.active=0.98,autocor.inactive=0.99999,inactive.freq=0.999,rate=10^-1),
                invariants=list(tree=tree,treeStates=tree.states))
        } else {
            write("Incorrect tree parameter specified",stdout())
        }        
    } else {
        write("treeParams must be a list",stdout())
    }
    tree.hmm=flexHMM::HMM$new(tree.emis,tree.trans,threads=threads)
    return(tree.hmm)
}


## Misc. functions.
computeWeights <- function(weight.params){
    mix.weights=numeric(length(weight.params)+1)
    for(z in 1:length(weight.params)){
        mix.weights[z]=(1-sum(mix.weights))*weight.params[z]
    }
    mix.weights[length(mix.weights)]=1-sum(mix.weights)
    return(mix.weights)
}
