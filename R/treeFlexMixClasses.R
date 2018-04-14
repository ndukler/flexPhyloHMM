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
    ## Compute the probabilities of 0/1 at the leafs
    self$computeLeafProbabilities()
    
    ## Update the probabilities of the enumerated states
    for(chain in 1:length(self$data)){
        self$emissionLogProb[[chain]] = updateTreeEnumeratedStates(self$invariants$leafProb[[chain]],
                                self$invariants$treeStates,
                                self$nstates,
                                self$invariants$absorbing.state)
    }
} ,overwrite=TRUE)

## Function that computes the probabilities of 0/1 at each leaf
treeEmission$set("public","computeLeafProbabilities",function(){
    if(is.null(self$invariants$leafProb)){
        self$invariants$leafProb=lapply(lapply(self$data,function(x) nrow(x[[1]])),function(x) matrix(1.0,nrow=x,ncol=2*length(self$invariants$tree$tip.label)))
    }
    backList=lapply(as.list(self$invariants$tree$tip.label),function(x) matrix(1.1,nrow=1,ncol=1))
    meansList=lapply(as.list(self$invariants$tree$tip.label),function(x) matrix(1.1,nrow=self$invariants$numComponents,ncol=1))
    weightsList=lapply(as.list(self$invariants$tree$tip.label),function(x) numeric(self$invariants$numComponents))
    for(z in 1:length(self$invariants$tree$tip.label)){
        s=self$invariants$tree$tip.label[z]
        backList[[z]][1,1]=self$params[self$getParamIndicies(paste0("mu.back.",s))]
        meansList[[z]][,1]=cumsum(c(self$params[self$getParamIndicies(paste0("mix.mu.base.", s))], self$params[self$getParamIndicies(paste0("mix.mu.", s))]))
        weightsList[[z]]=computeWeights(self$params[self$getParamIndicies(paste0("mix.weight.", s))])
    }        
    for (chain in 1:length(self$data)){
        self$invariants$leafProb[[chain]] = flexPhyloHMM:::computeLeafProb(self$data[[chain]], 
                                    self$invariants$scaleFactors[[chain]], self$invariants$deletionRanges[[chain]], 
                                    backList, meansList, weightsList, self$invariants$sampleNormFactors, 
                                    self$invariants$dispersionParams)
    }    
},overwrite=TRUE)

## Function that updates the probability of the absorbing state at every step
treeEmission$set("public","forcedEmissionUpdates",function(){
    ## Only do forced updates if the model contains an absorbing state
    if(self$invariants$absorbing.state){
        ## Check that the enumerated emission state probabilities have been calculated
        if(is.null(self$emissionLogProb)){
            self$updateEmissionProbabilities()
        }
        ## Compute the probability of each enumerated tree state and the inactive tree state without data
        trans=self$hmm$transition
        pruning.params=computeASRunInvariants(trans$invariants$tree,trans$params[trans$getParamIndicies("rate")],trans$params[trans$getParamIndicies("inactive.freq")])
        ## Convert tree states to leaf log-probabilities, where if the absorbing state is present, the probability of all states on the leaves are set to 1 to calcuate the partition function
        lp=matrix(nrow=trans$nstates,ncol=2*ncol(trans$invariants$treeStates))
        lp[,seq(1,ncol(lp),2)]=abs(trans$invariants$treeStates-1)
        lp[,seq(2,ncol(lp),2)]=abs(trans$invariants$treeStates)
        lp=lp[-self$nstates,]
        lp=log(lp)
        ## Compute the probability of all the enumerated leaf patterns
        enum.tree.log.prob=with(pruning.params,treeLL(lp, qConcat = Q.concat, traversal = tt$traversal, tips = tt$tips, logBaseFreq = lbf))
        ## Compute the probability of the tree absorbing state as the complement of all the other states 
        absorb.tree.log.prob=logMinusExp(0,enum.tree.log.prob) ## Absorbing state prob
        ## For each chain update the probability of the absorbing state
        for(chain in 1:length(self$invariants$leafProb)){
            ## Compute the liklihood over all leaf labelings of the  tree given probabalistic tip labels 
            data.log.prob=with(pruning.params,treeLL(self$invariants$leafProb[[chain]], qConcat = Q.concat, traversal = tt$traversal, tips = tt$tips, logBaseFreq = lbf))
            ## Update the absorbing state probabilities
            self$emissionLogProb[[chain]][,self$nstates]=absorbingStateLogProb(self$emissionLogProb[[chain]][,-self$nstates],data.log.prob,enum.tree.log.prob,absorb.tree.log.prob)   
        }
    }
},overwrite=TRUE)

treeEmission$set("public","checkInvariantValidity",function(invariants){
    ## Check that the number of species in the data actually equals the nuimber expected in the tree states
    nspec=length(self$data[[1]])
    if(ncol(invariants$treeStates)!=nspec){
        stop(paste("The number of species in the tree is", ncol(invariants$treeStates), "the number in the data is", nspec))
    }
    ## Check that the length of of dispersionParams
    if(!is.matrix(invariants$dispersionParams) ||  ncol(invariants$dispersionParams)!=2 || nrow(invariants$dispersionParams)!=length(invariants$tree$tip.label)){
        stop("dispersionParams must be a matrix with two columns and the same number of rows as there are species.")
    }
    if(any(is.na(invariants$dispersionParams)) | any(is.nan(invariants$dispersionParams)) | any(!is.finite(invariants$dispersionParams)) | any(invariants$dispersionParams<0)){
        stop("All dispersion parameters must be finite, real values, greater than 0.")
    }
},overwrite=TRUE)

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
    ## Create transition matrix
    self$transitionLogProb=matrix(-Inf,self$nstates,self$nstates)
    ## Compute the parameter values required to run the pruning algorithm
    pruning.params=computeASRunInvariants(self$invariants$tree,self$params[self$getParamIndicies("rate")],self$params[self$getParamIndicies("inactive.freq")])
    ## Convert tree states to leaf log-probabilities, where if the absorbing state is present, the probability of all states on the leaves are set to 1 to calcuate the partition function
    lp=matrix(nrow=self$nstates,ncol=2*ncol(self$invariants$treeStates))
    lp[,seq(1,ncol(lp),2)]=abs(self$invariants$treeStates-1)
    lp[,seq(2,ncol(lp),2)]=abs(self$invariants$treeStates)
    lp[lp>1]=1
    lp=log(lp)
    ## Compute the state marginal probabilities
    if(self$invariants$absorbing.state){
        self$invariants$treeLL=with(pruning.params,treeLL(lp[-self$nstates,], qConcat = Q.concat, traversal = tt$traversal, tips = tt$tips, logBaseFreq = lbf))
        norm=log(1-exp(self$invariants$treeLL[1])) ## The normalizing factor is the L_norm=L(all states)-L(inactive state)
        enum.state.log.prob=self$invariants$treeLL[-1]-norm ## p(state_enum)=L(state_enum)/L_norm
        absorb.state.log.prob=log(1-exp(logSumExp(enum.state.log.prob))) ## Absorbing state prob
        state.log.prob=c(enum.state.log.prob,absorb.state.log.prob)
    } else {
        self$invariants$treeLL=with(pruning.params,treeLL(lp, qConcat = Q.concat, traversal = tt$traversal, tips = tt$tips, logBaseFreq = lbf))
        norm=logSumExp(self$invariants$treeLL[-1]) ## The normalizing factor is the sum of all enumerated states
        state.log.prob=self$invariants$treeLL[-1]-norm ## p(state_enum)=L(state_enum)/L_norm
    }
    ## Solve for rho_active parameter as function of marginal state probabilities and p_inactive (probabilities are in log space)
    rho.inact=1-((1-self$params[self$getParamIndicies("autocor.active")])*exp(logSumExp(self$invariants$treeLL[-1])-exp(self$invariants$treeLL[1])))
    ## Set all sub-diagonal values in the first column
    self$transitionLogProb[-1,1]=log(1-self$params[self$getParamIndicies("autocor.active")])
    ## Set all elements of the diagonal to the appropriate autocorrelation
    diag(self$transitionLogProb)= log(self$params[self$getParamIndicies("autocor.active")])
    ## Set all values in the first row of the transition matrix
    self$transitionLogProb[1,-1]=log(1-rho.inact)+state.log.prob
    self$transitionLogProb[1,1]=log(rho.inact)
})

## Tree transition with no stationarity constraint
#' @export
treeTransition <- R6::R6Class("treeTransition",inherit=flexHMM::Transition)
## Compute tree transition probabilities (creates sparse transition matrix)
treeTransition$set("public","updateTransitionProbabilities",function(){
    ## Create transition matrix
    self$transitionLogProb=matrix(-Inf,self$nstates,self$nstates)
    ## Compute the parameter values required to run the pruning algorithm
    pruning.params=computeASRunInvariants(self$invariants$tree,self$params[self$getParamIndicies("rate")],self$params[self$getParamIndicies("inactive.freq")])
    ## Convert tree states to leaf log-probabilities, where if the absorbing state is present, the probability of all states on the leaves are set to 1 to calcuate the partition function
    lp=matrix(nrow=self$nstates,ncol=2*ncol(self$invariants$treeStates))
    lp[,seq(1,ncol(lp),2)]=abs(self$invariants$treeStates-1)
    lp[,seq(2,ncol(lp),2)]=abs(self$invariants$treeStates)
    lp[lp>1]=1
    lp=log(lp)
    ## Compute the state marginal probabilities
    if(self$invariants$absorbing.state){
        self$invariants$treeLL=with(pruning.params,treeLL(lp[-self$nstates,], qConcat = Q.concat, traversal = tt$traversal, tips = tt$tips, logBaseFreq = lbf))
        norm=log(1-exp(self$invariants$treeLL[1])) ## The normalizing factor is the L_norm=L(all states)-L(inactive state)
        enum.state.log.prob=self$invariants$treeLL[-1]-norm ## p(state_enum)=L(state_enum)/L_norm
        absorb.state.log.prob=log(1-exp(logSumExp(enum.state.log.prob))) ## Absorbing state prob
        state.log.prob=c(enum.state.log.prob,absorb.state.log.prob)
    } else {
        self$invariants$treeLL=with(pruning.params,treeLL(lp, qConcat = Q.concat, traversal = tt$traversal, tips = tt$tips, logBaseFreq = lbf))
        norm=logSumExp(self$invariants$treeLL[-1]) ## The normalizing factor is the sum of all enumerated states
        state.log.prob=self$invariants$treeLL[-1]-norm ## p(state_enum)=L(state_enum)/L_norm
    }    
    ## Set all sub-diagonal values in the first column
    self$transitionLogProb[-1,1]=log(1-self$params[self$getParamIndicies("autocor.active")])
    ## Set all elements of the diagonal to the appropriate autocorrelation
    diag(self$transitionLogProb)= log(self$params[self$getParamIndicies("autocor.active")])
    self$transitionLogProb[1,1]=log(self$params[self$getParamIndicies("autocor.inactive")])
    ## Set all values in the first row of the transition matrix
    self$transitionLogProb[1,-1]=log(1-self$params[self$getParamIndicies("autocor.inactive")])+state.log.prob
})


####
## Misc. functions
####

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

computeWeights <- function(weight.params){
    mix.weights=numeric(length(weight.params)+1)
    for(z in 1:length(weight.params)){
        mix.weights[z]=(1-sum(mix.weights))*weight.params[z]
    }
    mix.weights[length(mix.weights)]=1-sum(mix.weights)
    return(mix.weights)
}
