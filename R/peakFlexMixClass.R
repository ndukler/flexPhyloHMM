###
## All code for a negative binomial emission model where the active componant is a mixture model with means
## set by weights from discretized gamma distribution
###

## Create negative binomial emission model class

#' @export
flexMixNBEmission <- R6::R6Class("hmm.emission",inherit=flexHMM::Emission)

## Custom functions for checking validity of inputs, called by constructor
flexMixNBEmission$set("public","checkParamValidity",function(params){
    validParams=c("mix.weight","mix.mu","mu.back","mix.mu.base")
    if(sum(!names(params) %in% c(validParams))>0){
        stop(paste("There is an invalid parameter. Valid parameter names are:",paste(validParams,collapse=",")))
    }
    if(sum(validParams %in% names(params)) != length(validParams)){
        stop(paste("There is a missing parameter. Valid parameter names are:",paste(validParams,collapse=",")))
    }
})

flexMixNBEmission$set("public","checkInvariantValidity",function(invariants){
    validInvariants=c("scaleFactors","dispParams","sampleNormFactors")
    if(sum(!names(invariants) %in% c(validInvariants))>0){
        stop(paste("There is an invalid invariant. Valid invariant names are:",paste(validInvariants,collapse=",")))
    }
    if(sum(validInvariants %in% names(invariants)) != length(validInvariants)){
        stop(paste("There is a missing invariant. Valid invariant names are:",paste(validInvariants,collapse=",")))
    }
    if(length(invariants$scaleFactors)!=length(self$data)){
        stop(paste0("The length of scaleFactors(",length(invariants$scaleFactors),") must be the same as the length of the data(",length(self$data),")"))
    }
    for(i in 1:length(self$data)){
        if(!is.numeric(invariants$scaleFactors[[i]])){
            stop(paste("Element",i,"of scaleFactors is not numeric."))
        }
        if(nrow(self$data[[i]][[1]])!=length(invariants$scaleFactors[[i]])){
            stop(paste("Scale factor length and chain length are not the same on chain",i))
        }
    }
})

flexMixNBEmission$set("public","updateEmissionProbabilities",function(){
    if(sum(self$getParameterTable()$replicate!=-1)>0){ stop("model doesn't allow for replicated parameters yet")}
    ## iterate over all states
    self$emissionLogProb=foreach::foreach(i=1:length(self$data),.export=c("updatePeakFlexMixCpp","computeWeights")) %do% {
        ## Transform the mu's to make them identifiable
        mm=matrix(cumsum(c(self$params[self$getParamIndicies("mix.mu.base")],self$params[self$getParamIndicies("mix.mu")])),ncol=1)
        ## Construct the weights via stick-breaking
        mw=matrix(computeWeights(self$params[self$getParamIndicies("mix.weight")]),ncol=1)
        ## Return the logLiklihood
        return(updatePeakFlexMixCpp(self$data[[i]][[1]],
                                                       self$invariants$scaleFactors[[i]],
                                                       self$nstates,
                                                       matrix(self$params[self$getParamIndicies("mu.back")],ncol=1),
                                                       mm,
                                                       mw,
                                                       self$invariants$dispParams,
                                                       self$invariants$sampleNormFactors))
    }
})

## Mask out data with small scale factor
flexMixNBEmission$set("public","maskData",function(maskSignalThresh=0.05){
    for(chain in 1:length(self$data)){
        foo=(self$invariants$scaleFactors[[chain]]<=maskSignalThresh)
        self$data[[chain]][[1]][foo]=0
        self$invariants$scaleFactors[[chain]][foo]=0
    }
},overwrite=TRUE)


###
## HMM with negative binomial emissions, but explicitly only two states (peak and background) 
###

###
## All code for simple two state transition model
###

## Create two state transition model

#' @export
simpleTwoStateTransition <- R6::R6Class("simpleTwoStateTransition",inherit=flexHMM::Transition)

simpleTwoStateTransition$set("public","updateTransitionProbabilities",function(){
    ## iterate over all states
    self$transitionLogProb=log(matrix(c(self$params[self$getParamIndicies("s11",1)],1-self$params[self$getParamIndicies("s11",1)],
                                        1-self$params[self$getParamIndicies("s22",1)],self$params[self$getParamIndicies("s22",1)]),ncol=2,byrow=TRUE))
})

simpleTwoStateTransitionLog <- R6::R6Class("simpleTwoStateTransitionLog",inherit=flexHMM::Transition)

simpleTwoStateTransitionLog$set("public","updateTransitionProbabilities",function(){
    ## iterate over all states
    self$transitionLogProb=matrix(c(self$params[self$getParamIndicies("s11",1)],log(1-exp(self$params[self$getParamIndicies("s11",1)])),
                                        log(1-exp(self$params[self$getParamIndicies("s22",1)])),self$params[self$getParamIndicies("s22",1)]),ncol=2,byrow=TRUE)
})

## Misc. functions.

#' @export
computeWeights <- function(weight.params){
    mix.weights=numeric(length(weight.params)+1)
    for(z in 1:length(weight.params)){
        mix.weights[z]=(1-sum(mix.weights))*weight.params[z]
    }
    mix.weights[length(mix.weights)]=1-sum(mix.weights)
    return(mix.weights)
}
