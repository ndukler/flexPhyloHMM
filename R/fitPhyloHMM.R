rhoI_Ineq <- function(x,hmmObj,scale=1){
    if (length(x) != sum(!hmmObj$transition$fixed) + sum(!hmmObj$emission$fixed)) {
        stop("Length of x does not match the number of expected parameters")
    }
    ## Save parameter values
    eparms=hmmObj$emission$params
    tparms=hmmObj$transition$params
    
    ## Overwrite object parameters with new parameter values
    hmmObj$emission$params[!hmmObj$emission$fixed] = x[1:sum(!hmmObj$emission$fixed)]
    hmmObj$transition$params[!hmmObj$transition$fixed] = x[1:sum(!hmmObj$transition$fixed) + 
                                sum(!hmmObj$emission$fixed)]
    ## Check that 0<rhoI<1
    pruning.params = computeASRunInvariants(hmmObj$transition$invariants$tree, 
        hmmObj$transition$params[hmmObj$transition$getParamIndicies("rate")], hmmObj$transition$params[hmmObj$transition$getParamIndicies("inactive.freq")])
    lp = matrix(rep(c(0,-Inf),ncol(hmmObj$transition$invariants$treeStates)),nrow = 1)
    p.so=with(pruning.params, flexPhyloHMM:::treeLL(lp, qConcat = Q.concat, traversal = tt$traversal, 
        tips = tt$tips, logBaseFreq = lbf))
    rho.inact = 1 - ((1 - hmmObj$transition$params[hmmObj$transition$getParamIndicies("autocor.active")]) * 
        ((1 - exp(p.so))/exp(p.so)))
    ## Write constraint status out to log file
    fp=file.path(hmmObj$logDir,"parameterUpdateTimer.txt")
    write("---Evaluating Constraint status---",file=fp,append=TRUE)
    write(paste("Transition parameter vector", paste(hmmObj$transition$params,collapse=",")),file=fp,append=TRUE)
    write(paste0("rho.inact=", rho.inact),file=fp,append=TRUE)
    ## Restore original parameter values
    hmmObj$emission$params=eparms
    hmmObj$transition$params=tparms
    ## return rho.inact
    return(rho.inact) 
}

#' @export
fitPhyloHMM <- function(hmm,nthreads=1){
    ## Clear the parameter update watcher for this round
    fp=file.path(hmm$logDir,"parameterUpdateTimer.txt")
    write("***Beginning optimization***\n", file=fp,append=TRUE)
    start.time=Sys.time()
    write(paste0("Start time= ",start.time), file=fp,append=TRUE)
    nsites=sum(unlist(lapply(hmm$emission$data,function(x) nrow(x[[1]]))))
    ## Pass non-fixed parameters for optimization
    tryCatch({
        sink(file=file.path(hmm$logDir,"optim_log.txt"),append=TRUE)
        write(paste0(paste0(rep("-",30),collapse=""),"START",paste0(rep("-",30),collapse="")),stdout())
        final.params=Rsolnp::solnp(fun=updateAllParams, pars=c(hmm$emission$params[!hmm$emission$fixed],hmm$transition$params[!hmm$transition$fixed]),
            LB = c(hmm$emission$lowerBound[!hmm$emission$fixed],hmm$transition$lowerBound[!hmm$transition$fixed]),
            UB = c(hmm$emission$upperBound[!hmm$emission$fixed],hmm$transition$upperBound[!hmm$transition$fixed]),
            ineqfun=rhoI_Ineq, ineqLB=0.05,ineqUB=1,
            hmmObj=hmm,scale=-1/nsites,control=list(delta=10^-6,tol=10^-7,rho=5))
        hmm$convergence=final.params$convergence
        write(paste0(paste0(rep("-",30),collapse=""),"END",paste0(rep("-",30),collapse="")),stdout())
        sink()
    }, interrupt=function(i){            
        write(paste0(paste0(rep("-",30),collapse=""),"USER TERMINATED",paste0(rep("-",30),collapse="")),stdout())
        sink()
            warning("HMM optimization interrupted by user.")
    }, error = function(e){
        save(hmm,file=file.path(hmm$logDir,"hmm.bin"))
        sink()
        stop(e)
    })
    ## Update hmm to use final parameters and update logLiklihood
    if(sum(!hmm$emission$fixed)>0){
        hmm$emission$params[!hmm$emission$fixed]=final.params$par[1:sum(!hmm$emission$fixed)]
        hmm$emission$updateEmissionProbabilities()
    }
    if(sum(!hmm$transition$fixed)>0){
        hmm$transition$params[!hmm$transition$fixed]=final.params$par[1:sum(!hmm$transition$fixed)+sum(!hmm$emission$fixed)]
        hmm$transition$forcedTransitionUpdates()
    }
    ## Run any additional forced updates to the emission or transition matricies
    hmm$emission$forcedEmissionUpdates()
    hmm$transition$forcedTransitionUpdates()
    hmm$transition$updatePrior()
    ## Run forward algorithm
    hmm$forwardAlgorithm()
    ## Compute log-likelihood
    hmm$computeLogLiklihood()
    end.time=Sys.time()
    write(paste0("End time= ",end.time), file=fp,append=TRUE)
    write(end.time-start.time, file=fp,append=TRUE)
    write("------------------------------------------------------------------------------------------------------", file=fp,append="TRUE")    
}
