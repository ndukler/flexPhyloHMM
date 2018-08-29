#' @export
writePAF <- function (hmm, marginal = NULL, viterbi = NULL, marginalThresh = 0.5) {
    if (is.null(marginal) & is.null(viterbi)) {
        stop("Either marginal or viterbi path must be passed to function.")
    }
    if (!is.null(marginal) & !is.null(viterbi)) {
        stop("Cannot pass both viterbi and marginal to function")
    }
    if (!is.null(marginal)) {
        states = lapply(marginal, function(x) unlist(apply(x >= 
            marginalThresh, 1, which)))
        eList=lapply(states, function(x) x>1)
    } else {
        eList=lapply(viterbi, function(x) x>1)
    }
    fp=file.path(hmm$logDir,"conditionalAlleleProbabilities")
    calcPAF(eList,hmm$emission$invariants$leafProb,hmm$emission$invariants$tree$tip.label,fp,
            hmm$emission$invariants$bed$chrom,hmm$emission$invariants$bed$start,hmm$emission$invariants$binSize)    
}
