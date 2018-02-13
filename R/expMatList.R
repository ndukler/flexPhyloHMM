#' Compute exponentiated rate matrix for all branches
#'
#' @param rate a mutation rate
#' @param branch.lengths the branch lengths of interes
#' @param log if TRUE returns log transition matricies
#' @return a list of transition matricies, one entry for each specified branch length
#' @export
expMatList <- function(rate,base.freq.zero,branch.lengths,log=TRUE){
    bf=c(base.freq.zero,1-base.freq.zero)
    ## construct the scaled rate matrix
    Q=rate*matrix(c(-1/(2*bf[1]),1/(2*bf[1]),1/(2*bf[2]),-1/(2*bf[2])),nrow = 2,ncol=2,byrow = TRUE)
    Q.list=lapply(as.list(branch.lengths),function(x) as.matrix(Matrix::expm(Q*x)))
    if(log==TRUE)
        return(lapply(Q.list,log))
    else
        return(Q.list)
}
