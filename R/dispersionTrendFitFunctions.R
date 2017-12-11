## Sample uniformly from the given quantiles
uniQuantSample <- function(data,quantiles,n){
    breaks =quantile(data,sort(quantiles))
    adj.quantiles=c(quantiles[!duplicated(breaks)])
    breaks=breaks[!duplicated(breaks)]
    probTable=1/table(as.numeric(cut(data,breaks,include.lowest=TRUE)))
    val.quants=probTable[cut(data,breaks,include.lowest=TRUE)]
    return(sample(which(!is.na(val.quants)),size=n,replace=TRUE,prob=val.quants[!is.na(val.quants)]))
}

## Pre-estimate dispersion using all samples using a subsample of the data sampled uniformly from quantiles
## defined by prob

#' @export
estimateJointDispersion <- function(binSums,split,samples=10000,probs=c(0,0.5,0.75,0.9,0.95,0.99,0.999,1),normFactors){
    temp=list()
    cata=unlist(unique(getExperimentTable(binSums@expDes)[,..split]))
    for(z in 1:length(cata)){
        ## combine data from all chains
        dat=do.call("rbind",lapply(splitBw(binSums,f=split,to.matrix=TRUE),function(x) x[[cata[z]]]))
        ## Subsample data to fit dispersion parameters, draw an equal number of samples from each chain
        if(ncol(dat)>1){
            temp[[cata[z]]]=dat[uniQuantSample(data=rowMeans(dat),quantiles=probs,n=samples),]
        }
    }
    dTab=data.table::data.table(split=rep(names(temp),unlist(lapply(temp,ncol))),id=unlist(lapply(temp,colnames)),normFactors=normFactors[unlist(lapply(temp,colnames))])
    temp=do.call("cbind",temp)
    dispTrendFun=fitDispersionTrend(temp,dTab)
    return(dispTrendFun)
}

## Some functions that use DESeq2 as a shortcut to calculate dispersion per expression level
fitDispersionTrend <- function(x,dTab){
    if(length(unique(dTab$split))>1){
        dds <- DESeq2::DESeqDataSetFromMatrix(countData = x,
                                              colData =dTab,
                                              design= ~ split)
    } else {
        dds <- DESeq2::DESeqDataSetFromMatrix(countData = x,
                                              colData =dTab,
                                              design= ~ 1)
    }
    DESeq2::sizeFactors(dds)=dTab$normFactors
    dds <- DESeq2::estimateDispersions(dds)
    dispFun=DESeq2::dispersionFunction(dds)
    if(attr(dispFun,"fitType")=="local"){
        ddat=data.table::data.table(x=2^(seq(1,10,0.1)),y=dispFun(2^(seq(1,10,0.1))))
        m=nls(y~a+b/x,data=ddat,start=list(a=1,b=1),lower=list(a=10^-5,b=10^-5),algorithm="port")
        out=coef(m)
        names(out)=c('asymptDisp','extraPois') 
    }else {
        out=attributes(dispFun)$coefficients
    }
    return(out)
}
