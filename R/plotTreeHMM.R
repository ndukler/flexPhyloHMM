## Now some methods to plot HMM

#' @export
plot.tree.hmm <- function(hmm,viterbi=NULL,marginal=NULL,calls.bed=NULL,geneAnnotations=NULL,start=NA_real_,end=NA_real_,chain=1,dat.min=1,dat.max=50,plot.type="dot",show.bins=TRUE,include.scale.factors=TRUE,gene.type.filter=NULL){
              ## If start and end not set, use full data length as default
              if(is.na(start))
                  start=1
              if(is.na(end))
                  end=nrow(hmm$emission$emissionLogProb[[chain]])
              if(start>end){
                  stop("Start must be less than end")
              }
              if(nrow(hmm$emission$emissionLogProb[[chain]])<end){
                  stop("End is greater than length of requested chain")
              }
              if(start<1){
                  stop("Start must be greater than 0.")
              }             
              ## Check that gene models are in GRangesList if not null
              if(!is.null(geneAnnotations) & class(geneAnnotations)[1]!="GRangesList"){
                  ## stop("GeneAnnotations must be in GRangesList format for ggbio.")
              }
              ## Calculate genomic position mapping
              chrom=hmm$emission$invariants$bed[chain,]$chrom
              chrom.start=hmm$emission$invariants$bed[chain,]$start
              chrom.end=hmm$emission$invariants$bed[chain,]$end
              chrom.loci=seq(chrom.start+thmm$emission$invariants$binSize/2,chrom.end-thmm$emission$invariants$binSize/2,by=thmm$emission$invariants$binSize) ## create chrom translation so corrdinates are centered on bin
              chrom.lim=chrom.loci[c(start,end)]+c(-thmm$emission$invariants$binSize/2,thmm$emission$invariants$binSize/2)
              ## Create x-axis labels
              loci.breaks=labeling::extended(chrom.lim[1],chrom.lim[2], m = 5)
              if(show.bins){
                  bin.labs=scales::comma(round((loci.breaks-chrom.start)/hmm$emission$invariants$binSize+1))
                  bin.labs[bin.labs<=0]="OOR"
                  x.labs=paste0(scales::comma(loci.breaks),"\n(",bin.labs,")")
              } else {
                  x.labs=scales::comma(loci.breaks)
              }
              
              ## Set tick width and positions
              num.species=length(hmm$emission$invariants$tree$tip.label)
              x.tick.width=100*thmm$emission$invariants$binSize
              xt=data.table(
                  x=seq(chrom.start,chrom.end,by=x.tick.width),
                  y=0.5,
                  yend=num.species+0.5)
              yt=data.table(
                  x=chrom.loci[start],
                  xend=chrom.loci[end],
                  y=seq(0.5, num.species+0.5, by=1))
              
              ## Create plot list
              plist=list()
              ## If geneAnnotations are present
              if(!is.null(geneAnnotations)){
                  dbcon <- RSQLite::dbConnect(RSQLite::SQLite(),dbname=geneAnnotations)
                  dt=data.table::as.data.table(RSQLite::dbGetQuery(dbcon, "SELECT * FROM transcript
                                            INNER JOIN gene ON transcript._tx_id = gene._tx_id
                                           WHERE tx_chrom= :chr AND ((tx_start > :start AND tx_start < :end) OR (tx_end > :start AND tx_end < :end))",
                      params=data.frame(chr=chrom,start=chrom.lim[1],end=chrom.lim[2])))
                  if(nrow(dt)>0){
                      ## Extract total region covered by gene and clip to query region
                      dt=dt[,.(chrom=tx_chrom[1],start=max(chrom.lim[1],min(tx_start)),end=min(chrom.lim[2],max(tx_end)),gene_name=gene_name[1],strand=tx_strand[1],gene_type=gene_type[1]),by="gene_id"]
                      ## Filter on gene type if desired
                      if(!is.null(gene.type.filter)){
                          dt=dt[gene_type %in% gene.type.filter]
                      }
                      ## Convert to IRanges to check for overlaps
                      gtf=with(dt,IRanges::IRanges(start,end))
                      overlaps=data.table::as.data.table(IRanges::findOverlaps(gtf))
                      tracks=min(max(overlaps[,length(subjectHits),by="queryHits"]$V1,2),7)
                      
                      dt[,tracks:=1:tracks]
                      offset=rep(c(0.3,-0.3),nrow(dt))
                      dt[,name.pos:=offset[1:length(gene_name)],by="tracks"]
                      
                      g.gene <- ggplot2::ggplot(dt,ggplot2::aes(ymin=tracks - 0.16,ymax=tracks+0.16,xmin=start,xmax=end,color=gene_type,fill=gene_type))+
                      ggplot2::geom_rect()+
                          cowplot::theme_cowplot()+
                          ggplot2::geom_text(ggplot2::aes(x=(start+end)/2,y=tracks+name.pos,label=paste0(gene_name,"(",strand,")")))+
                          ggplot2::guides(fill=ggplot2::guide_legend(title="Gene Type"),color=FALSE)+
                          ggplot2::theme(line = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank() , title = ggplot2::element_blank())
                      plist[[1]]=g.gene
                  } else {
                      plist[[1]]=ggplot2::ggplot()
                  }

              }              
              ## Add elements one layer at a time
              ## First create call plots
                 ## If final calls are passed to plotter in bed form
              if(!is.null(calls.bed)){
                  ## pull out the calls that are in the region of interest
                  sub.calls=calls.bed[(start >= chrom.lim[1] & start <= chrom.lim[2]) | (end >= chrom.lim[1] & end  <= chrom.lim[2]) & !grepl("#",chrom) & chrom==chrom]
                  ## Translate coordinates to bin locations
                  sub.calls[,start.bin:=(start-chrom.start)/thmm$emission$invariants$binSize+1]
                  sub.calls[,end.bin:=(end-chrom.start)/thmm$emission$invariants$binSize+1]
                  ## Create projection onto species space
                  if(nrow(sub.calls)>0){
                      spec=data.table::data.table(hmm$emission$invariants$treeStates[sub.calls$state,, drop = FALSE])
                      data.table::setnames(spec,colnames(spec),hmm$emission$invariants$tree$tip.label)
                      spec=reshape2::melt(cbind(sub.calls[,.(start,end)],spec),id.vars=c("start","end"))
                      spec=spec[value==1]
                      spec[,ymin:=dat.min]
                      spec[,ymax:=dat.max]
                      data.table::setnames(spec,"variable",".id")
                  } else {
                      spec=data.table::data.table()
                  }
              } else {
                  spec=data.table::data.table()
              }
              
              ## Then create data plots
              if(!is.null(hmm$emission$emissionLogProb)){
                  ## Plot data
                  dat=rbindlist(lapply(hmm$emission$data[[chain]],function(x){
                      x=data.table::as.data.table(x)[start:end]
                      x[,index:=start:end]
                      x[,loci:=chrom.loci[index]]
                      melt(x,id.vars=c("index","loci"))
                  }),idcol=TRUE)
                  dat[,variable:=NULL]
                  dat[,type:="Raw data"]
                  if(is.numeric(dat.max)){
                      dat[value>dat.max,value:=dat.max]
                  }
                  dat.max=max(dat$value,dat.max,na.rm=TRUE)
                  ## dat.by=(log(dat.max)-log(dat.min))/5
                  ## dat.seq=round(exp(seq(log(dat.min),log(dat.max),dat.by)),2)
                  dat.by=(dat.max-dat.min)/5
                  dat.seq=round(seq(dat.min,dat.max,dat.by),2) 
                  if(plot.type=="heatmap"){
                      dat[,.id:=factor(as.character(.id),levels=names(hmm$emission$data[[chain]]))]
                      g.dat=ggplot2::ggplot()+geom_raster(data=dat,ggplot2::aes(x=loci,y=.id,fill=value),inherit.aes=FALSE)+
                          ggplot2::scale_fill_gradient(limits=c(dat.min,dat.max),labels=as.character(dat.seq),breaks=dat.seq,low="#C1E7FA", high="#062F67",na.value="white")+
                          ggplot2::ylab("Species")+
                          ggplot2::scale_x_continuous(name=chrom,breaks=loci.breaks,labels=x.labs,limits=c(chrom.lim[1],chrom.lim[2]))+
                          ggplot2::guides(fill=ggplot2::guide_colorbar(title="Trait Value"))+
                          cowplot::theme_cowplot()+
                              ggplot2::theme(axis.title.y=ggplot2::element_blank())
                      pdf(NULL)
                      g.dat <- gtable::gtable_add_cols(ggplot2::ggplotGrob(g.dat), unit(2,"cm"),0)
                      g.dat <- gtable::gtable_add_grob(g.dat, ggplot2::ggplotGrob(ggtree::ggtree(hmm$emission$invariants$tree)),
                                               t = 5, l=1, b=6, r=1)
                      dev.off()
                  } else if (plot.type=="dot") {
                      ## Add data.table to put line at dat.min
                      dat.line=data.table(x1=min(dat$loci),x2=max(dat$loci),y=dat.min,.id=factor(names(hmm$emission$data[[chain]])[-1],levels=rev(names(hmm$emission$data[[chain]])[-1])))
                      ## Convert species names to factors so they will plot in correct order
                      dat[,.id:=factor(as.character(.id),levels=rev(names(hmm$emission$data[[chain]])))]
                      g.dat=ggplot2::ggplot()+
                          ggplot2::geom_abline(data = dat.line, ggplot2::aes(intercept=y,slope=0), color="black", size=0.1)+
                          ggplot2::geom_point(data=dat,ggplot2::aes(x=loci,y=value,color=.id,fill=.id),inherit.aes=FALSE,size=0.5)+
                          ggplot2::facet_grid(.id~.)+
                          ggplot2::ylab("Species")+
                          ggplot2::scale_x_continuous(name=chrom,breaks=loci.breaks,labels=x.labs,limits=c(chrom.lim[1],chrom.lim[2]))+
                          cowplot::theme_cowplot()+
                          ggplot2::theme(axis.title.y=ggplot2::element_blank(),strip.text.y = ggplot2::element_blank())+
                          ggplot2::scale_y_continuous(name=ggplot2::element_blank(),breaks=c(dat.min,dat.max),labels=as.character(c(dat.min,dat.max)),limits=c(dat.min,dat.max))
                     ## Add calls overlay
                      if(nrow(spec)>1){
                          g.dat <- g.dat +
                              ggplot2::geom_rect(data=spec,ggplot2::aes(xmin=start,xmax=end,ymin=dat.min,ymax=dat.max),color="red",fill=NA)
                      }
                      pdf(NULL)
                      g.dat <- gtable::gtable_add_cols(ggplot2::ggplotGrob(g.dat), grid::unit(2,"cm"),0)
                   ##   g.dat <- gtable_add_grob(g.dat, ggplotGrob(ggtree::ggtree(hmm$emission$invariants$tree)),
                   ##                             t = 5, l=1, b=6, r=1)
                      dev.off()
                  }
                  plist[[2]]=g.dat
                  ## Plot scale factors
                  if(include.scale.factors){
                      sca=data.table::data.table(index=start:end,loci=chrom.loci[start:end],do.call(cbind,hmm$emission$invariants$scaleFactors[[chain]])[start:end,])
                      data.table::setnames(sca,colnames(sca)[c(-1,-2)],hmm$emission$invariants$tree$tip.label)
                      sca=data.table::melt(sca,id.vars=c('index','loci'))
                      sca[,variable:=factor(variable,levels=hmm$emission$invariants$tree$tip.label)]
                      sca[value>1,value:=1] ## hacky fix for earlylier off by one error, will fix later
                      ## Add deletion ranges
                      del=data.table::rbindlist(lapply(hmm$emission$invariants$deletionRanges[[chain]], function(x) data.table::data.table(s=x[,1]+1,e=x[,2]+1)),idcol=TRUE)
                      del[,.id:=factor(.id,levels=hmm$emission$invariants$tree$tip.label)]
                      del=del[with(del,IRanges::findOverlaps(IRanges::IRanges(s,e),IRanges::IRanges(start,end)))@from]
                      del[s<start,s:=start]
                      del[e>end,e:=end]
                      del[,s:=chrom.loci[s]]
                      del[,e:=chrom.loci[e]]
                      ## Convert del coordinates to 
                      g.scale=ggplot2::ggplot()+
                          ggplot2::geom_raster(data=sca,ggplot2::aes(x=loci,y=variable,fill=value),inherit.aes=FALSE)+
                          ggplot2::geom_rect(data=del,ggplot2::aes(xmin=s,xmax=e,ymin=as.numeric(.id)-0.5,ymax=as.numeric(.id)+0.5),inherit.aes=FALSE,color="red",fill="gray")+
                          viridis::scale_fill_viridis(limits=c(1/100,1),option="magma",na.value="white",direction=-1)+
                          ggplot2::scale_x_continuous(name=chrom,breaks=loci.breaks,limits=c(chrom.lim[1],chrom.lim[2]),labels=scales::comma)+
                          ggplot2::ylab("Species")+
                          ggplot2::guides(fill=ggplot2::guide_colorbar(title="Scale Factor"))+
                          cowplot::theme_cowplot()+
                          ggplot2::theme(axis.title.y=ggplot2::element_blank())
                  ## Add tree
                      pdf(NULL)
                      g.scale <- gtable::gtable_add_cols(ggplot2::ggplotGrob(g.scale), grid::unit(2,"cm"),0)
                      g.scale <- gtable::gtable_add_grob(g.scale, ggplot2::ggplotGrob(ggtree::ggtree(hmm$emission$invariants$tree)),
                                                         t = 5, l=1, b=6, r=1)
                      dev.off()
                      plist[[5]]=g.scale
                  }

              }
              if(!is.null(viterbi)){
                  vit=data.table::data.table(index=start:end,hmm$emission$invariants$treeStates[viterbi[[chain]][start:end],])
                  data.table::setnames(vit,colnames(vit)[-1],hmm$emission$invariants$tree$tip.label)
                  vit=data.table::melt(vit,id.vars='index')
                  vit[,loci:=chrom.loci[index]]
                  vit[,variable:=factor(variable,levels=hmm$emission$invariants$tree$tip.label)]
                  ## vit[,value:=as.factor(value)]
                  ## Add gridlines
                  g.path=ggplot2::ggplot()+
                      ggplot2::geom_raster(data=vit,ggplot2::aes(x=loci,y=variable,fill=as.factor(value)),inherit.aes=FALSE)+
                      ggplot2::geom_segment(data=xt,ggplot2::aes(x=x,y=y,xend=x,yend=yend),color="gray")+
                      ggplot2::geom_segment(data=yt,ggplot2::aes(x=x,y=y,xend=xend,yend=y),color="gray")+
                          ggplot2::scale_fill_manual(labels=c("Absent","Present"),values=c("white","#062F67"))+
                          ggplot2::ylab("Species")+
                          ggplot2::scale_x_continuous(name=chrom,breaks=loci.breaks,limits=c(chrom.lim[1],chrom.lim[2]),labels=scales::comma)+                    
                          ggplot2::guides(fill=ggplot2::guide_legend(title="Viterbi State"))+
                          cowplot::theme_cowplot()+
                          ggplot2::theme(legend.key = ggplot2::element_rect(colour = 'grey', fill = 'white', size = 1, linetype='solid'),axis.title.y=ggplot2::element_blank())
                  ## Add tree
                  pdf(NULL)
                  g.path <- gtable::gtable_add_cols(ggplot2::ggplotGrob(g.path), grid::unit(2,"cm"),0)
                  g.path <- gtable::gtable_add_grob(g.path, ggplot2::ggplotGrob(ggtree::ggtree(hmm$emission$invariants$tree)),
                                                    t = 5, l=1, b=6, r=1)
                  dev.off()
                  plist[[4]]=g.path

              }
              if(!is.null(marginal)){
                  mar=data.table::data.table(index=start:end,t(apply(marginal[[chain]][start:end,],1, function(x) x %*% hmm$emission$invariants$treeStates)))
                  data.table::setnames(mar,colnames(mar)[-1],hmm$emission$invariants$tree$tip.label)
                  mar=data.table::melt(mar,id.vars='index')
                  mar[,loci:=chrom.loci[index]]
                  mar[,variable:=factor(variable,levels=hmm$emission$invariants$tree$tip.label)]
                  mar[value>1,value:=1] ## deal with numeric precision placing values slightly above 1
                  g.mar=ggplot2::ggplot()+
                      ggplot2:::geom_raster(data=mar,ggplot2::aes(x=loci,y=variable,fill=value),inherit.aes=FALSE)+
                      ggplot2::scale_fill_gradient(limits=c(0,1),low="white", high="#062F67",na.value="white")+
                      ggplot2::geom_segment(data=xt,ggplot2::aes(x=x,y=y,xend=x,yend=yend),color="gray")+
                      ggplot2::geom_segment(data=yt,ggplot2::aes(x=x,y=y,xend=xend,yend=y),color="gray")+
                      ggplot2::ylab("Species")+
                      ggplot2::scale_x_continuous(name=chrom,breaks=loci.breaks,limits=c(chrom.lim[1],chrom.lim[2]),labels=scales::comma)+
                      ggplot2::guides(fill=ggplot2::guide_colorbar(title="Probability of Being\n Active"))+
                      cowplot::theme_cowplot()+
                      ggplot2::theme(axis.title.y=ggplot2::element_blank())
                   ## Add tree
                  pdf(NULL)
                  g.mar <- gtable::gtable_add_cols(ggplot2::ggplotGrob(g.mar), grid::unit(2,"cm"),0)
                  g.mar <- gtable::gtable_add_grob(g.mar, ggplot2::ggplotGrob(ggtree::ggtree(hmm$emission$invariants$tree)),
                                                     t = 5, l=1, b=6, r=1)
                  dev.off()
                  plist[[3]]=g.mar
              }
                         
              ## Remove unplotted plots
              keep=which(!unlist(lapply(plist,is.null)))
              plist=plist[keep]
              lab.pool=c("A","B","C","D")[1:sum(keep!=1)]
              labs=character(length(keep))
              labs[keep!=1]=lab.pool
              ## Arrange plots 
              g <- cowplot::plot_grid(plotlist=plist,ncol=1,labels = labs, align = "v",axis="lr",rel_heights = c(1,2,1,1,1)[keep])
              return(g)
          }
          
