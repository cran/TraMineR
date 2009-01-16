#######################################################
###Compare groups and returns most discriminating func
#######################################################

#seqecmpgroup.survival <- function(time, event, var, stest) {
#	library(survival)
#	test <- survdiff(Surv(time, event) ~ var, rho=stest)
#	pval <- 1-pchisq(test$chisq,length(test$n)-1)

#	pval
#}

#seqecmpgroup.chisq<-function(group,seqp){
#  chi<- chisq.test(group,seqp)
#  return(chi$statistic)
#}

seqecmpgroup<-function(subseq, group, method="chisq", pvalue.limit=NULL){
  seqecmpgroup.chisq<-function(index,group,seqmatrix,bonferroni,ntest){
    sp<-sum(seqmatrix[,index])
    if(sp>0&&sp<length(group)){
      suppressWarnings(chi<- chisq.test(group,seqmatrix[,index]))
      if(bonferroni) chi$p.value <- (1-(1-chi$p.value)^ntest)
      resid<-as.list(chi$residuals[,"1"])
      names(resid)<-paste("Resid",names(chi$residuals[,"1"]),sep=".")
      freq<-as.list(chi$observed[,"1"]/rowSums(chi$observed))
      names(freq)<-paste("Freq",names(chi$residuals[,"1"]),sep=".")
      return(data.frame(p.value=chi$p.value,statistic=chi$statistic,index=index,freq,resid))
      #return(data.frame(p.value=chi$p.value,statistic=chi$statistic,index=index))
    }
    vals<-numeric(length(levels(group)))
    vals2<-numeric(length(levels(group)))
    names(vals)<-paste("Resid",levels(group),sep=".")
    names(vals2)<-paste("Freq",levels(group),sep=".")
    vals2[]<-sp/length(group)
    return(data.frame(p.value=NA,statistic=NA,index=index,as.list(vals2),as.list(vals)))
#    return(data.frame(p.value=NA,statistic=NA,index=index))
  }
  if(is.null(pvalue.limit))pvalue.limit<- 2
  if(!is.subseqelist(subseq))stop("subseq should be a subseqelist")
  group<-factor(group)
  if(method=="bonferroni"){
    bonferroni<-TRUE
    method<-"chisq"
  }else{
    bonferroni<-FALSE
  }
  if(method=="chisq"){
    testfunc<-seqecmpgroup.chisq
    ntest<-length(subseq$subseq)
    seqmatrix<-seqeapplysub(subseq, method="presence")
    testfunc.arg<-list(group=group,bonferroni=bonferroni,ntest=ntest,seqmatrix=seqmatrix)
    decreasing<-FALSE
  }else{
    stop("This method is not (yet) implemented")
  }
	res<-data.frame()
	for (i in 1:length(subseq$subseq)) {
    testfunc.arg$index<-i
    stat<-do.call(testfunc,testfunc.arg)
    res<-rbind(res,stat)
	}
	subseqnum<-1:sum(res[,1]<=pvalue.limit)
  cres<-order(as.double(res[,1]),decreasing =decreasing)[subseqnum]#[!is.na(as.double(res$stat))]
  data<-data.frame(as.data.frame(subseq$data[cres,]),res[cres,])
  rownames(data)<-1:nrow(data)
  ret<-createsubseqelist(subseq$seqe,subseq$constraint,subseq$subseq[cres],data=data,type=method)
  ret$labels<-levels(group)
  class(ret)<-c("subseqelistchisq",class(ret))
	return(ret)
}



plot.subseqelistchisq<-function(x, ylim="uniform", rows=NA, cols=NA,
            residlevels=c(2,4), cpal=brewer.pal(1+2*length(residlevels),"RdBu"),legendcol=NULL,legend.cex=1,ptype="freq",...){

  if(!inherits(x,"subseqelistchisq"))stop("Subseq should be a result of seqecmpgroup")
  nplot<-length(x$labels)
#  print(ylim)

  #print(cpal)
  residbreaks<-c(-Inf,-sort(residlevels),sort(residlevels),Inf)
  lout <- TraMineR.setlayout(nplot, rows, cols, TRUE, "all")
  layout(lout$laymat, heights=lout$heights, widths=lout$widths)
 
  if(ptype=="resid"){
    baseIndex <- 4 +nplot
  }else{
     baseIndex<-4
  }
  if(ylim=="uniform"){
    ylim<-c(min(min(x$data[,(baseIndex+1):(baseIndex+nplot)]),0),max(x$data[,(baseIndex+1):(baseIndex+nplot)]))
  }
  for(i in 1:nplot){
    ccol<-as.character(cut(x$data[,4+nplot+i], breaks=residbreaks,labels=cpal))
    plot.subseqelist(x,freq=x$data[,baseIndex+i],col=ccol,main=x$labels[i],ylim=ylim,...)
  }
  savepar <- par(no.readonly = TRUE)
  par(mar = c(1, 1, 0.5, 1) + 0.1, xpd=FALSE)
  plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")
  title(main="Pearson residuals",cex=legend.cex)
  legncol<-length(c(paste("-",rev(residlevels)),"neutral",residlevels))
  if(is.null(legendcol)&&lout$legpos=="center"){
    legncol<-1
  }
  else if(!is.null(legendcol)&&legendcol)legncol<-1
  legend(lout$legpos,
			# inset=c(0,leg.inset),
			legend=c(paste("-",rev(residlevels)),"neutral",residlevels),
			fill=cpal,
			ncol=legncol,
			cex=legend.cex,
			bty="o"
      )
   par(savepar)
  
}