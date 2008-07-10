##########################################
###Plot a list of subsequence with support
##########################################

##seqefplot<-function(subseq,seqe, group=NULL,mfrow=c(2,3),ylim=c(0,1),...){
seqefplot<-function(subseq,seqe, group=NULL,mfrow=c(2,3),...){
  seqpmatrix<-seqeapplysub(subseq, seqe, method="presence")
  if(is.null(group)){
    support<-colSums(seqpmatrix)/nrow(seqpmatrix)
    osup<-order(support,decreasing=TRUE)
    sfreq<-support[osup]
    slegend<-colnames(seqpmatrix)[osup]
    barpos<-barplot(sfreq,names.arg=c(""),...)
    text(x=barpos,y=0.02,labels=slegend, srt=90, adj=c(0,0.5))
  }else{
    if(is.null(mfrow)){
      mfrow<-c(1,smax)
    }
    smax<-min(mfrow[1]*mfrow[2], (ncol(seqpmatrix)))
    par(mfrow=mfrow)
    if(is.factor(group)){
      lev<-levels(group)
    }else{
      lev<-names(table(group))
    }
    nlev<-length(lev)
    titles<-colnames(seqpmatrix)
    for(i in 1:smax){
      sfreq<-numeric(1)
      slegend<-character(1)
      nlev<-0
      for(j in 1:length(lev)){

        tab<-table(seqpmatrix[group==lev[j],i])
        tmpfreq<-as.double(tab[2])/(as.double(tab[1])+as.double(tab[2]))
        if(!is.na(tmpfreq)&tmpfreq>0){
          nlev<-nlev+1
          sfreq[nlev]<-tmpfreq
          slegend[nlev]<-paste(lev[j])
        }
      }
      if(nlev>0)  barplot(sfreq,names.arg=slegend,main=titles[i],...)
    }

  }

}